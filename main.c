#include <stdio.h>
#include <math.h>
#include "matrix.h"

#define LINE_SIZE 3000

#define M 2000
#define K 137

// Max iter. count
#define N 1000
#define TOL 1e-8

bool read_matrix(const char *filename, int m, int n, double *A) {
	FILE *file = fopen(filename, "r");
	if (file == NULL) {
		printf("File %s couldn't be opened.\n", filename);
		return false;
	}

	char line[LINE_SIZE];

	int i = 0;
	while (fgets(line, sizeof(line), file) != NULL && i < m) {
		char *lptr = line;
		for (int j = 0; j < n; j++) {
			int offset;
			if(sscanf(lptr, "%lf%n", &A[i * n + j], &offset) != 1)
			{
				printf("Error at: %d %d\n", i, j);
				return false;
			}
			lptr += offset;
		}
		i++;
	}
	fclose(file);
	if (i < m) {
		printf("Error: Not enough lines read.\n");
		return false;
	}
	return true;
}

bool write_matrix(const char *filename, int m, int n, double *A) {
	FILE *file = fopen(filename, "w");
	if (file == NULL) {
		printf("File %s couldn't be opened.\n", filename);
		return false;
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if(fprintf(file, "%.16lf ", A[i * n + j]) < 0) {
				printf("Error while writing number in file %s.\n", filename);
				return false;
			}
		}
		if(fprintf(file, "\n") < 0) {
			printf("Error while writing newline character in file %s.\n", filename);
			return false;
		}
	}
	return true;

}


double calculate_val(double *X, double *w, double *y) {
	double res = 0;
	for (int i = 0; i < M; i++) {
		res += exp(-y[i] * dm_dot_product(K, w, &X[i * K]));
	}
	return res;
}

double tmpv[K];

void calculate_derivative(double *X, double *w, double *y, double *res) {
	for (int i = 0; i < M; i++) {
		double lambda = -y[i] * exp(-y[i] * dm_dot_product(K, w, &X[i * K]));
		dm_scalar_mul(K, 1, lambda, &X[i * K], tmpv);
		dm_inline_add(K, 1, res, tmpv);
	}
}

double tmpm[M][K];
void calculate_hessian(double *X, double *X_t, double *w, double *y, double *H) {
	// Copying X into temporary matrix
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < K; j++) {
			tmpm[i][j] = X[i * K + j];
		}
	}

	// Multiplying each row according to the exponential
	for (int i = 0; i < M; i++) {
		double lambda = exp(-y[i] * dm_dot_product(K, w, &X[i * K]));
		dm_row_mul(M, K, i, lambda, tmpm);
	}
	dm_mul(K, M, K, X_t, tmpm, H);
}

int main() {
	double X[M][K], X_t[K][M], y[M], H[K][K], diff[K], step[K], w[K], bw[K];
	double best_val = 1e9;
	if(!read_matrix("/home/caganyanmaz/projects/cm-project/cm-sheet5-X.txt", M, K, X))
		return 0;
	if(!read_matrix("cm-sheet5-y.txt", M, 1, y))
		return 0;
	dm_transpose(M, K, X, X_t);

	// Setting up initial value (0 default)
	for (int i = 0; i < M; i++) {
		w[i] = 0;
	}
	
	// Using Newton's Method
	for (int i = 0; i < N; i++) {
		double val = calculate_val(X, w, y);
		if (val < best_val) {
			best_val = val;
			for (int j = 0; j < K; j++) {
				bw[j] = w[j];
			}
		}
		// I've shown the infimum is 0, so the optimal solution should get as close to 0 as possible (bacward error)
		printf("Error: %f\n", val);
		if (val < TOL) {
			printf("Optimal solution reached!\n");
			break;
		}
		calculate_derivative(X, w, y, diff);
		calculate_hessian(X, X_t, w, y, H);
		bool solved = adm_solve(K, K, H, step, diff);
		if (!solved) {
			printf("Couldn' solve the system of equations\n");
			break;
		}
		dm_inline_scalar_mul(K, 1, -1, step);
		dm_inline_add(K, 1, w, step);
	}
	printf("Best error rate: %f\n", best_val);
	int accurate_count = 0;
	for (int i = 0; i < M; i++) {
		if (y[i] * dm_dot_product(K, bw, &X[i]) > 0)
		{
			printf("Nice\n");
			accurate_count++;
		}
	}
	printf("%d\n", accurate_count);
	double accuracy = ((double)accurate_count) / M;
	printf("Accuracy: %f%\n", accuracy * 100);


	write_matrix("cm-sheet5-w.txt", K, 1, bw);


}
