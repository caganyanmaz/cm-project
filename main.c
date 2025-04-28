#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "matrix.h"

#define LINE_SIZE 3000

#define M 2000
#define K 137

// Max iter. count
#define N 1000
#define TOL 1e-4
#define EPSILON 1e-16


const double p = 0.5;
const double p_inv = 2;
const double delta = 1e-4;

double normal_random(double mean, double std) {
	double u = ((double)rand()) / RAND_MAX;
	double v = ((double)rand()) / RAND_MAX;
	double x = sqrt(-2 * log(u)) * cos(2 * M_PI * v) * std;
	return x;
}

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


double calculate_val(const double *X, const double *w, const double *y) {
	double res = 0;
	for (int i = 0; i < M; i++) {
		res += exp(-y[i] * dm_dot_product(K, w, &X[i * K]));
	}
	return res;
}

double tmpv[K];

void calculate_derivative(const double *X, const double *w, const double *y, double *res) {
	for (int i = 0; i < M; i++) {
		double lambda = -y[i] * exp(-y[i] * dm_dot_product(K, w, &X[i * K]));
		dm_scalar_mul(K, 1, lambda, &X[i * K], tmpv);
		dm_inline_add(K, 1, res, tmpv);
	}
}

double tmpm[M * K];
void calculate_hessian(const double *X, const double *X_t, const double *w, const double *y, double *H) {
	// Copying X into temporary matrix
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < K; j++) {
			tmpm[i * K + j] = X[i * K + j];
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
	double X[M * K], X_t[K * M], y[M], H[K * K], diff[K], step[K], w[K], tmp[K], tmpA[M * K];
	if(!read_matrix("/home/caganyanmaz/projects/cm-project/cm-sheet5-X.txt", M, K, X))
		return 0;
	if(!read_matrix("cm-sheet5-y.txt", M, 1, y))
		return 0;
	dm_transpose(M, K, X, X_t);

	// Setting up initial value (0 default)
	const char *weights_filename = "cm-sheet5-w.txt";
	if (access(weights_filename, F_OK) == 0) {
		printf("Loading pre-existing weights file...\n");
		read_matrix(weights_filename, K, 1, w);
		
	} else {
		printf("Initializing new weights...\n");
		for (int i = 0; i < M; i++) {
			w[i] = normal_random(0, 0.1);
		}
	}
	
	
	// Using Newton's Method
	for (int i = 0; i < N; i++) {
		double val = calculate_val(X, w, y);
		// I've shown the infimum is 0, so the optimal solution should get as close to 0 as possible (bacward error)
		printf("Error: %f\n", val);
		if (val < TOL) {
			printf("Optimal solution reached!\n");
			break;
		}
		// Finding the derivative, Hessian and the step at that point
		calculate_derivative(X, w, y, diff);
		calculate_hessian(X, X_t, w, y, H);
		memcpy(tmp, diff, sizeof(diff)); // Saving diff for further calculations
		assert(dm_check_equal(K, 1, tmp, diff, TOL));

		bool solved = adm_solve(K, K, H, step, tmp); // H and tmp are thrashed

		
		if (!solved) {
			printf("Couldn' solve the system of equations\n");
			break;
		}

		dm_inline_scalar_mul(K, 1, -1, step);
		// Applying Armijo rule until we get a good step (Armijo rule always works as we don't have stationary points and the function is convex)

		double prod = dm_dot_product(K, step, diff);
		assert(prod < 0);
		double alpha = 1;
		double new_val;
		do {
			dm_scalar_mul(K, 1, alpha, step, tmp);
			dm_inline_add(K, 1, tmp, w);
			new_val = calculate_val(X, tmp, y);
			alpha *= p;
		} while (new_val >= val + delta * alpha * p_inv * prod && alpha > EPSILON);
		if (alpha <= EPSILON)
		{
			printf("Armijo rule failed\n");
			break;
		}
		printf("%f %f\n", new_val, alpha);
		memcpy(w, tmp, sizeof(tmp));
		
		//dm_inline_add(K, 1, w, step);
	}
	int accurate_count = 0;
	for (int i = 0; i < M; i++) {
		if (y[i] * dm_dot_product(K, w, &X[i * K]) > 0)
		{
			accurate_count++;
		}
	}
	printf("%d\n", accurate_count);
	double accuracy = ((double)accurate_count) / M;
	printf("Accuracy: %f%\n", accuracy * 100);

	write_matrix(weights_filename, K, 1, w);


}
