#include <stdio.h>
#include "matrix.h"

#define LINE_SIZE 3000

#define M 2000
#define K 137

int read_matrix(const char *filename, int m, int n, double *A) {
	FILE *file = fopen(filename, "r");
	if (file == NULL) {
		printf("File %s couldn't be opened.\n", filename);
		return -1;
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
				return -1;
			}
			lptr += offset;
		}
		i++;
	}
	fclose(file);
	if (i < m) {
		printf("Error: Not enough lines read.\n");
		return -1;
	}
	return 0;
}


int main() {
	double X[M][K];
	double y[M];
	read_matrix("/home/caganyanmaz/projects/cm-project/cm-sheet5-X.txt", M, K, X);
	read_matrix("cm-sheet5-y.txt", M, 1, y);
	for (int i = 0; i < M; i++) {
		printf("%f\n", y[i]);
	}

}
