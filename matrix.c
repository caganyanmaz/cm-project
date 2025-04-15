
void dm_add(int m, int n, const double *A, const double *B, double *C) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			C[i * n + j] = A[i * n + j] + B[i * n + j];
		}
	}
}
void dm_mul(int m, int k, int n, const double *A, const double *B, double *C) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			C[i * n + j] = 0;
			for (int r = 0; r < k; r++) {
				C[i * n + j] += A[i * k + r] * B[r * n + j];
			}
		}
	}
}

void dm_row_swap(int m, int n, int I, int J, double *A, double *b);
void dm_row_mul(int m, int n, int I, double lambda, double *A, double *b); 
void dm_row_add(int m, int n, int I, int J, double lambda, double *A, double *b);

int dm_solve(int m, int n, double *A, double *b) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < 
	}
}

void dm_row_swap(int m, int n, int I, int J, double *A, double *b) {}
void dm_row_mul(int m, int n, int I, double lambda, double *A, double *b) {}
void dm_row_add(int m, int n, int I, int J, double lambda, double *A, double *b) {}


