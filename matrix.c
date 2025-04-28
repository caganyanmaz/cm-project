#include <assert.h>
#include <string.h>
#include <math.h>

#include "matrix.h"
#define EPSILON 1e-8

int min(int a, int b) {
	if (a < b) return a;
	return b;
}

void dm_add(int m, int n, const double *A, const double *B, double *C) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			C[i * n + j] = A[i * n + j] + B[i * n + j];
		}
	}
}

void dm_inline_add(int m, int n, double *A, const double *B) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			A[i * n + j] += B[i * n + j];
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

void dm_transpose(int m, int n, const double *A, double *B) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			B[j * m + i] = A[i * n + j];
		}
	}
}

double dm_dot_product(int n, const double *a, const double *b) {
	double res = 0;
	for (int i = 0; i < n; i++) {
		res += a[i] * b[i];
	}
	return res;
}

void dm_scalar_mul(int m, int n, double lambda, const double *A, double *B) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			B[i * n + j] = A[i * n + j] * lambda;
		}
	}
}

void dm_inline_scalar_mul(int m, int n, double lambda, double *A) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			A[i * n + j] *= lambda;
		}
	}
}

bool adm_solve(int m, int n, double *A, double *x, double *b) {
	int current_index = 0;
	for (int j = 0; j < min(m, n); j++) {
		int i = current_index;
		while (i < m && fabs(A[i * n + j]) <= EPSILON) i++;
		if (i == m)
			continue;
		adm_row_swap(m, n, i, current_index, A, b);
		adm_row_mul(m, n, current_index, 1 / A[current_index * n + j], A, b);
		for (i = current_index+1; i < m; i++) {
			adm_row_add(m, n, current_index, i, -A[i * n + j], A, b);
		}
		current_index++;
	}
	for (int i = current_index; i < m; i++) {
		if(fabs(b[i]) > EPSILON)
			return false;
	}
	// Unprocessed ones must be zero (not necessary, but I don't one (1923213123, -1232131.1238213813, 0.1231) as a result just because it's a possible solution)
	memset(x, 0, n * sizeof(double)); 

	for (int i = current_index-1; i >= 0; i--) {
		int pivot = 0;
		while (pivot < n && fabs(A[i * n + pivot] - 1) > EPSILON) pivot++;
		x[pivot] = b[i];
		for (int j = pivot + 1; j < n; j++) {
			x[pivot] -= x[j] * A[i * n + j];
		}
	}
	return true;
}

void adm_row_swap(int m, int n, int I, int J, double *A, double *b) {
	assert(0 <= I && I < m && 0 <= J && J < m);
	if (I == J)
		return;
	int tmp = b[I];
	b[I] = b[J];
	b[J] = tmp;
	for (int j = 0; j < n; j++) {
		tmp = A[I * n + j];
		A[I * n + j] = A[J * n + j];
		A[J * n + j] = tmp;
	}

}

void adm_row_mul(int m, int n, int I, double lambda, double *A, double *b) {
	assert(0 <= I && I < m);
	assert(fabs(lambda) > EPSILON);
	if (fabs(lambda - 1) <= EPSILON)
		return;
	b[I] *= lambda;
	for (int j = 0; j < n; j++) {
		A[I * n + j] *= lambda;
	}
}

void adm_row_add(int m, int n, int I, int J, double lambda, double *A, double *b) {
	assert(0 <= I && I < m && 0 <= J && J < m && I != J);
	if (fabs(lambda) < EPSILON)
		return;
	b[J] += b[I] * lambda;
	for (int j = 0; j < n; j++) {
		A[J * n + j] += A[I * n + j] * lambda;
	}
}

void dm_row_mul(int m, int n, int I, double lambda, double *A) {
	assert(0 <= I && I < m);
	if (fabs(lambda - 1) <= EPSILON)
		return;
	for (int j = 0; j < n; j++) {
		A[I * n + j] *= lambda;
	}
}

bool dm_check_equal(int m, int n, const double *a, const double *b, double epsilon) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (fabs(a[i * n + j] - b[i * n + j]) >= epsilon)
				return false;
		}
	}
	return true;
}


