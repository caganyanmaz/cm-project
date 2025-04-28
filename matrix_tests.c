#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include "matrix.h"
#define EPSILON 1e-16


double rd(double lo, double hi) {
	return (((double) rand()) / RAND_MAX) * (hi - lo) + lo;
}

bool d_equal(double a, double b) {
	return fabs(a - b) < 4 * EPSILON;
}

void test1() {
	double A[] =
	{
		1, 0,
		0, 1
	};
	
	double x[] = 
	{
		0,
		0
	};

	double b[] = 
	{
		5,
		3,
	};

	bool res = adm_solve(2, 2, A, x, b);
	assert(res);
	assert(d_equal(x[0], 5));
	assert(d_equal(x[1], 3));
}

void test2() {
	double A[] =
	{
		1, 0,
		1, 1
	};
	
	double x[] = 
	{
		0,
		0
	};

	double b[] = 
	{
		5,
		3,
	};

	bool res = adm_solve(2, 2, A, x, b);
	assert(res);
	assert(d_equal(x[0], 5));
	assert(d_equal(x[1], -2));
}

void test3() {
	double A[] =
	{
		2, 0,
		1, 1
	};
	
	double x[] = 
	{
		0,
		0
	};

	double b[] = 
	{
		5,
		3,
	};

	bool res = adm_solve(2, 2, A, x, b);
	assert(res);
	assert(d_equal(x[0], 2.5));
	assert(d_equal(x[1], 0.5));
}

void test4() {
	double A[] =
	{
		0, 0,
		0, 1
	};
	
	double x[] = 
	{
		0,
		0
	};

	double b[] = 
	{
		0,
		1,
	};

	bool res = adm_solve(2, 2, A, x, b);
	assert(res);
	assert(d_equal(x[0], 0));
	assert(d_equal(x[1], 1));
}


void test5() {
	double A[] =
	{
		1, 0,
		1, 1,
		2, 2,
	};
	
	double x[] = 
	{
		0,
		0
	};

	double b[] = 
	{
		5,
		3,
		6,
	};

	bool res = adm_solve(3, 2, A, x, b);
	assert(res);
	assert(d_equal(x[0], 5));
	assert(d_equal(x[1], -2));
}


void test6() {
	double A[] =
	{
		1, 0,
		1, 1,
		2, 2,
	};
	
	double x[] = 
	{
		0,
		0
	};

	double b[] = 
	{
		5,
		3,
		5,
	};

	bool res = adm_solve(3, 2, A, x, b);
	assert(!res);
}

void test7() {
	double A[] =
	{
		1, 0, 3,
		1, 1, 5,
	};
	
	double x[] = 
	{
		0,
		0,
		0
	};

	double b[] = 
	{
		5,
		3,
	};

	bool res = adm_solve(2, 3, A, x, b);
	assert(res);
	assert(d_equal(x[0], 5));
	assert(d_equal(x[1], -2));
	assert(d_equal(x[2], 0));
}

void test8() {
	double A[] =
	{
		1, 1,
		1, 0
	};
	double B[] =
	{
		1, 0,
		0, 1
	};
	double C[] = 
	{
		0, 0,
		0, 0
	};
	double *_B = B;
	double *_C = C;
	for (int i = 0; i < 10; i++) {
		dm_mul(2, 2, 2, A, _B, _C);
		double *tmp = _B;
		_B = _C;
		_C = tmp;
	}
	double b[] = 
	{
		1,
		0
	};
	double c[] =
	{
		0,
		0
	};
	dm_mul(2, 2, 1, B, b, c);
	assert(d_equal(c[0], 89));
	assert(d_equal(c[1], 55));
}

void test9() {
	int m = 10;
	int n = 8;
	double A[10][8], B[10][8], C[10][8];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = i * i + j * j;
			B[i][j] = i + 123 * j;
		}
	}
	dm_add(10, 8, A, B, C);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			assert(d_equal(C[i][j], i * i + j * j + i + 123 * j));
		}
	}
}

void test10() {
	for (int i = 0; i < 1000; i++) 
	{
		printf("%d\n", i);
		double A[] = {
			rd(-10, 10), rd(-10, 10),
			rd(-10, 10), rd(-10, 10)
		};
		double x[] = {
			rd(-10, 10), rd(-10, 10)
		};
		double b[] = {0, 0};
		dm_mul(2, 2, 1, A, x, b);
		double tmpA[4];
		double tmpb[2];
		double tmp[2];
		memcpy(tmpA, A, sizeof(A));
		memcpy(tmpb, b, sizeof(b));
		adm_solve(2, 2, tmpA, x, tmpb);
		dm_mul(2, 2, 1, A, x, tmp);
		if (!dm_check_equal(2, 1, tmp, b, 1e-4)) {

			printf("A:\n%f %f\n%f %f\n", A[0], A[1], A[2], A[3]);
			printf("tmpA:\n%f %f\n%f %f\n", tmpA[0], tmpA[1], tmpA[2], tmpA[3]);
			printf("x:\n%f %f\n", x[0], x[1]);
			printf("b:\n%f %f\n", b[0], b[1]);
			printf("tmpb:\n%f %f\n", tmpb[0], tmpb[1]);
			assert(false);
		}
	

	}
}

int main() {
	srand(0);
	test1();
	test2();
	test3();
	test4();
	test5();
	test6();
	test7();
	test8();
	test9();
	test10();
	printf("All correct!\n");
}

