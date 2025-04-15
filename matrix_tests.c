#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include "matrix.h"
#define EPSILON 1e-16


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

	bool res = dm_solve(2, 2, A, x, b);
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

	bool res = dm_solve(2, 2, A, x, b);
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

	bool res = dm_solve(2, 2, A, x, b);
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

	bool res = dm_solve(2, 2, A, x, b);
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

	bool res = dm_solve(3, 2, A, x, b);
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

	bool res = dm_solve(3, 2, A, x, b);
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

	bool res = dm_solve(2, 3, A, x, b);
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

int main() {
	test1();
	test2();
	test3();
	test4();
	test5();
	test6();
	test7();
	test8();
	test9();
	printf("All correct!\n");
}

