#ifndef _MATRIX_H_
#define _MATRIX_H_
#include <stdbool.h>


void dm_add(int m, int n, const double *A, const double *B, double *C);
void dm_inline_add(int m, int n, double *A, const double *B);
void dm_mul(int m, int k, int n, const double *A, const double *B, double *C);

void dm_transpose(int m, int n, const double *A, double *B);
double dm_dot_product(int n, const double *a, const double *b);

void dm_scalar_mul(int m, int n, double lambda, const double *A, double *B);
void dm_inline_scalar_mul(int m, int n, double lambda, double *A);

bool adm_solve(int m, int n, double *A, double *x, double *b);

void adm_row_swap(int m, int n, int I, int J, double *A, double *b);
void adm_row_mul(int m, int n, int I, double lambda, double *A, double *b); 
void adm_row_add(int m, int n, int I, int J, double lambda, double *A, double *b);

void dm_row_mul(int m, int n, int I, double lambda, double *A);


#endif
