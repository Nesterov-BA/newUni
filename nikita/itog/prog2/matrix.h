#ifndef __MATRIX_H_INCLUDED__
#define __MATRIX_H_INCLUDED__

int InputMatrix(int n, int s, double* a, int inputMode, FILE* fin);

double get_norm(double* a, int n);

void PrintMatrix(int n, double* a);

void PrintVector(int n, double* x);

void optimize(double*a, double*diag, double*side, int n);

#endif /* not __MATRIX_H_INCLUDED__ */
