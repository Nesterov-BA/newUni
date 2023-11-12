#ifndef __MATRIX_H_INCLUDED__
#define __MATRIX_H_INCLUDED__

#include "stdio.h"

double matrNorm(double* matrix, int size);

int InputMatrix(int n, int s, double* a, int inputMode, FILE* fin);

void PrintMatrix(int n, double* a);

double SolutionError(int n, double* a, double* x);

void PrintVector(int n, double* x);

void PrintNorm(int size, double* matrix);

#endif /* not __MATRIX_H_INCLUDED__ */
