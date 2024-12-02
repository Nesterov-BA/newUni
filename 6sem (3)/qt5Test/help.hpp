#include <stdio.h>
#include <iostream>
#include <cmath>

double f0(double x, double y);
double f1(double x, double y);
double f2(double x, double y);
double f3(double x, double y);
double f4(double x, double y);
double f5(double x, double y);
double f6(double x, double y);
double f7(double x, double y);
double d2f0(double x, double y, int arg);
double d2f1(double x, double y, int arg);
double d2f2(double x, double y, int arg);
double d2f3(double x, double y, int arg);
double d2f4(double x, double y, int arg);
double d2f5(double x, double y, int arg);
double d2f6(double x, double y, int arg);
double d2f7(double x, double y, int arg);
double max_matr(double *M, int nx, int ny);
double min_matr(double *M, int nx, int ny);
double max4(double a, double b, double c, double d);
double min4(double a, double b, double c, double d);
double max(double a, double b);
double min(double a, double b);
void sqMatrMult(double* left, double* right, int size, int pos);
void transponse(double* matrix, int size);
void print_matrix(double* matrix, int nCols, int nRows);
