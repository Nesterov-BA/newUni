#ifndef _FUN_H_
#define _FUN_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double Eigen_Value (const int m, const int N, const double h);
int Eigen_Vector (double *y, const int m, const int N);
double Scalar_Prod (const double *x1, const double *x2, const int N, const double h);
double Measure(double *matrix, const double *y, const double Lambda, const int N, const double h);

#endif