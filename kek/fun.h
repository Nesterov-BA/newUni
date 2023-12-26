#ifndef _FUN_H_
#define _FUN_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double egenValue ( const int m, const int N, const double h);
int egenVector ( double *y, const int m, const int N);
double prod (const double *u, const double *v, const int N, const double h);
double measure( double *A, const double *y, const double lambda, const int N, const double h); 


#endif
