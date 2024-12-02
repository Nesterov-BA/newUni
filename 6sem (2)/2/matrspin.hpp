#ifndef MATRSPIN
#define MATRSPIN

#include <stdio.h>
#include <iostream>
#include "functions.hpp"
#include "printmatrix.hpp"
#include <cmath>

void spinmatr(double* buff, double a, double b, int i, int j, int n);
//double* spinmatrprod(double* vec, int n);
void triang(double* A, int size);

void spinmatr(double* buff, double a, double b, int k, int l, int n)
{
	id(buff, n);
	//printf("%d %d \n", k, l);
	if(a*a+b*b>0){
		buff[l+l*n] = a/sqrt(a*a+b*b);
		buff[l+k*n] = b/sqrt(a*a+b*b);
		buff[k+k*n] = buff[l+l*n];
		buff[k+l*n] = -buff[l+k*n];
	}
}


void triang(double* A, int size)
{
	double* buff = new double[size*size];
	for (int j = 0; j < size-1; ++j)
	{
		for (int i = j+2; i < size; ++i)
		{
			spinmatr(buff, A[j+(j+1)*size], A[j+i*size], j+1, i, size);
			sqmatr_mult(buff, A, size, 1);
			transponse(buff, size);
			sqmatr_mult(A, buff, size, 0);	
			//print_matrix(size, A);
		}
	}
}

#endif