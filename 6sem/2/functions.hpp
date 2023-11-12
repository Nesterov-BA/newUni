#ifndef FUNC
#define FUNC
#include <stdio.h>
#include <iostream>
#include <cmath>
#include "printmatrix.hpp"

double dabs(double x);
double scalar_product(double* vector1, double* vector2, int n);
double vecmod(double* vec, int n);
void normalize(double* vec, int n);
void sqmatr_mult(double* A, double* B, int n, int pos);
void vectimesvec(double* buff, double* vec1, double* vec2, int n);
void id(double* A,int n);
void matreq(double* A, double* B, int n);
void transponse(double* A, int n);
double checkunit(double* A, int n);

double dabs (double x){
	if (x<0) {
		return -x;
	} else {
		return x;
	}
}

double scalar_product(double* vector1, double* vector2, int n)
{
	double buff = 0;
	for (int i = 0; i < n; ++i)
	{
		buff += vector1[i]*vector2[i];
	}
	return buff;
}

double vecmod(double* vec, int n)
{
	double mod = 0;
	mod = scalar_product(vec, vec, n);
	mod = sqrt(mod);
	return mod;
}

void normalize(double* vec, int n) 
{
	double a = vecmod(vec, n);
	if (a > 0)
	{
		for (int i = 0; i < n; ++i)
		{
			vec[i] /= a;
		}
	} //else
	//{
	//	printf("i break\n");
	//}
}

void sqmatr_mult(double* A, double* B, int n, int pos)
{
	
	double* buff = new double[n*n];
	double* row = new double [n];
	double* col = new double [n]; 
	for(int i = 0; i<n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			for (int m = 0; m < n; ++m)
			{
				row[m] = A[m+i*n];
				col[m] = B[j + m*n];
			}
			buff[j+i*n] = scalar_product(row, col, n);
			
		}
	}
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if(pos == 0) 
			{
				A[j+n*i] = buff[j+n*i];
			} else 
			{
				B[j+n*i] = buff[j+n*i];
			}
		}
	}
}


void vectimesvec(double* buff, double* vec1, double* vec2, int n)
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			buff[j+n*i] = vec1[i]*vec2[j];
		}
	}
}

void id(double* A, int n)
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if (i == j)
			{
				A[j+n*i] = 1;
			} 
			else 
			{
				A[j+n*i] = 0;
			}
		}
	}
}

void matreq(double* A, double* B, int n)
{
	for (int i = 0; i < n; ++i)
	{
		A[i] = B[i];
	}
}

void transponse(double* A, int n)
{
	double buff = 0;
	for (int i = 0; i < n-1; ++i)
	{
		for (int j = i+1; j < n; ++j)
		{
			buff = A[j+n*i];
			A[j+n*i] = A[i+n*j];
			A[i+n*j] = buff;
		}
	}
}

double checkunit(double* A, int n)
{
	double* buff = new double[n*n];
	double norm = 0;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			buff[j+n*i] = A[j+n*i];
		}
	}
	transponse(buff, n);
	sqmatr_mult(A, buff, n, 1);
	for (int i = 0; i < n; ++i)
	{
		buff[i+n*i] -= 1;
	}
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			norm += buff[i+n*j]*buff[i+n*j];
		}
	}
	return norm;
}


#endif