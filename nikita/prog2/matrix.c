#include <stdio.h>
#include <math.h>

#include "matrix.h"

#define MAX_OUTPUT_SIZE 5

static double f(int n, int s, int i, int j)
{
    if (s == 1)
    {
        return 1.0 / (i + j + 1.0);
    }
    if (s == 2)
    {
        return fabs((double)i - (double)j);
    }
    if (s == 3)
    {
        return (double)n - fmax((double)i, (double)j)+1;
    }
    if (s == 4)
    {
        return fmax((double)i, (double)j);
    }
	if (s == 5)
	{
		return 1/((double)i + (double)j + 1);
	}
    return 0;
    
}



void optimize(double*a, double*diag, double*side, int n)
{
	int i;
	for(i = 0; i < n; i++)
	{
		diag[i] = a[i * n + i];
		if (i < n - 1)
			side[i] = a[(i + 1)*n + i];
	}
}

double get_norm(double* a, int n)
{
	double out = 0.;
	for (int i = 0; i < n; i++)
	{
		int j;
		double s = 0.;
		double *ptr = a + i;
		for (j = 0; j < i; j++)
		{
			s += fabs (*ptr);
			ptr += n - j - 1;
		}
		for (; j < n; j++)
		s += fabs (ptr[j - i]);
		if (s > out)
		out = s;
	}
	return out;
}

int InputMatrix(int n, int s, double* a, int inputMode, FILE* fin)
{
	int i;
	int j;

	if (inputMode == 1)
	{
		for (i = 0; i < n; ++i)
			for (j = 0; j < n; ++j)
				if (fscanf(fin, "%lf", &a[i * n + j]) != 1)
					return -1;
	}
	else
	{
		for (i = 0; i < n; ++i)
			for (j = 0; j < n; ++j)
                                a[i * n + j] = f(n, s, i, j);
	}

	return 0;
}

void PrintMatrix(int n, double* a)
{
	int i;
	int j;
	int nPrint;

	nPrint = (n > MAX_OUTPUT_SIZE) ? MAX_OUTPUT_SIZE : n;

	for (i = 0; i < nPrint; ++i)
	{
		for (j = 0; j < nPrint; ++j)
			printf("%10.5g ", a[i * n + j]);
		printf("\n");
	}
}

void PrintVector(int n, double* x)
{
	int i;
	int nPrint;

	nPrint = (n > MAX_OUTPUT_SIZE) ? MAX_OUTPUT_SIZE : n;

	for (i = 0; i < nPrint; ++i)
		printf("%5.5g ", x[i]);
	printf("\n");
}
