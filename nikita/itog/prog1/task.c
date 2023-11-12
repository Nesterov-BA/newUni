#include <math.h>

#include "task.h"
#include "matrix.h"

/*int InvertMatrix(int n, double* a, double* x, int* index)
{
	int i;
	int j;
	int k;
	int indMax;
	double tmp;
	double max;

	for (i = 0; i < n; ++i) // единичная матрица
		for (j = 0; j < n; ++j)
			x[i * n + j] = (double)(i == j);

	for (i = 0; i < n; ++i)
		index[i] = i;

	for (i = 0; i < n; ++i)
	{
		max = fabs(a[i * n + i]); // макс по диаг
		indMax = i;

		for (j = i + 1; j < n; ++j)
			if (max < fabs(a[i * n + j]))
			{
				max = fabs(a[i * n + j]);
				indMax = j;
			}

		k = index[i];
		index[i] = index[indMax];
		index[indMax] = k;

		for (j = 0; j < n; ++j)
		{
			tmp = a[j * n + i];
			a[j * n + i] = a[j * n + indMax];
			a[j * n + indMax] = tmp;
		}

		if (fabs(a[i * n + i]) < 1e-100)
			return -1;

		tmp = 1.0 / a[i * n + i];
		for (j = i; j < n; ++j)
			a[i * n + j] *= tmp;

		for (j = 0; j < n; ++j)
			x[i * n + j] *= tmp;

		for (j = 0; j < i; ++j)
		{
			tmp = a[j * n + i];
			for (k = i; k < n; ++k)
				a[j * n + k] -= a[i * n + k] * tmp;

			for (k = 0; k < n; ++k)
				x[j * n + k] -= x[i * n + k] * tmp;
		}

		for (j = i + 1; j < n; ++j)
		{
			tmp = a[j * n + i];
			for (k = i; k < n; ++k)
				a[j * n + k] -= a[i * n + k] * tmp;

			for (k = 0; k < n; ++k)
				x[j * n + k] -= x[i * n + k] * tmp;
		}
	}

	return 0;
}*/
int InvertMatrix(int n, double* a, double* x, int* index, double norm, int mode)
{
	int i;
	int j;
	int k;
	int indMax;
	double tmp;
	double max;
	double temp2;
	int k1 = 0;
	if(norm < 1e-50)
		return -1;
	for (i = 0; i < n; ++i)
		for (j = 0; j < n; ++j)
			x[i * n + j] = (double)(i == j);

	for (i = 0; i < n; ++i)
		index[i] = i;

	for (i = 0; i < n; ++i)
	{
		max = fabs(a[i * n + i]);
		indMax = i;

		for (j = i + 1; j < n; ++j)
			if (max < fabs(a[i * n + j]))
			{
				max = fabs(a[i * n + j]);
				indMax = j;
			}
		//printf("max = %f, max*1e+30 = %f\n", max, max*1e+16);
		temp2 = fabs(max/norm);
		while(temp2 < 1 && temp2 > 0)
		{
			temp2 = temp2*10;
			k1++;
		}
		

		if(k1 > 15 && mode == 1)
		{
			printf("i = %d, k = %d\n", i, k1);
			return -1;
		}
		k1 = 0;
		temp2 = 0;
		k = index[i];
		index[i] = index[indMax];
		index[indMax] = k;

		for (j = 0; j < n; ++j)
		{
			tmp = a[j * n + i];
			a[j * n + i] = a[j * n + indMax];
			a[j * n + indMax] = tmp;
		}

		if (fabs(a[i * n + i]) < norm*1e-17 && mode == 1)
			return -1;

		tmp = 1.0 / a[i * n + i];
		for (j = i; j < n; ++j)
			a[i * n + j] *= tmp;

		for (j = 0; j < n; ++j)
			x[i * n + j] *= tmp;

		for (j = i + 1; j < n; ++j)
		{
			tmp = a[j * n + i];
			for (k = i; k < n; ++k)
				a[j * n + k] -= a[i * n + k] * tmp;

			for (k = 0; k < n; ++k)
				x[j * n + k] -= x[i * n + k] * tmp;
		}
	}

	for (k = 0; k < n; ++k)
		for (i = n - 1; i >= 0; --i)
		{
			tmp = x[i * n + k];
			for (j = i + 1; j < n; ++j)
				tmp -= a[i * n + j] * x[j * n + k];
			x[i * n + k] = tmp;
		}

	for (i = 0; i < n; ++i)
		for (j = 0; j < n; ++j)
			a[index[i] * n + j] = x[i * n + j];

	for (i = 0; i < n; ++i)
		for (j = 0; j < n; ++j)
			x[i * n + j] = a[i * n + j];

	return 0;
}
