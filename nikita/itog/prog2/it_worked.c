#include <math.h>
#include <stdio.h>
#include "task.h"
#include <stdlib.h>
#include <math.h>

#define MAX_OUTPUT_SIZE_2 5

double get_value(double* diag, double* side, int i, int j)
{
	if(i == j)
		return diag[i];
	if(i - j == 1)
		return side[j];
	if(i - j == -1)
		return side[i];
	return 0.;
}

void assign(double* diag, double* side, int i, int j, double value)
{
	if(i == j)
	{	
		diag[i] = value;
		return;
	}
	if (i - j == 1)
	{
		side[j] = value;
		return;
	}
	return;
}

double CalcNorm(int n, double* a)
{
	int i;
	int j;
	double tmp;
	double rezult;

	rezult = 0.0;
	for (i = 0; i < n; ++i)
	{
		tmp = 0.0;
		for (j = 0; j < n; ++j)
			tmp += fabs(a[i * n + j]);

		if (rezult < tmp)
			rezult = tmp;
	}

	return rezult;
}

void Shift(int n, double* a, int k, double s)
{
	int i;

	for (i = 0; i < k; ++i)
		a[i * n + i] -= s;
}

void SaveMatrix(int n, double* a)
{
	FILE* savefile;
	int i;
	int j;
	savefile = fopen("../checkpoint.txt", "w");
	if (!savefile)
	{
		printf("[X] CANNOT SAVE MATRIX\n");
	}
	fprintf(savefile, "%d\n", n);
	for (i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			fprintf(savefile, "%10.3g", a[i * n + j]);
		}
		fprintf(savefile, "\n");
	}
	fclose(savefile);
	printf("[V] SAVED!\n");
	return;
}


void PrintMatrix2(int n, double* a)
{
	int i;
	int j;
	int nPrint;

	nPrint = (n > MAX_OUTPUT_SIZE_2) ? MAX_OUTPUT_SIZE_2 : n;

	for (i = 0; i < nPrint; ++i)
	{
		for (j = 0; j < nPrint; ++j)
			printf("%10.3g ", a[i * n + j]);
		printf("\n");
	}
}


void left_multiplication(double* a, int n, int k)
{
	int j = 0;
	int l = 0;
	double tmp = 0.;

	for (j = k + 1; j < n; ++j)
	{
		tmp = 0.0;
		for (l = k + 1; l < n; ++l)
			tmp += a[l * n + k] * a[l * n + j];

		tmp *= 2.0;
		for (l = k + 1; l < n; ++l)
			a[l * n + j] -= tmp * a[l * n + k];
	}
}


void right_multiplication(double* a, int n, int k)
{
	double tmp;
	int j = 0;
	int l = 0;
	for (j = 0; j < n; ++j)
	{
		tmp = 0.0;
		for (l = k + 1; l < n; ++l)
			tmp += a[l * n + k] * a[j * n + l];

		tmp *= 2.0;
		for (l = k + 1; l < n; ++l)
			a[j * n + l] -= tmp * a[l * n + k];
	}
}


void cleaning(double* a, int n, int k, double norm)
{
	int j = 0;
	a[(k + 1) * n + k] = norm;
	a[k * n + (k + 1)] = norm;
	for (j = k + 2; j < n; ++j)
	{
		a[j * n + k] = 0.0;
		a[k * n + j] = 0.0;
	}
}


void TriDiag(double* a, int n)
{
	int j;
	int k;
	int i;
	double tmp1;
	double tmp2;
	double column_norm;
	double sk;
	double norm;
	double coefficient;

	int res = 0;
	for (i = 0; i < n - 2; ++i)
	{
		//=================================================
		// CREATE X
		//=================================================
		{
			tmp1 = 0.0;
			for (j = i + 2; j < n; ++j)
				tmp1 += a[j * n + i] * a[j * n + i];

			sk = tmp1;

			tmp2 = sqrt(a[(i + 1) * n + i] * a[(i + 1) * n + i] + sk);

			column_norm = tmp2;

			if (tmp2 < 1e-100)
			{
				a[(i + 1) * n + i] = 0.0;
				a[(i + 2) * n + i] = 0.0;

				continue;
			}

			if (tmp1 < 1e-100)
			{
				a[(i + 2) * n + i] = 0.0;

				continue;
			}

			// Create x(i) vector
			a[(i + 1) * n + i] -= column_norm;

			// Делим на норму a1 - ||a1||e1
			coefficient = 1.0 / sqrt(a[(i + 1) * n + i] * a[(i + 1) * n + i] + tmp1);
			for (j = i + 1; j < n; ++j)
				a[j * n + i] *= coefficient;	
		}
		//=================================================
		//=================================================

		left_multiplication(a, n, i);
		right_multiplication(a, n, i);
		cleaning(a, n, i, column_norm);
	}
}

void left_multiply_by_T_k_kpp(double* a, int n, int len, int k, double cos, double sin, double r)
{
	printf("L MULTIPLIED\n");
	int j;
	double x;
	double y;
	for (j = k + 1; j < len; ++j)
	{
		x = a[k * n + j];
		y = a[(k + 1) * n + j];

		a[k * n + j] = x * cos - y * sin;
		a[(k + 1) * n + j] = x * sin + y * cos;
	}
	a[k * n + k] = r;
	a[(k + 1) * n + k] = 0.0;
}
void right_multiply_by_T_k_kpp(double* a, int n, int len, int k, double cos, double sin, double r)
{	
	printf("R MULTIPLIED\n");
	int j;
	double x;
	double y;
	for (j = 0; j < len + 2; ++j)
	{
		if (abs(j - k) >= 2)
		{
			a[j * n + k] = 0;
		}
		else
		{
			x = a[j * n + k];
			y = a[j * n + k + 1];

			a[j * n + k] = x * cos - y * sin;
			a[j * n + k + 1] = x * sin + y * cos;
		}
	}
}


void QR_without_arrays(int n, double* a, int k)
{
	int i;
	int j;
	double x;
	double y;
	double r;

	double currcos;
	double currsin;
	double prevcos;
	double prevsin;

	for (i = 0; i < k - 1; ++i)
	{
		// COLLECT SIN COS
		{
		x = a[i * n + i];
		y = a[(i + 1) * n + i];

		// знаменатель для cos sin
		r = sqrt(x * x + y * y);

		if (r < 1e-100)
		{
			currcos = (a[i * n + i] > 0.0 ? 1.0 : -1.0);
			currsin = 0.0;
		}
		else
		{
			currcos = x / r;
			currsin = -y / r;
		}
		}
		left_multiply_by_T_k_kpp(a, n, k, i, currcos, currsin, r);
		if(i >= 1)
			right_multiply_by_T_k_kpp(a, n, k, i-1, prevcos, prevsin, r);
		prevcos = currcos;
		prevsin = currsin;
	}
	/*
	for (i = 0; i < k - 1; ++i)
	{
		right_multiply_by_T_k_kpp(a, n, k, i, cosPhi[i], sinPhi[i], r);
	}
	*/
	right_multiply_by_T_k_kpp(a, n, k, k-2, currcos, currsin, r);
}




void QR(int n, double* a, int k, double* cosPhi, double* sinPhi)
{
	int i;
	int j;
	double x;
	double y;
	double r;

	for (i = 0; i < k - 1; ++i)
	{
		// COLLECT SIN COS
		{
		x = a[i * n + i];
		y = a[(i + 1) * n + i];

		// знаменатель для cos sin
		r = sqrt(x * x + y * y);

		if (r < 1e-100)
		{
			cosPhi[i] = (a[i * n + i] > 0.0 ? 1.0 : -1.0);
			sinPhi[i] = 0.0;
		}
		else
		{
			cosPhi[i] = x / r;
			sinPhi[i] = -y / r;
		}
		}
		left_multiply_by_T_k_kpp(a, n, k, i, cosPhi[i], sinPhi[i], r);
	}
	for (i = 0; i < k - 1; ++i)
	{
		right_multiply_by_T_k_kpp(a, n, k, i, cosPhi[i], sinPhi[i], r);
	}
}

/*
void RQ(int n, double* a, int k, double* cosPhi, double* sinPhi)
{
	int i;
	int j;
	double x;
	double y;

	for (i = 0; i < k - 1; ++i)
	{
		for (j = 0; j < i + 2; ++j)
		{
			if (abs(j - i) >= 2)
			{
				a[j * n + i] = 0;
			}
			else
			{
				x = a[j * n + i];
				y = a[j * n + i + 1];

				a[j * n + i] = x * cosPhi[i] - y * sinPhi[i];
				a[j * n + i + 1] = x * sinPhi[i] + y * cosPhi[i];
			}
		}
	}
}
*/

void FindValues(int n, double* a, double* values, double eps, int* iterOut)
{
	double t;
	double shift;
	int i = 0;
	double D;
	int iter = 0;
	int k = 0;
	
 	//SaveMatrix(n, a);
	
	t = CalcNorm(n, a) * eps;
	double sin;
	double cos;

	for (k = n; k > 2; k--)
	{
		while(fabs(a[(k - 1) * n + k - 2]) > t)
		{
			shift = a[(k - 1) * n + k - 1];
			
			//QR(n, a, k, cosPhi, sinPhi);
			QR_without_arrays(n, a, k);

			Shift(n, a, k, shift);
			//Step(n, a, k);
			if (iter < 5)
			{
				printf("step=%d, k=%d\n", iter, k);
				PrintMatrix2(n, a);
			}
			Shift(n, a, k, -shift);
			++iter;
		}
	}
	// решение квадратного уравнения, полученного в явном виде
	// при поиске собственных значений оставшейся матрицы 2х2
	if (n > 1)
	{
		t = a[0 * n + 0] + a[1 * n + 1];
		D = a[0 * n + 0] * a[1 * n + 1] - a[0 * n + 1] * a[1 * n + 0];
		D = sqrt(t * t - 4.0 * D);

		a[0 * n + 0] = 0.5 * (t + D);
		a[1 * n + 1] = 0.5 * (t - D);
	}

	for (i = 0; i < n; ++i)
		values[i] = a[i * n + i];

	*iterOut = iter;
}
