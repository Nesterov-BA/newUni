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

double CalcNorm(int n, double* diag, double* side)
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
			tmp += fabs(get_value(diag, side, i, j));

		if (rezult < tmp)
			rezult = tmp;
	}

	return rezult;
}

//void Shift(int n, double* diag, int k, double s)
void Shift(double* diag, int k, double s)
{
	int i;
	for (i = 0; i < k; ++i)
		diag[i] -= s;
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

void PrintVector2(int n, double* x)
{
	int i;
	int nPrint;

	nPrint = (n > MAX_OUTPUT_SIZE_2) ? MAX_OUTPUT_SIZE_2 : n;

	for (i = 0; i < nPrint; ++i)
		printf("%5.5g ", x[i]);
	printf("\n");
}

void PrintTriDMatrix2(int n, double* diag, double* side)
{
	int i;
	int j;
	int nPrint;

	nPrint = (n > MAX_OUTPUT_SIZE_2) ? MAX_OUTPUT_SIZE_2 : n;

	for (i = 0; i < nPrint; ++i)
	{
		for (j = 0; j < nPrint; ++j)
			printf("%10.5g ", get_value(diag, side, i, j));
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
	//a[k * n + (k + 1)] = norm;
	for (j = k + 2; j < n; ++j)
	{
		a[j * n + k] = 0.0;
		a[k * n + j] = 0.0;
	}
}


void tridiag(double* a, int n)
{
	int j;
        //int k;
	int i;
	double tmp1;
	double tmp2;
	double column_norm;
	double sk;
        //double norm;
	double coefficient;

        //int res = 0;
	for (i = 0; i < n - 2; ++i)
	{
		
		// CREATE X
		
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

 
void left_multiply_by_T_k_kpp(double* diag, double* side, double* sideup, int k, double cos, double sin)
{
        //int j;
        //double x;
        //double y;
    
    
	double d1 = diag[k];
	double d2 = diag[k+1];
	double u1 = sideup[k];
	double u2 = sideup[k+1];
	double l1 = side[k];
//printf("k = %d\n", k);

    if(fabs(l1)<1e-6 && fabs(sin) < 1e-6)
    {
        l1=0.;
        sin = 0.;
    }
    
    if(fabs(d1)<1e-6 && fabs(cos) < 1e-6)
    {
        d1=0.;
        cos = 0.;
    }
      if(fabs(u1)<1e-6 && fabs(sin) < 1e-6)
    {
        u1=0.;
        sin = 0.;
    }
       if(fabs(d2)<1e-6 && fabs(cos) < 1e-6)
    {
        d2=0.;
        cos=0.;
    }
    
        if(fabs(d1)<1e-6 && fabs(sin) < 1e-6)
    {
        d1=0.;
        sin=0.;
    }
        if(fabs(l1)<1e-6 && fabs(cos) < 1e-6)
    {
        l1=0.;
        cos=0.;
    }
        if(fabs(u1)<1e-6 && fabs(cos) < 1e-6)
    {
        u1=0.;
        cos=0.;
    }
        if(fabs(d2)<1e-6 && fabs(sin) < 1e-6)
    {
        d2=0.;
        sin=0.;
    }
        if(fabs(u2)<1e-6 && fabs(cos) < 1e-6)
    {
        u2=0.;
        cos=0.;
    }
    
    
    //printf("d1 = %f, d2 = %f, u1 = %f, u2 = %f, l1 = %f, cos = %f, sin = %f\n ", d1, d2, u1, u2,l1, cos, sin);
	diag[k] = d1 * cos - l1 * sin;
	diag[k+1] = u1 * sin + d2 * cos;

	side[k] = d1 * sin + l1 * cos;
	sideup[k] = u1 * cos - d2 * sin;
	sideup[k+1] = u2 * cos;
    
}
//void right_multiply_by_T_k_kpp(double* diag, double* side, double* sideup, int n, int len, int k, double cos, double sin, double r)
void right_multiply_by_T_k_kpp(double* diag, double* side, double* sideup, int k, double cos, double sin)
{	
    

	double d1 = diag[k];
	double d2 = diag[k + 1];
	double u1 = sideup[k];
        //double u2 = sideup[k + 1];
        //double l1 = side[k];
	double l2 = side[k + 1];
	     if(fabs(l2)<1e-6 && fabs(cos) < 1e-6)
    {
        l2=0.;
        cos = 0.;
    }
    
    if(fabs(d1)<1e-6 && fabs(cos) < 1e-6)
    {
        d1=0.;
        cos = 0.;
    }
      if(fabs(u1)<1e-6 && fabs(sin) < 1e-6)
    {
        u1=0.;
        sin = 0.;
    }
       if(fabs(d2)<1e-6 && fabs(cos) < 1e-6)
    {
        d2=0.;
        cos=0.;
    }
    
        if(fabs(d1)<1e-6 && fabs(sin) < 1e-6)
    {
        d1=0.;
        sin=0.;
    }
       
        if(fabs(u1)<1e-6 && fabs(cos) < 1e-6)
    {
        u1=0.;
        cos=0.;
    }
        if(fabs(d2)<1e-6 && fabs(sin) < 1e-6)
    {
        d2=0.;
        sin=0.;
    }
    
	diag[k] = d1 * cos - u1 * sin;
	diag[k+1] = d2 * cos;
	side[k] = -d2 * sin;
	side[k+1] = l2 * cos;
	sideup[k] = d1 * sin + u1 * cos;
}


//void QR_without_arrays(int n, double* diag, double* side, double* sideup, int k)
void QR_without_arrays(double* diag, double* side, double* sideup, int k)
{
	int i;
        //int j;
	double x;
	double y;
	double r;

        double currcos = 1.;
        double currsin = 0.;
        double prevcos = 1.;
        double prevsin = 0.;

	// COPY VALUES
	for (i = 0; i < k; i++)
	{
		sideup[i] = side[i];
	}
	for (i = 0; i < k - 1; ++i)
	{
		// COLLECT SIN COS
		{
		x = get_value(diag, side, i, i);
		y = get_value(diag, side, i+1, i);

		// знаменатель для cos sin
		r = sqrt(x * x + y * y);

		if (r < 1e-100)
		{
			currcos = (get_value(diag, side, i, i) > 0.0 ? 1.0 : -1.0);
			currsin = 0.0;
		}
		else
		{
			currcos = x / r;
			currsin = -y / r;
		}
		}
                //left_multiply_by_T_k_kpp(diag, side, sideup, n, k, i, currcos, currsin, r);
                left_multiply_by_T_k_kpp(diag, side, sideup, i, currcos, currsin);
		if(i >= 1)
		{
                        right_multiply_by_T_k_kpp(diag, side, sideup, i-1, prevcos, prevsin);
		}
		prevcos = currcos;
		prevsin = currsin;
	}
        right_multiply_by_T_k_kpp(diag, side, sideup, k-2, currcos, currsin);
}

void FindValues(int n, double* diag, double* side, double* sideup, double* values, double eps, int* iterOut)
{
	double t;
	double shift;
	int i = 0;
	double D;
	int iter = 0;
	int k = 0;

 	//SaveMatrix(n, a);
        //double sin;
        //double cos;
	t = CalcNorm(n, diag, side) * eps;


	for (k = n; k > 2; k--)
	{
		while(fabs(get_value(diag, side, k-1, k-2)) > t)
		{
			shift = get_value(diag, side, k-1, k-1);
			
                        //Shift(n, diag, k, shift);
                        Shift(diag, k, shift);
                        QR_without_arrays(diag, side, sideup, k);

			//Shift(n, diag, k, -shift);
                        Shift(diag, k, -shift);
			++iter;
		}
	}
	// решение квадратного уравнения, полученного в явном виде
	// при поиске собственных значений оставшейся матрицы 2х2
	if (n > 1)
	{
		t = get_value(diag, side, 0, 0) + get_value(diag, side, 1, 1);
		D = get_value(diag, side, 0, 0) * get_value(diag, side, 1, 1) - get_value(diag, side, 0, 1) * get_value(diag, side, 1, 0);
		D = sqrt(t * t - 4.0 * D);

		assign(diag, side, 0, 0, 0.5 * (t + D));
		assign(diag, side, 1, 1, 0.5 * (t - D));
	}

	for (i = 0; i < n; ++i)
		values[i] = diag[i];
	
	*iterOut = iter;
}
