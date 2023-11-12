#include <math.h>
#include <stdio.h>
#include "task.h"

#define MAX_OUTPUT_SIZE_2 5

static double OldCalcNorm(int n, double* a)
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

static void OldAlmost_triangular(int n, double* a)
{
	int i;
	int j;
	int k;
	double tmp1;
	double tmp2;
	double column_norm;
	double sk;
    // k <-> i
    //
	for (i = 0; i < n - 2; ++i)
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
		tmp1 = 1.0 / sqrt(a[(i + 1) * n + i] * a[(i + 1) * n + i] + tmp1);
		for (j = i + 1; j < n; ++j)
			a[j * n + i] *= tmp1;

		// Умножаем матрицу на U(x1) слева
		for (j = i + 1; j < n; ++j)
		{
			tmp1 = 0.0;
			for (k = i + 1; k < n; ++k)
				tmp1 += a[k * n + i] * a[k * n + j];

			tmp1 *= 2.0;
			for (k = i + 1; k < n; ++k)
				a[k * n + j] -= tmp1 * a[k * n + i];
		}

		// Умножаем матрицу на U(x1) справа
		for (j = 0; j < n; ++j)
		{
			tmp1 = 0.0;
			for (k = i + 1; k < n; ++k)
				tmp1 += a[k * n + i] * a[j * n + k];

			tmp1 *= 2.0;
			for (k = i + 1; k < n; ++k)
				a[j * n + k] -= tmp1 * a[k * n + i];
		}

		// Приводим столбец к почти треугольному виду
		a[(i + 1) * n + i] = tmp2;
		for (j = i + 2; j < n; ++j)
			a[j * n + i] = 0.0;
	}
}

void OldPrintMatrix_2(int n, double* a)
{
	int i;
	int j;
	int nPrint;

	nPrint = (n > MAX_OUTPUT_SIZE_2) ? MAX_OUTPUT_SIZE_2 : n;

	for (i = 0; i < nPrint; ++i)
	{
		printf("| ");
		for (j = 0; j < nPrint; ++j)
			printf("%10.3g ", a[i * n + j]);
		printf("|\n");
	}
}


static double OldGetShift(int n, double* a, int k)
{
	return a[(k - 1) * n + k - 1];
}

static void OldShift(int n, double* a, int k, double s)
{
	int i;

	for (i = 0; i < k; ++i)
		a[i * n + i] -= s;
}

static void OldQR(int n, double* a, int k, double* cosPhi, double* sinPhi)
{
	int i;
	int j;
	double x;
	double y;
	double r;

	for (i = 0; i < k - 1; ++i)
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

		// применяем матрицу поворота
		for (j = i + 1; j < k; ++j)
		{
			x = a[i * n + j];
			y = a[(i + 1) * n + j];

			a[i * n + j] = x * cosPhi[i] - y * sinPhi[i];
			a[(i + 1) * n + j] = x * sinPhi[i] + y * cosPhi[i];
		}

		a[i * n + i] = r;
		a[(i + 1) * n + i] = 0.0;
	}
}

static void OldRQ(int n, double* a, int k, double* cosPhi, double* sinPhi)
{
	int i;
	int j;
	double x;
	double y;

	for (i = 0; i < k - 1; ++i)
		for (j = 0; j < i + 2; ++j)
		{
			x = a[j * n + i];
			y = a[j * n + i + 1];

			a[j * n + i] = x * cosPhi[i] - y * sinPhi[i];
			a[j * n + i + 1] = x * sinPhi[i] + y * cosPhi[i];
		}
}

void OldFindValues(int n, double* a, double* values, double eps, double* cosPhi, double* sinPhi, int* iterOut)
{
	int i;
	int k;
	int iter;
	double t;
	double D;
	double shift;

	iter = 0;

	t = OldCalcNorm(n, a) * eps;
	
	// QR to almost Triangular
	//Almost_triangular(n, a);
	//printf("AFTER REF\n");
	//PrintMatrix_2(n, a);

	for (k = n; k > 2; --k)
		while (fabs(a[(k - 1) * n + k - 2]) > t)
		{
			// сдвиг
			shift = a[(k - 1) * n + k - 1];
			OldShift(n, a, k, shift);

			// QR разложение для почи треугольной матрицы
			// Q лежит в син кос, R лежит в А
			OldQR(n, a, k, cosPhi, sinPhi);

			// Умножаем матрицы R(cosphi, sinphi) и Q=А
			OldRQ(n, a, k, cosPhi, sinPhi);

			OldShift(n, a, k, -shift);

			++iter;
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
