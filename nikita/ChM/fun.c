#include "fun.h"

double Eigen_Value (const int m, const int N, const double h)
{
    return 4./h/h * sin(M_PI * m/ (2. *N - 1) * sin(M_PI * m /(2. *N - 1)));
}

int Eigen_Vector (double *y, const int m, const int N)
{
    for (int k = 0; k <= N; k++)
    {
        y[k] = sin(2. *M_PI * m * k/(2. *N - 1));
    }

    return 0;
}

double Scalar_Prod (const double *x1, const double *x2, const int N, const double h)
{
    double s = 0;

    for (int k = 0; k <= N - 1; k++)
    {
        s = s + x1[k]*x2[k]*h;
    }

    s = s - 3*x1[N]*x2[N]*h;

    return s;
}

double Measure(double *matrix, const double *y, const double Lambda, const int N, const double h)
{
    int k;
    double Measure;

    matrix[0] = 0;
    for (k = 1; k <= N-1; k++)
    {
        matrix[k] = (y[k+1] - 2*y[k] + y[k-1])/h/h + Lambda*y[k];
    }

    matrix[N] = (-3*y[N-1] + y[N-2]) /h/h + Lambda*y[N-1];

    Measure = Scalar_Prod(matrix, matrix, N, h);
    Measure = sqrt(Measure)/Lambda;
    return Measure;
}