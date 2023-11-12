#ifndef MYCOEFF
#define MYCOEFF

#include <cmath>
#include "help.hpp"

using namespace std;

void points(double* c, int n, double a, double b)
{
    double step = (b-a)/(n-1);
    for(int i = 0; i < n; i++)
        c[i] = a + step*i;
}

void Fill_F(double* F, int n, double* points, double(*f)(double))
{
    for(int i = 0; i < n; i++)
    {
        F[i] = f(points[i]);
    }
}

void thomas_algorithm(double* a,
                      double* b,
                      double* c,
                      double* d,
                      double* f, int N) {

    double* c_star = new double[N-1];
    double* d_star = new double[N];
    double m;

    c_star[0] = c[0]/b[0];
    d_star[0] = d[0]/b[0];

    for(int i = 1; i < N; i++)
    {
        m = b[i] - a[i-1]*c_star[i-1];
        if(i < N-1)
            c_star[i] = c[i]/m;
        d_star[i] = (d[i] - a[i-1]*d_star[i-1])/m;
    }

    f[N-1] = d_star[N-1];
    for(int i = N-2; i > -1; i--)
        f[i] = d_star[i] - c_star[i]*f[i+1];
    delete[] c_star;
    delete[] d_star;
}

void coefficents_eval(double* values, int size, double* coeff, double step, double d2Start, double d2End)
{
    double* v = new double[size+1];
    double* a = new double[size];
    double* b = new double[size+1];
    double* c = new double[size];
    double* d = new double[size+1];

    b[0] = 4;
    c[0] = b[0];
    a[size-1] =b[0];
    b[size] = b[0];
    d[0] = values[0]*8 + step*step*d2Start;
    d[size] = values[size-1]*8 + step*step*d2End;
    for(int row = 1; row < size; row++)
    {
        a[row-1] =  1;
        b[row] = 6;
        c[row] = 1;
        d[row] = 4*(values[row-1] + values[row]);
        v[row] = 0;
    }

    thomas_algorithm(a, b, c, d, v, size+1);

    for (int i = 0; i < size; i++)
    {
        coeff[3*i] = v[i];

        coeff[3*i + 1] = (2*values[i] - v[i+1] - v[i])/(step*step);

        coeff[3*i + 2] = 2*(v[i+1] + v[i] - 2*values[i])/(step*step);
    }
}

double approx_eval(double point, double* coeff, int intervalNum, double* points)
{
    double res = coeff[3*intervalNum] + coeff[3*intervalNum+1]*(point - points[intervalNum]) + coeff[3*intervalNum+2]*pow(point - points[intervalNum],2);
    return res;
}



void evalArr(double* arr, int nInt, int n, double* coeff, double* points, double* intPoints)
{
    for(int i = 0; i < n-1; i++)
    {
        for(int j = 0; j < nInt; j++)
        {
            if(intPoints[i*nInt + j] < (points[i] + points[i+1])/2)
                arr[i*nInt + j] = approx_eval(intPoints[i*nInt + j], coeff, i, points);
            else
                arr[i*nInt + j] = approx_eval(intPoints[i*nInt + j], coeff, i+1, points);
        }
    }
    arr[(n-1)*nInt] = approx_eval(intPoints[(n-1)*nInt], coeff, n-1, points);
}


#endif
