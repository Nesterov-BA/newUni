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

void coefficents_eval(double* Points, double* values, int size, double* coeff, double d2Start, double d2End)
{
    double* d = new double[size];
    for (int i = 1; i < size-1; i++)
    {
        d[i] = ((Points[i+1] - Points[i])*diffFun2(values[i-1], values[i], Points[i-1], Points[i]) + (Points[i] - Points[i-1])*diffFun2(values[i], values[i+1], Points[i], Points[i+1]))/(Points[i+1]-Points[i-1]);
    }
    d[0] = (3*diffFun2(values[0], values[1], Points[0], Points[1]) - d2Start*(Points[1]-Points[0])/2 - d[1])/2;
    d[size-1] = (3*diffFun2(values[size-2], values[size-1], Points[size-2], Points[size-1]) + d2End*(Points[size-1]-Points[size-2])/2 - d[size-2])/2;

    for (int i = 0; i < size-1; i++)
    {
        coeff[4*i] = values[i];
        coeff[4*i + 1] = d[i];
        printf("d[%d] = %f\n \n", i, d[i]);
        coeff[4*i + 2] = (3*diffFun2(values[i], values[i+1], Points[i], Points[i+1]) - 2*d[i] - d[i+1])/(Points[i+1] - Points[i]);
        coeff[4*i + 3] = (d[i] + d[i+1] - 2*diffFun2(values[i], values[i+1], Points[i], Points[i+1]))/pow(Points[i+1] - Points[i], 2);
    }
    delete[] d;
}

double approx_eval(double point, double* coeff, int intervalNum, double* points)
{
    double res = coeff[4*intervalNum] + coeff[4*intervalNum+1]*(point - points[intervalNum]) + coeff[4*intervalNum+2]*pow(point - points[intervalNum],2) + coeff[4*intervalNum+3]*pow(point - points[intervalNum],3);
    return res;
}

void evalArr(double* arr, int nInt, int n, double* coeff, double* points, double* intPoints)
{
    for(int i = 0; i < n-1; i++)
    {
        for(int j = 0; j < nInt; j++)
        {
            arr[i*nInt + j] = approx_eval(intPoints[i*nInt + j], coeff, i, points);
        }
    }
    arr[(n-1)*nInt] = approx_eval(intPoints[(n-1)*nInt], coeff, n-2, points);
}


#endif
