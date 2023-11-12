#ifndef MYCOEFF
#define MYCOEFF

#include <cmath>
#include "valInit.hpp"

using namespace std;

void coefficents_eval(double* Points, double* values, int size, int func, double* coeff)
{
    double* d = new double[size];
    for (int i = 1; i < size-1; i++)
	{
		d[i] = ((Points[i+1] - Points[i])*diffFun2(func, Points[i-1], Points[i]) - (Points[i] - Points[i-1])*diffFun2(func, Points[i], Points[i+1]))/(Points[i+1]-Points[i-1]);
	}
	d[0] = (3*diffFun2(func, Points[0], Points[1]) - d2Fun(func, Points[0])*(Points[1]-Points[0])/2 - d[1])/2;
	d[size-1] = (3*diffFun2(func, Points[size-2], Points[size-1]) + d2Fun(func, Points[size-1])*(Points[size-1]-Points[size-2])/2 - d[size-2])/2;

    for (int i = 0; i < size; i++)
    {
        coeff[4*i] = values[i];
        coeff[4*i + 1] = d[i];
        coeff[4*i + 2] = (3*diffFun2(func, Points[i], Points[i+1]) - 2*d[i] - d[i+1])/(Points[i+1] - Points[i]);
        coeff[4*i + 3] = (d[i] + d[i+1] - 2*diffFun2(func, Points[i], Points[i+1]))/pow(Points[i+1] - Points[i], 2);
    }
}

double approx_eval(double point, double* coeff, int intervalNum, double* points)
{
    double res = coeff[4*intervalNum] + coeff[4*intervalNum+1]*(point - points[intervalNum]) + coeff[4*intervalNum+1]*pow(point - points[intervalNum],2) + coeff[4*intervalNum+1]*pow(point - points[intervalNum],3);
    return res;
}


#endif
