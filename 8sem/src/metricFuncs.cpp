#include "metricFuncs.hpp"

double l1Norm(double* vector, int size)
{
    double sum = 0;
    for(int i = 0; i < size; i++)
    {
        sum += fabs(vector[i]);
    }
    return sum;
}

double lInfNorm(double* vector, int size)
{
    double max = 0;
    for(int i = 0; i < size; i++)
    {
        if(fabs(vector[i]) > max)
            max = fabs(vector[i]);
    }
    return max;
}       

double l1Metric(double* vector1, double* vector2, int size)
{
    double* vector = new double[size];
    for(int i = 0; i < size; i++)
    {
        vector[i] = vector1[i] - vector2[i];
    }
    return l1Norm(vector, size);
}