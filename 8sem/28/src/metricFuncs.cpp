#include "metricFuncs.hpp"
#include "iostream"
#include <cstdio>
#include <vector>

double l1Norm(double* vector, int size)
{
    double sum = 0;
    for(int i = 0; i < size; i++)
    {
        sum += fabs(vector[i]);
    }
    return sum;
}
double l2Norm(std::vector<double> vector, int size)
{
    double norm = 0;
    for (int i = 0; i<size; i++)
        norm += vector[i]*vector[i];
    return sqrt(norm);
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

void normalize(std::vector<double> vector, int size)
{
    for(int i = 0; i < size; i++)
        vector[i] /= l2Norm(vector, size);
}

void printMatrix(double** matrix, int size)
{
    for(int i = 0; i < size; i++)
    {
        printf("|");
        for (int j = 0; j < size-1; j++)
        {
            printf("%lf ", matrix[i][j]);
        }
        printf("%lf|\n", matrix[i][size-1]);
    }
}

void printMatrix(double** matrix, int size, double* coordinates)
{
    printf("Matrix at coordinates %lf, %lf:\n", coordinates[0], coordinates[1]);
    printMatrix(matrix, size);
}
