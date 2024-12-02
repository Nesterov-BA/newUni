#include "math.h"
#include "metricFuncs.hpp"
#include "rungeKutta.hpp"
#include "gaussian.hpp"
#include <cstdio>
#include <vector>

void jacobiMatrix(double p1, double p2, std::vector<double> (*func)(double, double), double** matrix)
{
    double eps = 1.e-7;
    std::vector<double> res(2);
    std::vector<double> resX(2);
    std::vector<double> resY(2);
    res = func(p1, p2);
    resX = func(p1 + eps, p2);
    resY = func(p1, p2 + eps);
    resX[0] -= res[0];
    resX[0] /= eps;
    matrix[0][0] = resX[0];
    resX[1] -= res[1];
    resX[1] /= eps;
    matrix[0][1] = resX[1];
    resY[0] -= res[0];
    resY[0] /= eps;
    matrix[1][0] = resY[0];
    resY[1] -= res[1];
    resY[1] /= eps;
    matrix[1][1] = resY[1];
}

void probe(double* start, std::vector<double> (*function)(double, double))
{
    double tempStart[2];
    double minStart[2];
    int count = 0;
    double step = 2;
    double startsX[10];
    double startsY[10];
    double minNorm = 1;
    std::vector<double> res(2);
    for(int i = 0; i < 10; i++)
    {
        startsX[i] = start[0] + i*step;
    }
    for (int j = 0; j < 10; j ++)
    {
        startsY[j] = start[0] + j*step;
    }
    for(int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j ++)
        {
            tempStart[0] = startsX[i];
            tempStart[1] = startsY[j];
            printf("Starting from: %lf, %lf\n", tempStart[0], tempStart[1]);
            if(findMinimum(tempStart, function) == 0)
            {
                printf("Found minimum!");
                start[0] = startsX[i];
                start[1] = startsY[j];
                break;
            }
            res = function(tempStart[0], tempStart[1]);
            if(l2Norm(res, 2) < minNorm)
            {
                minNorm = l2Norm(res, 2);
                minStart[0] = tempStart[0];
                minStart[1] = tempStart[1];
            }
            printf("End at: %lf, %lf\n", tempStart[0], tempStart[1]);
            printf("\n");
        }
    }
    printf("Minimum: %lf", minNorm);
    printf("Minimum a: %lf, %lf\n", tempStart[0], tempStart[1]);
}

int findMinimum(double* start, std::vector<double> (*function)(double, double))
{
    double correction[2];
    double** jacobi = new double*[2];
    jacobi[0] = new double[2];
    jacobi[1] = new double[2];
    double temp[2];
    double coefficent = -1;
    std::vector<double> errorVector(2);
    double currNorm;
    int count = 0;
    while(true)
    {
        count++;
        coefficent = -1;
        errorVector = function(start[0], start[1]);
        correction[0] = errorVector[0];
        correction[1] = errorVector[1];
        currNorm = l2Norm(errorVector, 2);
        jacobiMatrix(start[0], start[1], function, jacobi);
        gauss(jacobi, correction);
        std::vector<double> res = function(start[0] + correction[0]*coefficent, start[1] + correction[1]*coefficent);
        while(l2Norm(function(start[0] + correction[0]*coefficent, start[1] + correction[1]*coefficent), 2) > currNorm)
        {
            coefficent /= 2;
            res = function(start[0] + correction[0]*coefficent, start[1] + correction[1]*coefficent);
            if(fabs(coefficent) < 1.e-7)
            {
                printf("Small Coefficent\n");
                printf("Res = (%lf, %lf)\n", errorVector[0], errorVector[1]);
                return -1;
            }
        }

        start[0] += correction[0]*coefficent;
        start[1] += correction[1]*coefficent;
        if(l2Norm(function(start[0], start[1]), 2) < 1.e-4)
        {
            printf("Small Error!");
            return 0;
            break;
        }
        if(count > 15)
        {
            printf("Big Count!\n");
            printf("Res = (%lf, %lf)\n", res[0], res[1]);
            return -1;
            break;
        }
    }
}

std::vector<double> error(double p1, double p2, function* functions, double alpha)
{
    double start[4] = {p1, p2, 1, alpha};
    double end[4];
    std::vector<double> res(2);
    double finish = 2;
    solutionUpToTime(start, end, functions, finish);
    res[0] = end[2];
    res[1] = end[3] + alpha;
    return res;
}
std::vector<double> revError(double p1, double p2, function* functions, double alpha)
{
    double start[4] = {p1, p2, 0, -alpha};
    double end[4];
    std::vector<double> res(2);
    double finish = 2;
    solutionUpToTime(start, end, functions, finish);
    res[0] = end[2] - 1;
    res[1] = end[3] - alpha;
    return res;
}
