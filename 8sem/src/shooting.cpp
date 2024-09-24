#include "math.h"
#include "rungeKutta.hpp"
#include "normCalculation.hpp"
void shooting(double* start, double alpha, function* functions, double finish)
{
    double eps = 1.e-3;
    double end[4];
    for(int i = 0; i < 4; i++)
    {
        end[i] = 0;
    }
    double errVector[2];
    errVector[0] = end[2];
    errVector[1] = end[3] + 1;
    while(l1Norm(errVector, 2) > eps)
        solutionUpToTime(start, end, functions, finish);
        start[0] +=eps;
        solutionUpToTime(start, end, functions, finish);
}

double error(double p1, double p2, function* functions)
{
    double start[4] = {p1, p1, 1, 1};
    double end[4];
    double finish = 2;
    solutionUpToTime(start, end, functions, finish);
    return fabs(end[2]) + fabs(end[3] + 1);
}