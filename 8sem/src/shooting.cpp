#include "math.h"
#include "rungeKutta.hpp"
#include "normCalculation.hpp"

double error(double p1, double p2, function* functions, double alpha)
{
    double start[4] = {p1, p1, alpha, 1};
    double end[4];
    double finish = 2;
    solutionUpToTime(start, end, functions, finish);
    return sqrt(end[2]*end[2] + (end[3]+alpha)*(end[3]+alpha));
}


void shooting(double* start, double alpha, function* functions, double finish)
{
    double eps = 1.e-7;
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