#include "math.h"
#include "rungeKutta.hpp"

void shooting(double* start, double alpha, function* functions, double finish)
{
    double jacobiMatrix[2][2];
    double eps = 1.e-5;
    double end[4];
    for(int i = 0; i < 4; i++)
    {
        end[i] = 0;
    }
    double errVector[2];
    errVector[0] = end[2];
    errVector[1] = end[3] + 1;
    while(l1Norm(errVector, 2) > eps)
    {
        solutionUpToTime(start, end, functions, finish);
        start[0] +=eps;
        solutionUpToTime(start, end, functions, finish);
    }
}

double error(double p1, double p2, function* functions, double alpha)
{
    double start[4] = {p1, p2, alpha, 1};
    double end[4];
    double finish = 2;
    solutionUpToTime(start, end, functions, finish);
    return sqrt(end[2]*end[2] + (end[3]+alpha)*(end[3]+alpha));
}

void jacobiMatrix(double p1, double p2, function* functions, double** matrix)
{
   double eps = 1.e-5;
   double start[4] = {p1, p2, alpha, 1};
   double res1, res2;
   double end[4];

   solutionUpToTime(start, end, functions, finish);
   res1 = end[2];
   res2 = end[3];

   start[0] += eps;
   solutionUpToTime(start, end, functions, finish);
   matrix[0][0] = end[2] - res1;
   matrix[1][0] = end[3] - res2;

   start[0] -= eps;
   start[1] += eps;
   solutionUpToTime(start, end, functions, finish);
   matrix[0][1] = end[2] - res1;
   matrix[1][1] = end[3] - res2;

   matrix[0][0] /= eps;
   matrix[1][0] /= eps;
   matrix[0][1] /= eps;
   matrix[1][1] /= eps;
}
