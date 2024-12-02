#include "metricFuncs.hpp"
#include "graph.hpp"
#include "shooting.hpp"
#include "rungeKutta.hpp"
#include "testFuncs.hpp"
// #include "normCalculation.hpp"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

double tolerance;
double alpha;
double finish;

double dp1(double p1, double p2, double x1, double x2);
double dp2(double p1, double p2, double x1, double x2);
double dx1(double p1, double p2, double x1, double x2);
double dx2(double p1, double p2, double x1, double x2);
double errorFunc(double x2);

double dp1(double p1, double p2, double x1, double x2)
{
    p1 = p1;
    p2 = p2;
    x1 = x1;
    return -48/(1 + alpha*x2*x2);
}
double dp2(double p1, double p2, double x1, double x2)
{
    p2 = p2;
    return -96*alpha*x1*x2/((1 + alpha*x2*x2)*(1 + alpha*x2*x2)) - p1;
}
double dx1(double p1, double p2, double x1, double x2)
{
    p1 = p1;
    p2 = p2;
    x1 = x1;
    return x2;
}
double dx2(double p1, double p2, double x1, double x2)
{
    p1 = p1;
    x1 = x1;
    x2 = x2;
   // return -1;
    return p2;
}
function functions[] =
    {
        dp1,
        dp2,
        dx1,
        dx2
    };

vector<double >errorFunc(double p2, double x1)
{
    vector<double> result = error(p2, x1, functions);
    return result;
}

int main(int argc, char** argv)
{
    if(argc > 1)
        tolerance = atof(argv[1]);
    if(argc > 2)
        alpha = atof(argv[2]);
    double start[] = {0, -24, 10, 0};
    double start2[] = {-24,10};
    double* end = new double[4];
    finish = 1;
    solutionUpToTime(start, end, functions, finish);
    findMinimum(start2, errorFunc);
    printf("start: %lf, %lf\n", start2[0], start2[1]);
//    probe(start2, errorFunc);
    printf("End: %lf, %lf, %lf, %lf\n", end[0], end[1], end[2], end[3]);
/*
    printf("Testing function (x^2, y^2)\n");
    findMinimum(start2, x2y2);
    start2[0] = 5;
    start2[1] = 4;
    printf("Testing function (exp(x^2)-2, exp(y^2)-2)\n");
    findMinimum(start2, expofx2y2);
*/
    return 0;
}
