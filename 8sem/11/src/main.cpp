#include "rungeKutta.hpp"
#include "metricFuncs.hpp"
#include "gaussian.hpp"
#include "gradient.hpp"
#include "graph.hpp"
// #include "normCalculation.hpp"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

double tolerance;
double alpha;
double finish;

double dp1(double p1, double p2, double x1, double x2, double t);
double dp2(double p1, double p2, double x1, double x2, double t);
double dx1(double p1, double p2, double x1, double x2, double t);
double dx2(double p1, double p2, double x1, double x2, double t);
std::vector<double> errorFunc(double p1, double p2);
std::vector<double> revErrorFunc(double p1, double p2);

double dp1(double p1, double p2, double x1, double x2, double t)
{
    p1 = p1;
    x1 = x1;
    x2 = x2;
    return p2/(1 + alpha*t*t);
}
double dp2(double p1, double p2, double x1, double x2, double t)
{
    p2 = p2;
    x1 = x1;
    x2 = x2;
    t = t;
    return -p1;
}
double dx1(double p1, double p2, double x1, double x2, double t)
{
    p1 = p1;
    p2 = p2;
    x1 = x1;
    t = t;
    return x2;
}
double dx2(double p1, double p2, double x1, double x2, double t)
{
    p1 = p1;
    x2 = x2;
    return p2 - x1/(1 + alpha*t*t);
}
function functions[] =
    {
        dp1,
        dp2,
        dx1,
        dx2
    };

std::vector<double> errorFunc(double p1, double p2)
{
    std::vector<double> result = error(p1, p2, functions);
    return result;
}
double errNorm(double p1, double p2)
{
    return l2Norm(errorFunc(p1, p2), 2);
}

int main(int argc, char* argv[])
{
    finish = M_PI/2;
    double start[4];
    double end[4];
    double integral = 0;
    tolerance = 1.e-11;
    if(argc > 1)
    {
         tolerance = atof(argv[1]);
    }
    alpha = 0;
    if(argc > 2)
        alpha = atof(argv[2]);

    printf("Tolerance: %.e, alpha: %lf\n", tolerance, alpha);
    start[0] = 2;
    start[1] = 0;
    start[2] = 0;
    start[3] = 1;
    solutionUpToTime(start, end, functions, finish);
    double start2[] = {2, 1};
    printf("End: %.e, %lf, %.e %.e\n", end[0], end[1], end[2], end[3] + (M_PI/2));
//    findMinimum(start2, errorFunc);
    probe(start2, errorFunc);
    return 0;
}