#include "rungeKutta.hpp"
#include "metricFuncs.hpp"
#include "gaussian.hpp"
#include "gradient.hpp"
#include "graph.hpp"
// #include "normCalculation.hpp"
#include <cstdio>
#include <iostream>
#include <vector>

double tolerance;
double finish;
double alpha;

double dp1(double p1, double p2, double x1, double x2);
double dp2(double p1, double p2, double x1, double x2);
double dx1(double p1, double p2, double x1, double x2);
double dx2(double p1, double p2, double x1, double x2);
std::vector<double> errorFunc(double p1, double p2);
std::vector<double> revErrorFunc(double p1, double p2);

double dp1(double p1, double p2, double x1, double x2)
{
    p1 = p1;
    p2 = p2;
    x1 = x1;
    return -48*(1 + alpha*x2*x2);
}
double dp2(double p1, double p2, double x1, double x2)
{
    p2 = p2;
    return 96*alpha*x1*x2*(1 + alpha*x2*x2)*(1 + alpha*x2*x2) - p1;
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
    x2 = x2;
    x1 = x1;
    return p2;
}
function functions[] =
    {
        dp1,
        dp2,
        dx1,
        dx2
    };
double revdp1(double p1, double p2, double x1, double x2)
{
    return p2*p2*x1;
}
double revdp2(double p1, double p2, double x1, double x2)
{
    return p1;
}
double revdx1(double p1, double p2, double x1, double x2)
{
    return -x2;
}
double revdx2(double p1, double p2, double x1, double x2)
{
    return -p2*x1*x1;
}
function revFunctions[] =
{
        revdp1,
        revdp2,
        revdx1,
        revdx2
};
std::vector<double> errorFunc(double p1, double p2)
{
    std::vector<double> result = error(p1, p2, functions, alpha);
    return result;
}
double errNorm(double p1, double p2)
{
    return l2Norm(errorFunc(p1, p2), 2);
}
std::vector<double> revErrorFunc(double p1, double p2)
{
    return revError(p1, p2, revFunctions, alpha);
}
double revErrorNorm(double p1, double p2)
{
    return l2Norm(revErrorFunc(p1, p2), 2);
}

int main(int argc, char* argv[])
{
    finish = 1;
    double start[4];
    double end[4];
    double initialP1[40];
    double initialP2[40];
    FILE *startFile = fopen("start.csv", "r");
    FILE *errors = fopen("errors.csv", "w");
    if(fscanf(startFile, "%lf,%lf,%lf,%lf", &start[0], &start[1], &start[2], &start[3]) != 4)
    {
        printf("Error reading start file");
        return 1;
    }
    printf("Start = %lf, %lf, %lf, %lf\n", start[0], start[1], start[2], start[3]);

    tolerance = 1.e-11;
    if(argc > 1)
    {
         tolerance = atof(argv[1]);
    }
    if(argc > 2)
        alpha = atof(argv[2]);

    for(int i = 0; i < 40; i++)
    {
        initialP1[i] = 0.1*i - 2;
        initialP2[i] = 0.1*i - 2;
    }
    solutionUpToTime(start, end, functions, finish);

    printf("\n");
    return 0;
}
