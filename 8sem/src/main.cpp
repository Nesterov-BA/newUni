#include "rungeKutta.hpp"
#include "gaussian.hpp"
#include "gradient.hpp"
#include "graph.hpp"
// #include "normCalculation.hpp"
#include <cstdio>
#include <iostream>

double dp1(double p1, double p2, double x1, double x2);
double dp2(double p1, double p2, double x1, double x2);
double dx1(double p1, double p2, double x1, double x2);
double dx2(double p1, double p2, double x1, double x2);
double errorFunc(double p1, double p2);

double dp1(double p1, double p2, double x1, double x2)
{
    return -p2*p2*x1;
}
double dp2(double p1, double p2, double x1, double x2)
{
    return -p1;
}
double dx1(double p1, double p2, double x1, double x2)
{
    return x2;
}
double dx2(double p1, double p2, double x1, double x2)
{
    return p2*x1*x1;
}
function functions[] =
    {
        dp1,
        dp2,
        dx1,
        dx2
    };
double errorFunc(double p1, double p2)
{
    return error(p1, p2, functions, alpha);
}

int main(int argc, char* argv[])
{
    finish = 2;
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

    solutionUpToTime(start, end, functions, finish);
    double startRev[4];
    double endRev[4] = {10,10,0,-alpha};
    solutionReverse(startRev, endRev, functions, finish);
    printf("Error vector:(%lf,%lf)\n", fabs(end[2]), fabs(end[3]+1));
    for(int i = 0; i < 40; i++)
    {
        initialP1[i] = 0.1*i - 2;
        initialP2[i] = 0.1*i - 2;
    }
    // double alpha = argc > 2 ? atof(argv[2]) : 1;
    // for(int i = 0; i < 40; i++)
    // {
    //      for(int j = 0; j < 40 ; j++)
    //      {
    //          fprintf(errors, "%lf,%lf,%lf\n", initialP1[i], initialP2[j], error(initialP1[i], initialP2[j], functions, alpha));
    //      }
    // }
    double p[2] = {20, -15};
    double grad[2];
   // gradient(errorFunc, p, grad);
    //printf("Gradient: %lf, %lf\n", grad[0], grad[1]);

//    newton2d(errorFunc, p);
    //double bounds[4] = {21, 23, -18, -15};
    //graph("errorGraph.csv", errorFunc, bounds);
    //gradient_descent(errorFunc, p) ;
    // // /////////////////////////////////
    // // // double vec[3] = {1, 2, 3};\ //
    // // // gaussian(vec);              //
    // // // printf("Gaussian:\n");      //
    // // // for(int i = 0; i < 3; i++)  //
    // // // {                           //
    // // //     printf("%lf ", vec[i]); //
    // // // }                           //
    /////////////////////////////////
    printf("\n");
   return 0;
}
