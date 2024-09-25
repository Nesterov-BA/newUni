#include "rungeKutta.hpp"
#include "shooting.hpp"
// #include "normCalculation.hpp"
#include <iostream>

double dp1(double p1, double p2, double x1, double x2);
double dp2(double p1, double p2, double x1, double x2);
double dx1(double p1, double p2, double x1, double x2);
double dx2(double p1, double p2, double x1, double x2);

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

int main(int argc, char* argv[]) 
{
   double finish = 2;
   double* start = new double[4];
   double* end = new double[4];
   double initialP1[40];
   double initialP2[40];
   FILE *startFile = fopen("start.csv", "r");
   FILE *errors = fopen("errors.csv", "w");
   if(fscanf(startFile, "%lf,%lf,%lf,%lf", &start[0], &start[1], &start[2], &start[3]) != 4)
   {
        printf("Error reading start file");
        return 1;
   }
   function functions[] = 
   {
        dp1,
        dp2,
        dx1,
        dx2
   };
   tolerance = 1.e-11;
   if(argc > 1)
   {
        tolerance = atof(argv[1]);
   }
//    start[0] = 1;
//    start[1] = 0;
//    start[2] = 1;
//    start[3] = 0;
   solutionUpToTime(start, end, functions, finish);
   printf("Error vector:(%lf,%lf)\n", fabs(end[2]), fabs(end[3]+1));
   for(int i = 0; i < 40; i++)
   {
        initialP1[i] = 0.1*i - 2;
        initialP2[i] = 0.1*i - 2;
   }
   double alpha = argc > 2 ? atof(argv[2]) : 1;
   for(int i = 0; i < 40; i++)
   { 
         fprintf(errors, "%lf,%lf\n", initialP1[i],error(initialP1[i], 0, functions, alpha));
   }
   return 0;
}
