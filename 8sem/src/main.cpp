#include "rungeKutta.hpp"
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

int main(int argc, char** argv) 
{
   function functions[] = 
   {
        dp1,
        dp2,
        dx1,
        dx2
   };
   double p1 = 1;
   double p2 = 0;
   double x1 = 1;
   double x2 = 0;
   for(int i = 0; i < 4; i++)
        printf("%lf\n", functions[i](p1, p2, x1, x2));
    return 0;
}
    