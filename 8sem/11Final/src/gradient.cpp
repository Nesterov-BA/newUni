#include "metricFuncs.hpp"
#include "gradient.hpp"
#include <cstdio>
    double eps = 1.e-6;
void gradient(double (*func)(double, double), double *start, double *gradient)
{
    double value = func(start[0], start[1]);
    gradient[0] = -value + func(start[0] + eps, start[1]);
    gradient[1] = -value + func(start[0],  start[1] + eps);
    gradient[0] /= eps;
    gradient[1] /= eps;
    printf("Gradient in func: %lf, %lf\n", gradient[0], gradient[1]);
}

double derivative(double func(double, double), double* start, double* direction)
{
    double eps = 1.e-5;
    return (func(start[0] + eps*direction[0], start[1] + eps*direction[1]) - func(start[0], start[1]))/(eps*l1Norm(direction,2));
}

void newton2d(double func(double, double), double* start)
{
    double grad[2];
    int count = 0;
    double value;
    double deriv;
    FILE* newtonFile = fopen("newton.csv", "w");
    fprintf(newtonFile, "p1, p2, val\n");
    while(count < 20)
    {
        value = func(start[0], start[1]);
        gradient(func, start, grad);
        deriv = derivative(func, start, grad);
        printf("\nNow at: %lf, %lf\n", start[0], start[1]);
        printf("Derivative: %lf, Value: %lf\n", deriv, value);
        printf("Gradient:%lf, %lf\n", grad[0], grad[1]);
        printf("Difference: %lf - %lf =%.10f \n", func(start[0] + eps*grad[0], start[1] + eps*grad[1]),func(start[0] , start[1] ), func(start[0] + eps*grad[0], start[1] + eps*grad[1]) - func(start[0], start[1]));
        fprintf(newtonFile, "%lf, %lf, %lf\n", start[0], start[1], value);
        start[0] -= 0.5*grad[0]*value/deriv;
        start[1] -= 0.5*grad[1]*value/deriv;
        count++;
        if (value < 1.e-5)
            break;
    }
}

void gradient_descent(double (*func)(double, double), double *start)
{
    double gradVector[2];
    FILE* gradientFile = fopen("gradient.csv", "w");
    fprintf(gradientFile, "p1, p2, val\n");
    int count = 0;
    while (count < 2000)
    {
        gradient(func, start, gradVector);
        fprintf(gradientFile, "%lf, %lf, %lf\n", start[0], start[1], func(start[0], start[1]));
        start[0] -= gradVector[0]/5;
        start[1] -= gradVector[1]/5;
        count++;
        if(l1Norm(gradVector, 2) < 1.e-3)
            break;
    }
}
