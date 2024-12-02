#include <iostream>
#include <math.h>
#include <fstream>
#include "rungeKutta.hpp"



double alpha;
double dx(double x, double y);
double dy(double x, double y);

double dx(double x, double y)
{
    x = x;
    return y;
}
double dy(double x, double y)
{
    return /*alpha*/alpha*(1-x*x)*y-x;
}

int main (int argc, char** argv)
{
    double xStart = 2;
    double yStart = 0;
    double xEnd, yEnd;
    double cycleL, cycleR;
    alpha = 0.1;
    if(argc > 1)
        alpha = atof(argv[1]);
    solutionUpToTime(xStart, yStart, dx, dy, 10, &cycleL, &cycleR, &xEnd, &yEnd);
    // checkCycle(dx, dy);
    // Runge_Kutta4StepVariedSimple(xStart, yStart, dx, dy, &step, &xEnd, &yEnd, &errSum);
    // printf("Step = %f\n", step);
    //findCycle(xStart, yStart, dx, dy, &xEnd, &yEnd);
    //fasterFindCycle(xStart, yStart, dx, dy, &xEnd, &yEnd);

    return 0;
}
