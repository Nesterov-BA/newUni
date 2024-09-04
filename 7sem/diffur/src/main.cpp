#include <iostream>
#include <math.h>
#include <fstream>
#include "rungeKutta.hpp"

using namespace std;



double dx(double x, double y);
double dy(double x, double y);

double dx(double x, double y)
{
    return y;
}
double dy(double x, double y)
{
    return /*alpha*/alpha*(1-x*x)*y-x;
}

int main (int argc, char** argv)
{
    if(argc != 3 || argc < 3)
    {
        cout << "Usage: " << argv[0] << " <xStart> <yStart>" << endl;
        return 1;
    }
    double xStart = atof(argv[1]);
    double yStart = atof(argv[2]);
    double xEnd, yEnd;
    printf("Set alpha:\n");
    std::cin >> alpha;
    int numberOfPoints = 0;
    // solutionUpToTime(xStart, yStart, dx, dy, 10, &cycleL, &cycleR, &xEnd, &yEnd);
    // checkCycle(dx, dy);
    // Runge_Kutta4StepVariedSimple(xStart, yStart, dx, dy, &step, &xEnd, &yEnd, &errSum);
    // printf("Step = %f\n", step);
    findCycle(xStart, yStart, dx, dy, &xEnd, &yEnd);
    cout << "Number of points: " << numberOfPoints << endl;
    

    return 0;
}
