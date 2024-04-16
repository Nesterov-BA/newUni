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
    return /*alpha*/10*(1-x*x)*y-x;
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
    double cycleL, cycleR;
    int numberOfPoints = 0;
    // solutionUpToTime(xStart, yStart, dx, dy, 50, &cycleL, &cycleR);
    findCycle(xStart, yStart, dx, dy, &xEnd, &yEnd);
    cout << "Number of points: " << numberOfPoints << endl;
    cout << "Cycle time (less): " << cycleL << endl;
    cout << "Cycle time (more): " << cycleR << endl;

    return 0;
}