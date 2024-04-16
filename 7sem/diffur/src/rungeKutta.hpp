#include <bits/types/FILE.h>
#include <cmath>
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;
void clear();
void Runge_Kutta4Classic(double startX, double startY, double f(double, double), double g(double, double), double h, double* endX, double* endY);
void Runge_Kutta4StepVaried(double startX, double startY, double f(double, double), double g(double, double), double* newStep, double* endX, double* endY, double* errorSum);
void solutionUpToTime(double xStart, double yStart, double f(double, double), double g(double, double), double finish, double* cycleTimeLess, double* cycleTimeMore);
void findCycle(double xStart, double yStart, double f(double, double), double g(double, double), double* xEnd, double* yEnd);


int maxNumberOfCycles = 5;
bool limitedCycles = true;

double tolerance = 1.e-9;
double maxStep = 1;
double minStep = 0.1;
double factor = 0.9;
double* c = new double[3];
double* b = new double[3];
double* a = new double[6];
double* d = new double[4];



double a10 = 1./4.;
double a20 = -189./800.;
double a21 = 729./800.;
double a30 = 214./891.;
double a31 = 1./33.;
double a32 = 650./891.;
    

double c2 = 1./4.;
double c3 = 27./40.;
double c4 = 1.;
    
double d1 = 533./2106.;
double d2 = 0;
double d3 = 800./1053.;
double d4 = -1./78.;



void clear()
{
    ofstream file("data.txt");
    file<<"";
    file.close();
}

void Runge_Kutta4Classic(double startX, double startY, double f(double, double), double g(double, double), double h, double* endX, double* endY, double* error)
{
    // cout << "hello" << endl;
    double k1X, k1Y, k2X, k2Y, k3X, k3Y, k4X, k4Y;
    k1X = f(startX,                            startY);
    k1Y = g(startX,                            startY);
    k2X = f(startX + h * a10*k1X,              startY + h *  a10*k1Y);
    k2Y = g(startX + h * a10*k1X,              startY + h *  a10*k1Y);
    k3X = f(startX + h * (a20*k1X + a21*k2X),  startY + h * (a20*k1Y + a21*k2Y));
    k3Y = g(startX + h * (a20*k1X + a21*k2X),  startY + h * (a20*k1Y + a21*k2Y));
    double k1 = h*(a30*k1X + a31*k2X + a32*k3X);
    double k2 = h*(a30*k1Y + a31*k2Y + a32*k3Y);
    k4X = f(startX + k1,  startY +  k2);
    k4Y = g(startX + k1,  startY +  k2);

    double higherOrder1 = h*(d1*k1X + d2*k2X + d3*k3X + d4*k4X);
    double higherOrder2 = h*(d1*k1Y + d2*k2Y + d3*k3Y + d4*k4Y);
    // if(h < 100 && h > 0.01)
    //     {
    //         cout << "error in func = " << k1 - higherOrder1 << " + " << k2 - higherOrder2 << endl;
    //         cout << "startX, startY " << startX <<", " << startY << ", " << a21 << endl;
    //     }
    *error = fabs(k1 - higherOrder1) + fabs(k2 - higherOrder2);
    *endX = startX + k1;
    *endY = startY + k2;
}
   




void Runge_Kutta4StepVaried(double startX, double startY, double f(double, double), double g(double, double), double* newStep, double* endX, double* endY, double* errorSum)
{
    double err;
    double step = *newStep;
    Runge_Kutta4Classic(startX, startY, f, g, step, endX, endY, &err);
    //     cout << x2 << " " << y2 << " " << *endX << " " << *endY << endl;
    //  if(step < 100 && step > 0.01)
        //  cout << "err = " << err << " step = " << step << endl;
    *newStep = step*fmin(maxStep, fmax(minStep, factor*pow(tolerance/err, 0.25)));
    
    if(err > tolerance)
    {
        Runge_Kutta4StepVaried(startX, startY, f, g, newStep, endX, endY, errorSum);
    }
    else
    {
        *errorSum += err;
    }

}
//стартуем из (0,3), прогоняем до момента, пока не станет (0, -x_0)

void solutionUpToTime(double xStart, double yStart, double f(double, double), double g(double, double), double finish, double* cycleTimeLess, double* cycleTimeMore, double* xEnd, double* yEnd)
{
    ofstream file("data.txt");
    int count = 0;
    double time = 0;
    double step = 1;
    double tempX, tempY;
    double errorSum = 0;
    file << xStart << " " << yStart << " " << time << endl;
    while(time < finish)
    {
        Runge_Kutta4StepVaried(xStart, yStart, f, g, &step, &tempX, &tempY, &errorSum);
        time += step;
            if (yStart * tempY < 0 && ++count == 2) 
            {   
                *cycleTimeLess = time - step;
                *cycleTimeMore = time;

            }
        xStart = tempX;
        yStart = tempY;
        file << xStart << " " << yStart << " " << time << endl;
    }
    *xEnd = xStart;
    *yEnd = yStart;

}
void findCycle(double xStart, double yStart, double f(double, double), double g(double, double), double* xEnd, double* yEnd)
{
    xStart = 1;
    yStart = 0;
    double chordStart, chordEnd, chordTemp; 
    double cycleTimeLess, cycleTimeMore;
    double yEnd1, yEnd2;
    double xEnd1, xEnd2;
    double xStartTemp, yStartTemp;
    int count = 0;
    
    xStartTemp = xStart;
    yStartTemp = yStart;

    solutionUpToTime(xStart, yStart, f, g, 100, &cycleTimeLess, &cycleTimeMore, &xEnd1, &yEnd1);
    while(count < 10)
    {
        chordStart = cycleTimeLess;
        chordEnd = cycleTimeMore;
        while (chordEnd - chordStart > tolerance) 
        {
            solutionUpToTime(xStartTemp, yStartTemp, f, g, chordStart, &cycleTimeLess, &cycleTimeMore, &xEnd1, &yEnd1);
            solutionUpToTime(xStartTemp, yStartTemp, f, g, chordEnd, &cycleTimeLess, &cycleTimeMore, &xEnd2, &yEnd2);
            chordTemp = chordStart - (chordEnd-chordStart)*yEnd1/(yEnd2-yEnd1);
            chordStart = chordEnd;
            chordEnd = chordTemp;

            // x3 = x1 - (x2-x1)*y1[1]/(y2[1]-y1[1]);
        }
        cout << count << endl;
        cout << chordStart << " " << chordEnd << endl;
        if(fabs(xStartTemp - xEnd2) < tolerance*100)
        {
            *xEnd = xEnd2;
            *yEnd = yEnd2;
            cout << "Cycle found from " << xStartTemp << "," << yStartTemp << " to " << xEnd2 << "," << yEnd2 << endl;
            break;
        }
        xStartTemp = xEnd2;
        count++;
    }


}