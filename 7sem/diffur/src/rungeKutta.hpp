#include <bits/types/FILE.h>
#include <fstream>
#include <iostream>
#include "normCalculation.hpp"

using namespace std;
void clear();
void solutionUpToTime(double xStart, double yStart, double f(double, double), double g(double, double), double finish, double* cycleTimeLess, double* cycleTimeMore);
void findCycle(double xStart, double yStart, double f(double, double), double g(double, double), double* xEnd, double* yEnd);
void Runge_Kutta4ClassicSimple(double startX, double startY, double step, double f(double, double), double g(double, double), double* endX, double* endY);
void Runge_Kutta4StepVariedSimple(double startX, double startY, double f(double, double), double g(double, double), double* step, double* endX, double* endY, double* errorSum, double* globalError, double* globalErrorRegular);



int maxNumberOfCycles = 5;
bool limitedCycles = true;
double tolerance = 1.e-11;
double maxStep = 3;
double minStep = 0.00001;
double factor = 0.8;
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

void solutionUpToTime(double xStart, double yStart, double f(double, double), double g(double, double), double finish, double* cycleTimeLess, double* cycleTimeMore, double* xEnd, double* yEnd)
{
    ofstream file("data.txt");
    ofstream fileErr("errorLog.txt");
    ofstream fileErrReg("errorRegular.txt");
    ofstream numberOfPoints("numberOfPoints.txt");
    int count = 0;
    int numOfPoints = 0;
    double time = 0;
    double step = 0.1;
    double tempX, tempY;
    double errorSum = 0;
    double globalError = 0;
    double globalErrorRegular = 0;
    file << xStart << " " << yStart << " " << time << endl;
    while(time < finish)
    {
        if(time + 2*step > finish)
        {
            step = (finish - time)/2;
            Runge_Kutta4ClassicSimple(xStart, yStart, step, f, g, &tempX, &tempY);
            Runge_Kutta4ClassicSimple(tempX, tempY, step, f, g, &tempX, &tempY);
            numOfPoints+=2;

        }
        else
        {
            Runge_Kutta4StepVariedSimple(xStart, yStart, f, g, &step, &tempX, &tempY, &errorSum, &globalError, &globalErrorRegular);
            numOfPoints++;
        }

        time += 2*step;
        
        if (yStart * tempY < 0) 
        {
            printf("%.10f * %.10f, timeL = %f, timeR = %f, count = %d, errorSum = %e\n", yStart, tempY, time-2*step, time, count, errorSum);
        }
        if (yStart * tempY < 0 && ++count == 2) 
        {   
            *cycleTimeLess = time - 2*step;
            *cycleTimeMore = time;
        }
        fileErr << globalError << " " << logNormCalc(xStart, yStart, tempX, tempY, step) << " " << time << endl;
        fileErrReg << globalErrorRegular << " " << regNormCalc(xStart, yStart, tempX, tempY, step) << " " << time << endl;
        xStart = tempX;
        yStart = tempY;
        file << xStart << " " << yStart << " " << time << endl;
        
    }
    *xEnd = xStart;
    *yEnd = yStart;
    printf("Tolerance: %e Error: %e, numOfPoints: %d, globalError: %e, globalErrorReg: %e\n", tolerance, errorSum, numOfPoints, globalError, globalErrorRegular);
    numberOfPoints << numOfPoints;
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
    chordEnd = 10;
    xStartTemp = xStart;
    yStartTemp = yStart;

    while(count < 100)
    {
        solutionUpToTime(xStartTemp, yStartTemp, f, g, 50, &cycleTimeLess, &cycleTimeMore, &xEnd1, &yEnd1);
        chordStart = cycleTimeLess;
        chordEnd = cycleTimeMore;
        while (fabs(chordEnd - chordStart) > tolerance) 
        {
            solutionUpToTime(xStartTemp, yStartTemp, f, g, chordStart, &cycleTimeLess, &cycleTimeMore, &xEnd1, &yEnd1);
            solutionUpToTime(xStartTemp, yStartTemp, f, g, chordEnd, &cycleTimeLess, &cycleTimeMore, &xEnd2, &yEnd2);
            printf("(%.10f, %.10f), (%.10f, %.10f)\n", xEnd1, yEnd1, xEnd2, yEnd2);
            chordTemp = chordStart - (chordEnd-chordStart)*yEnd1/(yEnd2-yEnd1);
            chordStart = chordEnd;
            chordEnd = chordTemp;
            std::cout << "Cycle time: " << chordStart << " " << chordEnd << ", difference = " << -chordStart + chordEnd << endl;
        }
        printf("\nCycle end \n");

        if(fabs(xStartTemp - xEnd2) < tolerance*100)
        {
            *xEnd = xEnd2;
            *yEnd = yEnd2;
            cout << "Cycle found from " << xStartTemp << "," << yStartTemp << " to " << xEnd2 << "," << yEnd2 << endl;
            std::cout << "Cycle time: " << chordStart << " " << chordEnd << ", difference = " << -chordStart + chordEnd << endl;
            
            solutionUpToTime(xStartTemp, yStartTemp, f, g, chordStart, &cycleTimeLess, &cycleTimeMore, &xEnd1, &yEnd1);
            
            break;
        }
        xStartTemp = xEnd2;
        count++;
    }
    printf("\nCycle end, final count = %d\n", count);

}


void Runge_Kutta4ClassicSimple(double startX, double startY, double step, double f(double, double), double g(double, double), double* endX, double* endY)
{
    double k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y;
    
    k1x = f(startX, startY);
    k1y = g(startX, startY);

    k2x = f(startX + step/2*k1x, startY + step/2*k1y);
    k2y = g(startX + step/2*k1x, startY + step/2*k1y);
    
    k3x = f(startX + step/2*k2x, startY + step/2*k2y);
    k3y = g(startX + step/2*k2x, startY + step/2*k2y);

    k4x = f(startX + step*k3x, startY + step*k3y);
    k4y = g(startX + step*k3x, startY + step*k3y);

    *endX = startX + step/6*(k1x + 2*k2x + 2*k3x + k4x);
    *endY = startY + step/6*(k1y + 2*k2y + 2*k3y + k4y);

}

void Runge_Kutta4StepVariedSimple(double startX, double startY, double f(double, double), double g(double, double), double* step, double* endX, double* endY, double* errorSum, double* globalError, double* globalErrorRegular)
{
    double err = 1;
    double tempEndX, tempEndY;
    double stepFac;
    while (err > tolerance)
    {
        Runge_Kutta4ClassicSimple(startX, startY, 2*(*step), f, g, endX, endY);
        tempEndX = *endX;
        tempEndY = *endY;
        
        Runge_Kutta4ClassicSimple(startX, startY, *step, f, g, endX, endY);
        Runge_Kutta4ClassicSimple(*endX, *endY, *step, f, g, endX, endY);

        double dX = fabs(*endX+(*endX-tempEndX)/(16-1));
        double dY = fabs(*endY+(*endY-tempEndY)/(16-1));

        err = fabs((*endX-tempEndX)/dX) > fabs((*endY-tempEndY)/dY) ? fabs((*endX-tempEndX)/dX) : fabs((*endY-tempEndY)/dY);
        err /= (16-1);
        stepFac = factor*pow((tolerance/err),(double)1/(4+1));
        stepFac = minStep > stepFac ? minStep : stepFac;
        stepFac = maxStep < stepFac ? maxStep : stepFac;
        *step = *step*stepFac;
    }
    
    *globalError *= exp(logNormCalc(startX, startY, *endX, *endY, *step));
    *globalError += err;
    *globalErrorRegular *= exp(regNormCalc(startX, startY, *endX, *endY, *step));
    *globalErrorRegular += err;
    *errorSum += err;
  
}

void checkCycle(double f(double, double), double g(double, double))
{
    double xStart = 1;
    double yStart = 0;
    double xEnd, yEnd, cycleTimeLess, cycleTimeMore;
    double xEnd1, yEnd1, xEnd2, yEnd2;
    solutionUpToTime(xStart, yStart, f, g, 50, &cycleTimeLess, &cycleTimeMore, &xEnd, &yEnd);
    printf("\ncycle1: %.10f, %.10f\n", cycleTimeLess, cycleTimeMore);
    double tStart = cycleTimeLess;
    double tEnd = cycleTimeMore;
    solutionUpToTime(xStart, yStart, f, g, tStart, &cycleTimeLess, &cycleTimeMore, &xEnd1, &yEnd1);
    solutionUpToTime(xStart, yStart, f, g, tEnd, &cycleTimeLess, &cycleTimeMore, &xEnd2, &yEnd2);
    printf("ends: %.10f, %.10f\n", yEnd1, yEnd2);
}

