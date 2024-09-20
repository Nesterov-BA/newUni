#include "rungeKutta.hpp"

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


void Runge_Kutta4StepVariedSimple(double* start, double* end, function* functions, double finish, double* step)
{
    double err = 1;
    double* tempEnd = new double[4];
    double* errVector = new double[4];
    double stepFac;
    while (err > tolerance)
    {
        Runge_Kutta4ClassicSimple(start, end, functions, 2*(*step));
        
        Runge_Kutta4ClassicSimple(start, tempEnd, functions, *step);
        Runge_Kutta4ClassicSimple(tempEnd, tempEnd, functions, *step);

        for(int i = 0; i < 4; i++)
        {
            errVector[i] = fabs(tempEnd[i]+(tempEnd[i]-end[i])/(16-1));
            errVector[i] = fabs(end[i] - tempEnd[i]) / errVector[i];
        }

        err = lInfNorm(errVector, 4);
        err /= (16-1);
        stepFac = factor*pow((tolerance/err),(double)1/(4+1));
        stepFac = minStep > stepFac ? minStep : stepFac;
        stepFac = maxStep < stepFac ? maxStep : stepFac;
        *step = *step*stepFac;
    }
}

void Runge_Kutta4ClassicSimple(double* start, double* end, function* functions, double step)
{
    double* k1 = new double[4];
    double* k2 = new double[4];
    double* k3 = new double[4];
    double* k4 = new double[4];
    
    for(int i = 0; i < 4; i++)
    {
        k1[i] = functions[i](start[0], start[1], start[2], start[3]);
    }
    for(int i = 0; i < 4; i++)
    {
        k2[i] = functions[i](start[0] + step*k1[i]/2, start[1] + step*k1[i]/2, start[2] + step*k1[i]/2, start[3] + step*k1[i]/2);
    }
    for(int i = 0; i < 4; i++)
    {
        k3[i] = functions[i](start[0] + step*k2[i]/2, start[1] + step*k2[i]/2, start[2] + step*k2[i]/2, start[3] + step*k2[i]/2);
    }
    for(int i = 0; i < 4; i++)
    {
        k4[i] = functions[i](start[0] + step*k3[i], start[1] + step*k3[i], start[2] + step*k3[i], start[3] + step*k3[i]);
    }
    for(int i = 0; i < 4; i++)
    {
        end[i] = start[i] + step/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
}

// void solutionUpToTime(double* start, double* end, double p1(double), double finish)
// {
//     ofstream file("data.txt");
//     ofstream fileErr("errorLog.txt");
//     ofstream fileErrReg("errorRegular.txt");
//     ofstream numberOfPoints("numberOfPoints.txt");
//     int count = 0;
//     int numOfPoints = 0;
//     double time = 0;
//     double step = 0.1;
//     double tempX, tempY;
//     double errorSum = 0;
//     double globalError = 0;
//     double globalErrorRegular = 0;
//     file << xStart << " " << yStart << " " << time << endl;
//     while(time < finish)
//     {
//         if(time + 2*step > finish)
//         {
//             step = (finish - time)/2;
//             Runge_Kutta4ClassicSimple(xStart, yStart, step, f, g, &tempX, &tempY);
//             Runge_Kutta4ClassicSimple(tempX, tempY, step, f, g, &tempX, &tempY);
//             numOfPoints+=2;

//         }
//         else
//         {
//             Runge_Kutta4StepVariedSimple(xStart, yStart, f, g, &step, &tempX, &tempY, &errorSum, &globalError, &globalErrorRegular);
//             numOfPoints++;
//         }

//         time += 2*step;
//         fileErr << globalError << " " << logNormCalc(xStart, yStart, tempX, tempY, step) << " " << time << endl;
//         fileErrReg << globalErrorRegular << " " << regNormCalc(xStart, yStart, tempX, tempY, step) << " " << time << endl;
//         xStart = tempX;
//         yStart = tempY;
//         file << xStart << " " << yStart << " " << time << endl;
//     }
//     *xEnd = xStart;
//     *yEnd = yStart;
//     printf("Tolerance: %e Error: %e, numOfPoints: %d, globalError: %e, globalErrorReg: %e\n", tolerance, errorSum, numOfPoints, globalError, globalErrorRegular);
//     numberOfPoints << numOfPoints;


// }
