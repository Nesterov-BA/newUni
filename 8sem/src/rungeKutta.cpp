#include "rungeKutta.hpp"

int maxNumberOfCycles = 5;
bool limitedCycles = true;
// tolerance = 1.e-11;
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


void Runge_Kutta4StepVariedSimple(double* start, double* end, function* functions, double* step)
{
    double err = 1;
    double* tempEnd = new double[4];
    double* errVector = new double[4];
    double stepFac;
    while (err > tolerance)
    {
        Runge_Kutta4ClassicSimple(start, end, functions, *step);
        
        Runge_Kutta4ClassicSimple(start, tempEnd, functions, *step/2);
        Runge_Kutta4ClassicSimple(tempEnd, tempEnd, functions, *step/2);

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
        k2[i] = functions[i](start[0] + step*k1[0]/2, start[1] + step*k1[1]/2, start[2] + step*k1[2]/2, start[3] + step*k1[3]/2);
    }
    for(int i = 0; i < 4; i++)
    {
        k3[i] = functions[i](start[0] + step*k2[0]/2, start[1] + step*k2[1]/2, start[2] + step*k2[2]/2, start[3] + step*k2[3]/2);
    }
    for(int i = 0; i < 4; i++)
    {
        k4[i] = functions[i](start[0] + step*k3[0], start[1] + step*k3[1], start[2] + step*k3[2], start[3] + step*k3[3]);
    }
    for(int i = 0; i < 4; i++)
    {
        end[i] = start[i] + step/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
}

void solutionUpToTime(double* start, double* end, function* functions, double finish)
 {
     ofstream file("data.csv");
     ofstream fileErr("errorLog.txt");
     ofstream fileErrReg("errorRegular.txt");
     ofstream numberOfPoints("numberOfPoints.txt");
     int count = 0;
     int numOfPoints = 0;
     double time = 0;
     double step = 0.1;
     double* tempEnd = new double[4];
     double errorSum = 0;
     double globalError = 0;
     double globalErrorRegular = 0;
     file << "p1,p2,x1,x2,time" << endl;

     for (int i = 0; i < 4; i++) {
         file << start[i] << ",";
     }
     file << time << endl;
     while(time < finish)
     {
         if(time + step > finish)
         {
             step = (finish - time)/2;
             Runge_Kutta4ClassicSimple(start, tempEnd, functions, step);
             Runge_Kutta4ClassicSimple(tempEnd, tempEnd, functions, step);
             numOfPoints+=2;
             time += 2*step;
         }
         else
         {
             Runge_Kutta4StepVariedSimple(start, tempEnd, functions, &step);
             time += step;
             numOfPoints++;
         }

         //fileErr << globalError << " " << logNormCalc(xStart, yStart, tempX, tempY, step) << " " << time << endl;
         //fileErrReg << globalErrorRegular << " " << regNormCalc(xStart, yStart, tempX, tempY, step) << " " << time << endl;
         for (int i = 0; i < 4; i++)
         {
             file << tempEnd[i] << ",";
             start[i] = tempEnd[i];
         }
         file << time << endl;
     }
     for(int i = 0; i < 4; i++)
     {
         end[i] = tempEnd[i];
     }
     printf("Number of points: %d, average step: %e\n", numOfPoints, time/numOfPoints);
 }
