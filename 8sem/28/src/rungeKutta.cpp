#include "rungeKutta.hpp"
#include <cstdio>
#include <system_error>

double maxStep = 3;
double minStep = 0.0001;
double factor = 0.8;





void Runge_Kutta4StepVariedSimple(double* start, double* end, function* functions, double* step, double time, double* globalErr)
{
    double err = 1;
    double* tempEnd = new double[4];
    double* errVector = new double[4];
    double stepFac;
    double li, lii;
    while (err > tolerance)
    {
        Runge_Kutta4ClassicSimple(start, end, functions, *step, time);
        
        Runge_Kutta4ClassicSimple(start, tempEnd, functions, *step/2, time);
        Runge_Kutta4ClassicSimple(tempEnd, tempEnd, functions, *step/2, time + *step/2);

        for(int i = 0; i < 4; i++)
        {
            errVector[i] = fabs(tempEnd[i]+(tempEnd[i]-end[i])/(16-1));
            errVector[i] = fabs(end[i] - tempEnd[i]) / errVector[i];
        }

        err = lInfNorm(errVector, 2);
        err /= (16-1);
        if(fabs(err) < tolerance)
            break;
        stepFac = factor*pow((tolerance/err),(double)1/(4+1));
        stepFac = minStep > stepFac ? minStep : stepFac;
        stepFac = maxStep < stepFac ? maxStep : stepFac;
        *step = *step*stepFac;
    }
    li = abs(alpha*start[0]*sin(alpha*start[0]) + (alpha*start[0]*start[0]/4 - 0.5)*cos(alpha*start[0] + 0.5));
    lii = abs(alpha*end[0]*sin(alpha*end[0]) + (alpha*end[0]*end[0]/4 - 0.5)*cos(alpha*end[0] + 0.5));
    li = max(li, lii);
    li *= *step;
    *globalErr *= exp(li);
    *globalErr += err;
}
//t \in (0, T)
// F(t) = f(T-t)
// G(t) = g(T-t)
// f'(t) = m(f(t), g(t))
// g'(t) = n(f(t), g(t))
// F'(t) = -f'(T-t) = -m(f(T-t), g(T-t)) = -m(F(t), G(t))
// G'(t) = -g'(T-t) = -n(f(T-t), g(T-t)) = -n(F(t), G(t))


void Runge_Kutta4ClassicSimple(double* start, double* end, function* functions, double step, double time)
{
    double* k1 = new double[4];
    double* k2 = new double[4];
    double* k3 = new double[4];
    double* k4 = new double[4];
    time = time;
    
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
     int numOfPoints = 0;
     double initialStart[4];
     double time = 0;
     double step = 0.1;
     double* tempEnd = new double[4];
     double globalError = 0;
     file << "p1,p2,x1,x2,time" << endl;
     for (int i =0; i < 4; i++)
     {
         initialStart[i] = start[i];
     }
     for (int i = 0; i < 4; i++) {
         file << start[i] << ",";
     }
     file << time << endl;
     while(time < finish)
     {
         if(time + step > finish)
         {
             step = (finish - time)/2;
//             dorPri5(start, tempEnd, err, functions, step, time);
 //            dorPri5(tempEnd, end, err, functions, step, time);
             Runge_Kutta4ClassicSimple(start, tempEnd, functions, step, time);
             Runge_Kutta4ClassicSimple(tempEnd, tempEnd, functions, step, time + step);
             numOfPoints+=2;
             time += 2*step;
         }
         else
         {
//             dorPri5Varied(start, tempEnd, functions, &step, time);
             Runge_Kutta4StepVariedSimple(start, tempEnd, functions, &step, time, &globalError);
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
         start[i] = initialStart[i];
         end[i] = tempEnd[i];
     }
     //printf("Number of points: %d, average step: %e\n", numOfPoints, time/numOfPoints);
 }
void solutionUpToTime(double* start, double* end, function* functions, double finish, double* globalError)
 {
     ofstream file("data.csv");
     int numOfPoints = 0;
     double initialStart[4];
     double time = 0;
     double step = 0.1;
     double* tempEnd = new double[4];
     file << "x1,x2,time" << endl;
     for (int i =0; i < 4; i++)
     {
         initialStart[i] = start[i];
     }
     for (int i = 0; i < 4; i++) {
         file << start[i] << ",";
     }
     file << time << endl;
     while(time < finish)
     {
         if(time + step > finish)
         {
             step = (finish - time)/2;
//             dorPri5(start, tempEnd, err, functions, step, time);
 //            dorPri5(tempEnd, end, err, functions, step, time);
             Runge_Kutta4ClassicSimple(start, tempEnd, functions, step, time);
             Runge_Kutta4ClassicSimple(tempEnd, tempEnd, functions, step, time + step);
             numOfPoints+=2;
             time += 2*step;
         }
         else
         {
//             dorPri5Varied(start, tempEnd, functions, &step, time);
             Runge_Kutta4StepVariedSimple(start, tempEnd, functions, &step, time, globalError);
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
         start[i] = initialStart[i];
         end[i] = tempEnd[i];
     }
     //printf("Number of points: %d, average step: %e\n", numOfPoints, time/numOfPoints);
 }
