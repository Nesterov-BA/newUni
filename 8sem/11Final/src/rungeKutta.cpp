#include "rungeKutta.hpp"
#include <cstdio>
#include <system_error>

double maxStep = 1;
double minStep = 0.0001;
double factor = 0.8;

double c[6] = {0.2, 0.3, 0.8, 0.888888, 1, 1};
double a[21] =
{
0.2,
0.075,    0.225,
0.977777, -3.733333,  3.555555,
2.952598, -11.595793, 9.822892, -0.290809,
2.846275, -10.757575, 8.906422, 0.278409, -0.273531,
0.091145, 0,          0.449236, 0.651041, -0.322376, 0.130952
};
double bErr[7] = {0.089913, 0, 0.453489, 0.614062, -0.271512, 0.089047, 0.025};


void dorPri5Varied(double* start, double* end, function* functions, double* step, double time)
{
    double err[4];
    double errNorm = 1;
    double stepFac;
    while(errNorm > tolerance)
    {
        dorPri5(start, end, err, functions, *step, time);
        for (int i = 0; i < 4; i++) {
            err[i] -= end[i];
        }
        errNorm = l1Norm(err, 4);
        stepFac = factor*pow((tolerance/errNorm),(double)1/(5+1));
        stepFac = minStep > stepFac ? minStep : stepFac;
        stepFac = maxStep < stepFac ? maxStep : stepFac;
        *step = *step*stepFac;
    }
    dorPri5(start, end, err, functions, *step, time);
}

void Runge_Kutta4StepVariedSimple(double* start, double* end, function* functions, double* step, double time)
{
    double err = 1;
    double* tempEnd = new double[4];
    double* errVector = new double[4];
    double stepFac;
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
}
//t \in (0, T)
// F(t) = f(T-t)
// G(t) = g(T-t)
// f'(t) = m(f(t), g(t))
// g'(t) = n(f(t), g(t))
// F'(t) = -f'(T-t) = -m(f(T-t), g(T-t)) = -m(F(t), G(t))
// G'(t) = -g'(T-t) = -n(f(T-t), g(T-t)) = -n(F(t), G(t))
void dorPri5(double* start, double* end, double* err, function* functions, double step, double time)
{
    double k1[4];
    double k2[4];
    double k3[4];
    double k4[4];
    double k5[4];
    double k6[4];
    double k7[4];
    double iniStart[4];
    double oldSt[4] = {start[0], start[1], start[2], start[3]};
    for(int i = 0; i < 4; i++)
    {
        k1[i] = functions[i](start[0], start[1], start[2], start[3], time);
    }

    for(int i = 0; i < 4; i++)
    {
        iniStart[i] = start[i] + a[0]*step*k1[i];
    }
    for(int i = 0; i < 4; i++)
    {
        k2[i] = functions[i](iniStart[0], iniStart[1], iniStart[2], iniStart[3], time + c[0]*step);
    }
    for(int i = 0; i < 4; i++)
    {
        iniStart[i] = start[i] + step*(a[1]*k1[i] + a[2]*k2[i]);
    }
    for(int i = 0; i < 4; i++)
    {
        k3[i] = functions[i](iniStart[0], iniStart[1], iniStart[2], iniStart[3], time + c[1]*step);
    }
    for(int i = 0; i < 4; i++)
    {
        iniStart[i] = start[i] + step*(a[3]*k1[i] + a[4]*k2[i] + a[5]*k3[i]);
    }
    for(int i = 0; i < 4; i++)
    {
        k4[i] = functions[i](iniStart[0], iniStart[1], iniStart[2], iniStart[3], time + c[2]*step);
    }

    for(int i = 0; i < 4; i++)
    {
        iniStart[i] = start[i] + step*(a[6]*k1[i] + a[7]*k2[i] + a[8]*k3[i] + a[9]*k4[i]);
    }
    for(int i = 0; i < 4; i++)
    {
        k5[i] = functions[i](iniStart[0], iniStart[1], iniStart[2], iniStart[3], time + c[3]*step);
    }
    for(int i = 0; i < 4; i++)
    {
        iniStart[i] = start[i] + step*(a[10]*k1[i] + a[11]*k2[i] + a[12]*k3[i] + a[13]*k4[i] + a[14]*k5[i]);
    }
    for(int i = 0; i < 4; i++)
    {
        k6[i] = functions[i](iniStart[0], iniStart[1], iniStart[2], iniStart[3], time + c[4]*step);
    }
    for(int i = 0; i < 4; i++)
    {
        iniStart[i] = start[i] + step*(a[15]*k1[i] + a[16]*k2[i] + a[17]*k3[i] + a[18]*k4[i] + a[19]*k5[i] + a[20]*k6[i]);
    }
    for(int i = 0; i < 4; i++)
    {
        k7[i] = functions[i](iniStart[0], iniStart[1], iniStart[2], iniStart[3], time + c[5]*step);
        end[i] = iniStart[i];
    }
    for(int i = 0; i < 4; i++)
    {
        err[i] = start[i] + step*(bErr[0]*k1[i] + bErr[1]*k2[i] + bErr[2]*k3[i] + bErr[3]*k4[i] + bErr[4]*k5[i] + bErr[5]*k6[i] + bErr[6]*k7[i]);
    }
    for (int i = 0; i < 4; i++)
    {
        start[i] = oldSt[i];
    }
}

void Runge_Kutta4ClassicSimple(double* start, double* end, function* functions, double step, double time)
{
    double* k1 = new double[4];
    double* k2 = new double[4];
    double* k3 = new double[4];
    double* k4 = new double[4];
    
    for(int i = 0; i < 4; i++)
    {
        k1[i] = functions[i](start[0], start[1], start[2], start[3], time);
    }
    for(int i = 0; i < 4; i++)
    {
        k2[i] = functions[i](start[0] + step*k1[0]/2, start[1] + step*k1[1]/2, start[2] + step*k1[2]/2, start[3] + step*k1[3]/2, time + step/2);
    }
    for(int i = 0; i < 4; i++)
    {
        k3[i] = functions[i](start[0] + step*k2[0]/2, start[1] + step*k2[1]/2, start[2] + step*k2[2]/2, start[3] + step*k2[3]/2, time + step/2);
    }
    for(int i = 0; i < 4; i++)
    {
        k4[i] = functions[i](start[0] + step*k3[0], start[1] + step*k3[1], start[2] + step*k3[2], start[3] + step*k3[3], time + step);
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
             Runge_Kutta4StepVariedSimple(start, tempEnd, functions, &step, time);
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
void solutionUpToTime(double* start, double* end, function* functions, double finish, string filename, double* integral)
 {
     ofstream file(filename);
     int numOfPoints = 0;
     double initialStart[4];
     double time = 0;
     double step = 0.1;
     double* tempEnd = new double[4];
     *integral = 0;
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
             Runge_Kutta4StepVariedSimple(start, tempEnd, functions, &step, time);
             time += step;
             *integral += (tempEnd[1]*tempEnd[1] + start[1]*start[1])*step/2;
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
