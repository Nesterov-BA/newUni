#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


double leftBoundaryX = 0;
double rightBoundaryX = 1;
double leftBoundaryT = 0;
double rightBoundaryT = 2;


int numberOfXPoints = 50;
int numberOfTPoints = 50;

double** u = new double*[numberOfTPoints];

double** w = new double*[numberOfTPoints];

double wFunc(double x)
{
    return 1 - 16*x*x * (1 - x) * (1 - x);
}

void allocateMemory()
{
    for (int i =  0; i < numberOfTPoints; ++i) {
        u[i] = new double[numberOfXPoints];
    }


    for (int i =  0; i < numberOfTPoints; ++i) {
        w[i] = new double[numberOfXPoints];
    }
}
double* x = new double[numberOfXPoints];
double* t = new double[numberOfTPoints];

double dx = (rightBoundaryX - leftBoundaryX) / (numberOfXPoints - 1);
double dt = (rightBoundaryT - leftBoundaryT) / (numberOfTPoints - 1);

void createPoints()
{
    for (int i = 0; i < numberOfXPoints; i++) {
        x[i] = leftBoundaryX + i * dx;
    }
    for (int i = 0; i < numberOfTPoints; i++) {
        t[i] = leftBoundaryT + i * dt;
    }

    for(int i = 0; i < numberOfTPoints; i++)
    {
        for(int j = 0; j < numberOfXPoints; j++)
        {
            w[i][j] = wFunc(x[j]);
        }
    }
}

// d2u/t2 - d2u/dx2 + 2wu = 0
//  (u[i+1][j] - 2u[i][j] + u[i-1][j]) / dt*dt + (u[i][j+1] - 2u[i][j] + u[i][j-1]) / dx*dx  + 2w[i][j]*u[i][j] = 0  
// u[i+1][j] = -dt*dt *((u[i][j+1] - 2u[i][j] + u[i][j-1]) / dx*dx  + 2w[i][j]*u[i][j]) + 2u[i][j] - u[i-1][j]]

// (u[i+1][0] - u[i][0] )/dt= (u[i][1] - u[i][0])/dx + sin(t[i])
// (u[i+1][numberOfXPoints-1] - u[i][numberOfXPoints-1] )/dt= (u[i][numberOfXPoints-1] - u[i][numberOfXPoints-2])/dx
void networkMethod() 
{
    allocateMemory();
    createPoints();
    for (int j =  0; j < numberOfTPoints; ++j) {
        u[0][j] =  0;
        u[1][j] = 0; // u is zero at t=0
    }
    for(int i = 2; i < numberOfTPoints - 1; ++i)
    {
        u[i+1][0] = dt * ((u[i][1] - u[i][0])/dx + sin(t[i])) + u[i][0];
    
        for(int j = 1; j < numberOfXPoints - 1; ++j)
            u[i+1][j] = -dt*dt *((u[i][j+1] - 2*u[i][j] + u[i][j-1]) / (dx*dx) + 2*w[i][j]*u[i][j]) + 2*u[i][j] - u[i-1][j];
    
        u[i+1][numberOfXPoints-1] = dt * ((u[i][numberOfXPoints-1] - u[i][numberOfXPoints-2])/dx) + u[i][numberOfXPoints-1];
    }
}


void writeToDataFile()
{
    ofstream dataFile;
    dataFile.open("data.txt");
    networkMethod();
    for(int i = 0; i < numberOfTPoints; ++i)
    {
        for(int j = 0; j < numberOfXPoints; ++j)
        {
            dataFile << x[j] << " " << t[i] << " " << u[i][j] << endl;
        }
    }
    dataFile.close();
}



