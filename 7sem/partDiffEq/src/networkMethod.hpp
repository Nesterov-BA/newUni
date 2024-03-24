#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;


double leftBoundaryX = 0;
double rightBoundaryX = 1;
double leftBoundaryT = 0;
double rightBoundaryT = 5;


int numberOfXPoints = 50;
int numberOfTPoints = rightBoundaryT*3*numberOfXPoints;

double** u = new double*[numberOfTPoints];

double** w = new double*[numberOfTPoints];

double* f = new double[numberOfTPoints];

double** dudt = new double*[numberOfTPoints];

double** dudx = new double*[numberOfTPoints];

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

    for (int i =  0; i < numberOfTPoints; ++i) {
        dudt[i] = new double[numberOfXPoints];
    }

    for (int i =  0; i < numberOfTPoints; ++i) {
        dudx[i] = new double[numberOfXPoints];
    }
}
double* x = new double[numberOfXPoints];
double* t = new double[numberOfTPoints];

double dx = (rightBoundaryX - leftBoundaryX) / (numberOfXPoints - 1);
double dt = (rightBoundaryT - leftBoundaryT) / (numberOfTPoints - 1);

double k = dt/dx;

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

void networkMethod() 
{
    allocateMemory();
    createPoints();
    for (int j =  0; j < numberOfXPoints; ++j) {
        u[0][j] =  0;
        u[1][j] = 0;
    }
    for(int i = 2; i < numberOfTPoints - 1; ++i)
    {
        u[i+1][0] = (2 * dt * dt / (dt + dx)) * ((u[i][1] - u[i][0])/dx + sin(t[i]) + u[i-1][0] / (2*dt) - dx * (2*w[i][0]*u[i][0] + (u[i-1][0] - 2 * u[i][0]) / (dt*dt))/2) ;
    

        for(int j = 1; j < numberOfXPoints - 1; ++j)
            u[i+1][j] = dt*dt *((u[i][j+1] - 2*u[i][j] + u[i][j-1]) / (dx*dx) - 2*w[i][j]*u[i][j]) + 2*u[i][j] - u[i-1][j];
    
        u[i+1][numberOfXPoints-1] = (2 * dt * dt / (dt + dx)) * ((u[i][numberOfXPoints-1] - u[i][numberOfXPoints - 2])/dx + u[i-1][numberOfXPoints-1] / (2*dt) - dx * (2*w[i][numberOfXPoints-1]*u[i][numberOfXPoints-1] + (u[i-1][numberOfXPoints-1] - 2 * u[i][numberOfXPoints-1]) / (dt*dt))/2);
    }
}




void filldudt()
{
    cout << "dudt" << endl;
    cout << 1/dt << endl;
    for(int i = 0; i < numberOfXPoints; ++i)
    {
        dudt[0][i] = 0;
    }
    
    for(int i = 1; i < numberOfTPoints; ++i)
    {
        for(int j = 0; j < numberOfXPoints; ++j)
        {
            dudt[i][j] = (u[i][j] - u[i-1][j])/dt;
        }
    }

}

void filldudx()
{
    cout << "dudx" << endl;
    for(int i = 0; i < numberOfTPoints; ++i)
    {
        dudx[i][0] = dudt[i][0] - sin(t[i]);
        dudx[i][numberOfXPoints-1] = -dudt[i][numberOfXPoints-1];

        for(int j = 1; j < numberOfXPoints - 1; ++j)
        {
            dudx[i][j] = (u[i][j+1] - u[i][j-1])/(2*dx);
        }
    }
}

void integrate()
{
    double integralValue = 0;
    for(int i = 0; i < numberOfTPoints; ++i)
    {
        integralValue = 0;
        for (int j = 0; j < numberOfXPoints-1; ++j)
        {
            integralValue += (dudt[i][j]*dudt[i][j] + dudt[i][j+1]*dudt[i][j+1])*dx/2;
            integralValue += (dudx[i][j]*dudx[i][j] + dudx[i][j+1]*dudx[i][j+1])*dx/2;
            integralValue += (w[i][j]*u[i][j]*u[i][j] + w[i][j+1]*u[i][j+1]*u[i][j+1])*dx;
        }
        f[i] = integralValue;
    }
}


void writeToDataFile()
{
    ofstream dataFile;
    ofstream dudtFile;
    ofstream dudxFile;
    ofstream fFile;
    dataFile.open("data.txt");
    dudtFile.open("dudt.txt");
    dudxFile.open("dudx.txt");
    fFile.open("f.txt");
    networkMethod();
    filldudt();
    filldudx();
    integrate();
    cout << u[numberOfTPoints/2][numberOfXPoints/2] << endl;
    for(int i = 0; i < numberOfTPoints; ++i)
    {
        for(int j = 0; j < numberOfXPoints; ++j)
        {
            dataFile << x[j] << " " << t[i] << " " << u[i][j] << endl;
            dudtFile << x[j] << " " << t[i] << " " << dudt[i][j] << endl;
            dudxFile << x[j] << " " << t[i] << " " << dudx[i][j] << endl;
        }
        fFile << t[i] << " " << f[i] << endl;
    }
    dataFile.close();
    dudtFile.close();
    dudxFile.close();
}



