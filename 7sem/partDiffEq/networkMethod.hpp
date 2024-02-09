#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


double leftBoundaryX = 0;
double rightBoundaryX = 1;
double leftBoundaryT = 0;
double rightBoundaryT = 1;




void pointsOfInterpolation(int n, int m, double* xPoints, double* tPoints) 
{
    double dx = (rightBoundaryX - leftBoundaryX) / (n - 1);
    double dt = (rightBoundaryT - leftBoundaryT) / (m - 1);
    for (int i = 0; i < n; i++) {
        x[i] = leftBoundaryX + i * dx;
    }
    for (int i = 0; i < m; i++) {
        t[i] = leftBoundaryT + i * dt;
    }
}

void networkMethod(int n, int m, double* xPoints, double* tPoints) 
{
    double dx = (rightBoundaryX - leftBoundaryX) / (n - 1);
    double dt = (rightBoundaryT - leftBoundaryT) / (m - 1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            u[i][j] = 0;
        }
    }
    for (int i = 1; i < n - 1; i++) {
        for (int j = 1; j < m - 1; j++) {
            u[i][j] = (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1]) / 4;
        }
    }
}





