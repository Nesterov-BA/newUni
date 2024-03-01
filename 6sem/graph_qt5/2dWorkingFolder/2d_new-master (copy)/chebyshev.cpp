#include "chebyshev.hpp"

#include "help.hpp"
#include "math.h"
#include <cstdio>

#define PI 3.1415926535897932384626433832795

void points(double *cx, double *cy, int nx, int ny, double a, double b, double c, double d)
{
    double stepX = (b-a)/(nx-1);
    double stepY = (d-c)/(ny-1);
    for(int i = 0; i<nx; i++)
    {
        cx[i] = a+i*stepX;
    }

    for(int i = 0; i<ny; i++)
    {
        cy[i] = c+i*stepY;
    }
}


void Fill_F(double *F, int nx, int ny, double *cx, double *cy, double (*f)(double, double))
{
    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            F[j*nx + i] = f(cx[i], cy[j]);
            //  printf("F(%f, %f) = %f\n", cx[i], cy[j], F[j*nx+i]);
        }
    }
}

void Fill_d2F(double *d2FX, double *d2FY, int nx, int ny,  double *cx, double *cy, double (*d2f)(double, double, int))
{
    for(int i = 0; i < nx; i++)
    {
        d2FY[i] = d2f(cx[i], cy[0], 1);
        d2FY[i+nx] = d2f(cx[i], cy[ny-1], 1);
    }
    for(int i = 0; i < ny; i++)
    {
        d2FX[i] = d2f(cx[0], cy[i], 0);
        d2FX[i+ny] = d2f(cx[nx-1], cy[i], 0);
    } // first amd last rows of interpolation points left to right or bottom to top respectively
}

void evalGamma(double* gamma, double stepX, double stepY, double* F, int nX, int nY)
{
    double* A = new double[16];
    double* B = new double[16];
    double* tempF = new double[16];

    double cX = 1/(stepX*stepX);
    double cY = 1/(stepY*stepY);

    for(int i = 0; i < 8; i++)
    {
        A[i] = 0;
        B[i] = 0;
    }
    A[0] = 1;
    B[0] = 1;
    A[5] = 1;
    B[5] = 1;

    A[8] = -3*cX/stepX;
    B[8] = -3*cY/stepY;
    A[9] = -2*cX;
    B[9] = -2*cY;
    A[10] = 3*cX/stepX;
    B[10] = 3*cY/stepY;
    A[11] = -cX;
    B[11] = -cY;
    A[12] = 2*cX/stepX;
    B[12] = 2*cY/stepY;
    A[13] = cX;
    B[13] = cY;
    A[14] = -2*cX/stepX;
    B[14] = -2*cY/stepY;
    A[15] = cX;
    B[15] = cY;
    transponse(B, 4);
    for(int i1 = 0; i1 < nY; i1++)
    {
        for(int j1 = 0; j1 < nX; j1++)
        {
            for(int i = 0; i<16; i++)
            {
                tempF[i] = F[i1*nX*16 + j1*16 +i];
            }
            // printf("A = \n");
            // print_matrix(A, 4, 4);
            // printf("B = \n");
            // print_matrix(B, 4, 4);
            // printf("tempF = \n");
            // print_matrix(tempF, 4, 4);

            sqMatrMult(tempF, B, 4, 0);

            sqMatrMult(A, tempF, 4, 1);
            // printf("tempFAfter = \n");
            // print_matrix(tempF, 4, 4);

            for(int i = 0; i<16; i++)
            {
                gamma[i1*nX*16 + j1*16 +i] = tempF[i];
            }
        }
    }
}

void dCoeff(double* d, double* values, int size, double step, double d2ValBeg, double d2ValEnd)
{
    for (int i = 1; i < size-1; i++)
    {
        d[i] = (2*values[i] - values[i-1] - values[i+1])/(2*step);
    }
    d[0] = ((3*(values[1] - values[0])/step) - d2ValBeg*step/2 - d[1])/2;
    d[size-1] = (3*(values[size-1] - values[size-2])*step + d2ValEnd*(step)/2 - d[size-2])/2;

}

void evalFx(double* F, double* Fx, int nX, int nY , double step, double* d2Fy)
{
    double* fCol = new double[nY];
    double* tempD = new double[nY];
    for(int colCount = 0; colCount < nX; colCount++)
    {
        for(int i = 0; i < nY; i++)
            fCol[i] = F[i*nX + colCount];
        dCoeff(tempD, fCol, nY, step, d2Fy[colCount], d2Fy[colCount + nX]);
        for(int i = 0; i < nY; i++)
            Fx[i*nX+colCount] = tempD[i];
    }
    delete [] fCol;
    delete [] tempD;
}

void evalFy(double* F, double* Fy, int nX, int nY , double step, double* d2Fx)
{
    double* fRow = new double[nX];
    double* tempD = new double[nX];
    for(int rowCount = 0; rowCount < nY; rowCount++)
    {
        for(int i = 0; i < nX; i++)
            fRow[i] = F[rowCount*nX + i];
        dCoeff(tempD, fRow, nX, step, d2Fx[rowCount], d2Fx[rowCount + nY]);
        for(int i = 0; i < nX; i++)
            Fy[rowCount*nX + i] = tempD[i];
    }
    delete [] fRow;
    delete [] tempD;
}

void evalFxy(double* Fy, double* Fxy, int nX, int nY , double step, double* d2Fy)
{
    evalFx(Fy, Fxy, nX, nY, step, d2Fy);
}

void evalFij(double* F, double* Fx, double* Fy, double* Fxy, double* tensorF, int nX, int nY)
{
    for(int i = 0; i < nY-1; i++)
    {
        for(int j = 0; j < nX-1; j++)
        {
            tensorF[i*nX*16 + j*16] = F[i*nX + j];
            tensorF[i*nX*16 + j*16 + 2] = F[i*nX + j+1];
            tensorF[i*nX*16 + j*16 + 8] = F[(i+1)*nX +j];
            tensorF[i*nX*16 + j*16 + 10] = F[(i+1)*nX +j+1];

            tensorF[i*nX*16 + j*16 + 1] = Fy[i*nX +j];
            tensorF[i*nX*16 + j*16 + 3] = Fy[i*nX +j +1];
            tensorF[i*nX*16 + j*16 + 9] = Fy[(i+1)*nX +j];
            tensorF[i*nX*16 + j*16 + 11] = Fy[(i+1)*nX +j+1];

            tensorF[i*nX*16 + j*16 + 4] = Fx[i*nX +j];
            tensorF[i*nX*16 + j*16 + 6] = Fx[i*nX +j +1];
            tensorF[i*nX*16 + j*16 + 12] = Fx[(i+1)*nX +j];
            tensorF[i*nX*16 + j*16 + 14] = Fx[(i+1)*nX +j+1];

            tensorF[i*nX*16 + j*16 + 5] = Fxy[i*nX +j];
            tensorF[i*nX*16 + j*16 + 7] = Fxy[i*nX +j +1];
            tensorF[i*nX*16 + j*16 + 13] = Fxy[(i+1)*nX +j];
            tensorF[i*nX*16 + j*16 + 15] = Fxy[(i+1)*nX +j+1];
        }
    }
}

void basis(double x, double y, double xi, double yj, double* arr)
{
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            arr[i*4 + j] = pow(x - xi,j)*pow(y - yj,i);
        }
    }
}

void interpolatetedVal(double* gamma, int xInt, int yInt, double* intVal, double xStart, double xEnd, double yStart, double yEnd, int intPointsNum, int nX, int nY)
{
    double* basFuncs = new double[16];
    basis(xStart, yStart, xStart, yStart, basFuncs);
            for(int k = 0; k < 16; k++)
                intVal[yInt*nX + xInt] += gamma[xInt*nY*16 + yInt*16 + k]*basFuncs[k];

    delete [] basFuncs;
}

