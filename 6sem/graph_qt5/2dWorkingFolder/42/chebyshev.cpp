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
    double* A = new double[9];
    double* B = new double[9];
    double* tempF = new double[9];

    double cX = 1/(stepX*stepX);
    double cY = 1/(stepY*stepY);


    A[0] = 1;
    B[0] = 1;
    A[1] = 0;
    B[1] = 0;
    A[2] = 0;
    B[2] = 0;
    A[3] = -3/stepX;
    B[3] = -3/stepY;
    A[4] = 4/stepX;
    B[4] = 4/stepY;
    A[5] = -1/stepX;
    B[5] = -1/stepY;
    A[6] = 2*cX;
    B[6] = 2*cY;
    A[7] = -4*cX;
    B[7] = -4*cY;
    A[8] = 2*cX;
    B[8] = 2*cY;
    // printf("A = \n");
    // print_matr(A, 3, 3);
    transponse(B, 3);
    for(int i1 = 0; i1 < nY; i1++)
    {
        for(int j1 = 0; j1 < nX; j1++)
        {
            for(int i = 0; i<9; i++)
            {
                tempF[i] = F[i1*nX*9 + j1*9 +i];
            }
            // printf("A = \n");
            // print_matrix(A, 4, 4);
            // printf("B = \n");
            // print_matrix(B, 4, 4);
            // printf("(%d, %d)tempF = \n", j1, i1);
            // print_matr(tempF, 3, 3);

            sqMatrMult(tempF, B, 3, 0);

            sqMatrMult(A, tempF, 3, 1);
            // printf("tempFAfter = \n");
            // print_matr(tempF, 3, 3);

            for(int i = 0; i<9; i++)
            {
                gamma[i1*nX*9 + j1*9 +i] = tempF[i];
            }
        }
    }
    // print_matr(gamma, nX, nY*9);
}

void thomas_algorithm(double* a,
                      double* b,
                      double* c,
                      double* d,
                      double* f, int N) {

    double* c_star = new double[N-1];
    double* d_star = new double[N];
    double m;

    c_star[0] = c[0]/b[0];
    d_star[0] = d[0]/b[0];

    for(int i = 1; i < N; i++)
    {
        m = b[i] - a[i-1]*c_star[i-1];
        if(i < N-1)
            c_star[i] = c[i]/m;
        d_star[i] = (d[i] - a[i-1]*d_star[i-1])/m;
    }

    f[N-1] = d_star[N-1];
    for(int i = N-2; i > -1; i--)
        f[i] = d_star[i] - c_star[i]*f[i+1];
    delete[] c_star;
    delete[] d_star;
}

void vCoeff(double* v, double* values, int size, double step, double d2ValBeg, double d2ValEnd)
{
    double* a = new double[size];
    double* b = new double[size+1];
    double* c = new double[size];
    double* d = new double[size+1];

    b[0] = 4;
    c[0] = b[0];
    a[size-1] =b[0];
    b[size] = b[0];
    d[0] = values[0]*8/* + step*step*d2ValBeg*/;
    d[size] = values[size-1]*8/* + step*step*d2ValEnd*/;
    for(int row = 1; row < size; row++)
    {
        a[row-1] =  1;
        b[row] = 6;
        c[row] = 1;
        d[row] = 4*(values[row-1] + values[row]);
        v[row] = 0;
    }

    thomas_algorithm(a, b, c, d, v, size+1);
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
}



void evalFx(double* F, double* Fx, int nX, int nY , double step, double* d2Fx)
{
    double* fRow = new double[nX];
    double* tempV = new double[nX+1];
    // printf("memory leak not before Fx\n");
    for(int i = 0; i < nY; i++)
    {
        for(int j = 0; j < nX; j++)
        {
            fRow[j] = F[i*nX + j];
        }
        vCoeff(tempV, fRow, nX, step, d2Fx[i], d2Fx[i+nY]);
        // printf("vCoeff finished\n");
        // print_vect(tempV, nX+1);
        // printf("\n");
        for(int j = 0; j < nX+1; j++)
            Fx[j*nY + i] = tempV[j];
    }
    // print_matr(Fx, nX+1, nY);
    delete [] fRow;
    delete [] tempV; 


}

void evalFy(double* F, double* Fy, int nX, int nY , double step, double* d2Fy)
{
    double* fCol = new double[nY];
    double* tempV = new double[nY+1];

    for(int i = 0; i < nX; i++)
    {
        for(int j = 0; j < nY; j++)
        {
            fCol[j] = F[j*nX + i];
        }
        vCoeff(tempV, fCol, nY, step, d2Fy[i], d2Fy[i+nX]);
        for(int j = 0; j < nY+1; j++)
            Fy[i*(nY+1) + j] = tempV[j];
    }
    // printf("\n");
    // print_matr(Fy, nX, nY+1);
    // printf("\n");
    delete [] fCol;
    delete [] tempV; 


}

void evalFxy(double* Fy, double* Fxy, int nX, int nY , double step, double* d2Fx)
{
    double* fCol = new double[nX];
    double* tempV = new double[nX+1];
    for(int i = 0; i < nY + 1; i++)
    {
        for(int j = 0; j < nX; j++)
            fCol[j] = Fy[j*(nY+1) + i];
        if(i < nX)
            vCoeff(tempV, fCol, nX, step, d2Fx[i], d2Fx[i+nX]); 
        else
            vCoeff(tempV, fCol, nX, step, d2Fx[nX-1], d2Fx[2*nX-1]);

        // print_vect(tempV, nX+1);
        for(int j = 0; j < nX+1; j++)
            Fxy[j*(nY+1) + i] = tempV[j];
    }
    // print_matr(Fxy, nX+1, nY+1);
    delete[] fCol;
    delete[] tempV;
}

void evalFij(double* F, double* Fx, double* Fy, double* Fxy, double* tensorF, int nX, int nY)
{
    for(int j = 0; j < nY; j++)
    {
        for(int i = 0; i < nX; i++)
        {
            tensorF[j*nX*9 + i*9]     = Fxy[i*(nY+1) + j];
            tensorF[j*nX*9 + i*9 + 1] = Fx[i*nY + j];
            tensorF[j*nX*9 + i*9 + 2] = Fxy[i*(nY+1) + j+1];
            tensorF[j*nX*9 + i*9 + 3] = Fy[i*(nY+1) + j];
            tensorF[j*nX*9 + i*9 + 4] = F[j*nX + i];//
            tensorF[j*nX*9 + i*9 + 5] = Fy[i*(nY+1) + j+1];
            tensorF[j*nX*9 + i*9 + 6] = Fxy[(i+1)*(nY+1) + j];
            tensorF[j*nX*9 + i*9 + 7] = Fx[(i+1)*(nY) +j];
            tensorF[j*nX*9 + i*9 + 8] = Fxy[(i+1)*(nY+1) + j+1];

            
        }
    }
}

void basis(double x, double y, double xi, double yj, double* arr)
{
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            arr[i*3 + j] = pow(x - xi,j)*pow(y - yj,i);
        }
    }
}

void interpolatetedVal(double* gamma, int xInt, int yInt, double* intVal, double x, double y, double ksiX, double ksiY, int nX, int nY)
{
    double* basFuncs = new double[9];
    intVal[yInt*nX + xInt] = 0;
    basis(x, y, ksiX, ksiY, basFuncs);
            for(int k = 0; k < 9; k++)
                intVal[yInt*nX + xInt] += gamma[yInt*nX*9 + xInt*9 + k]*basFuncs[k];
    // printf("At (%d, %d) we get %f\n\n", xInt, yInt, intVal[xInt*nY + yInt]);
    if(xInt == 0 && yInt == 0)
    {
        // for(int k = 0; k < 9; k++)
        //     std::cout << gamma[yInt*nX*9 + xInt*9 + k] << "  " << basFuncs[k] << "\n";
    }
    delete [] basFuncs;
}

