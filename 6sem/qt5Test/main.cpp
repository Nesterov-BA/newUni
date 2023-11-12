#include "help.hpp"
#include "chebyshev.hpp"
#include <cstdio>
#include <cstdlib>
#include <math.h>

int main(int argc, const char *argv[]){
    double a,b,c,d;
    int nx, ny;
    double stepX, stepY;
    int k, n;
    int nInt;
    
    a = -0.25;
    c = -0.25;
    b = 0.25;
    d = 0.25;
    if(fabs(a-b)<1e-6 || fabs(c-d)<1e-6)
        return -2;
    //fclose(ff);
    if(argv[1] != NULL)
        nx=atoi(argv[1]);
    // if(!ok)
        // return -4;
    // if(nx<5)
    //     return -5;
    if(argv[2] != NULL)
        ny=atoi(argv[2]);
    // if(!ok)
        // return -6;
    // if(ny<5)
    //     return -7;
    if(argv[3] != NULL) 
        k=atoi(argv[3]);
    // if(!ok)
        // return -8;
    if(k<0 || k>7)
        return -9;
    stepX =(b-a)/nx;
    stepY = (d-c)/ny;
    nInt = 3;

    double max_z, min_z, absmax;

    double(*f)(double, double); //function
    double(*d2f)(double, double, int); //

    double* cx = new double[nx];
    double* cy = new double[ny];
    double* F = new double[nx*ny];
    double* intValues = new double[((nx-1)*nInt+1)*((ny-1)*nInt+1)];
    double* intX = new double[((nx-1)*nInt+1)];
    double* intY = new double[((ny-1)*nInt+1)];
    double* realValues = new double[((nx-1)*nInt+1)*((ny-1)*nInt+1)];
    double* Fx = new double[nx*ny];
    double* Fy = new double[nx*ny];
    double* Fxy = new double[nx*ny];
    double* TF = new double[nx*ny*16];
    double* Gamma = new double[nx*ny*16];
    double* d2Fx = new double[2*ny];
    double* d2Fy = new double[2*nx];


    f = f0;
    d2f = d2f0;

    points(cx, cy, nx, ny, a, b, c, d);
    points(intX, intY, (nx-1)*nInt+1, (ny-1)*nInt+1, a, b, c, d);
    // for(int i = 0; i < ny; i++)
    // {
    //     for (int j = 0; j < nx; j++) 
    //     {
    //         printf("(%f, %f) ", cx[j], cy[i]);
    //     }
    //     printf("\n");
    // }
    Fill_F(F, nx, ny, cx, cy, f);
    printf("F[] = %f\n", F[(ny-1)*nx + nx-1]);
    // print_matrix(F, nx, ny);
    Fill_F(realValues, (nx-1)*nInt+1, (nx-1)*nInt+1, intX, intY, f);
    max_z=max_matr(F,nx,ny);
    min_z=min_matr(F,nx,ny);
    absmax=max(fabs(max_z),fabs(min_z));
    Fill_d2F(d2Fx, d2Fy, nx, ny, cx, cy, d2f);
    evalFx(F, Fx, nx, ny, stepX, d2Fy);
    //  print_matrix(Fx, nx, ny);
    evalFy(F, Fy, nx, ny, stepY, d2Fx);
    //  print_matrix(Fy, nx, ny);
    evalFxy(F, Fxy, nx, ny, stepX, d2Fy);
    // print_matrix(Fxy, nx, ny);
    evalFij(F, Fx, Fy, Fxy, TF, nx, ny);
      print_matrix(TF, nx*4, ny*4);
    evalGamma(Gamma, stepX, stepY, TF, nx, ny);
    

    for(int i = 0; i < nx-1; i++)
    {
        for(int j = 0; j < ny-1; j++)
        {
            interpolatetedVal(Gamma, i, j, intValues, cx[i], cx[i+1], cy[j], cy[j+1], nInt, nx, ny);
        }
    }
    printf("bigSize = %d\n", ((nx-1)*nInt+1)*((ny-1)*nInt+1));
    print_matrix(intValues, (nx-1)*nInt+1, (ny-1)*nInt+1);
    double max = 0;
    for(int i = 0; i < ((nx-1)*nInt+1)*((ny-1)*nInt+1); i++)
    {
         if(intValues[i] - realValues[i] > max)
             max = intValues[i] - realValues[i];        
        printf("%f, i = %d\n", intValues[i] - realValues[i], i);
    }
    
    printf("max = %f\n", max);
    delete[] cx ;
    delete[] cy ;
    delete[] F ;
    delete[] intValues;
    delete[] intX  ;
    delete[] intY;
    delete[] realValues;
    delete[] Fx ;
    delete[] Fy;
    delete[] Fxy;
    delete[] TF;
    delete[] Gamma;
    delete[] d2Fx;
    delete[] d2Fy;
    
    return 0;
}
