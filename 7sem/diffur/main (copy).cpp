#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

#define EPS 1.e-11

void f(double *x, double *k,double h);
void coef_calc(double h, double* x, double* k0,double* k1,double* k2,double* k3,double* k4,double* k5);
double error_calc(double* k0, double* k1, double* k2, double* k3, double* k4, double* k5);
void Runge_Kutta (double x_0);


double  a10 = 1./2., 
        a20 = 1./4.,
        a21 = 1./4.,
        a30 = 0.,
        a31 = -1.,
        a32 = 2.,
        a40 = 7./27.,
        a41 = 10./27.,
        a42 = 0.,
        a43 = 1./27.,
        a50 = 28./625.,
        a51 = -125./625.,
        a52 = 546./625.,
        a53 = 54./625.,
        a54 = -378./625.;

double  r0 = -42./336.,
        r1 = 0.,
        r2 = -224./336.,
        r3 = -21./336.,
        r4 = 162./336.,
        r5 = 125./336.;

double  p0 = 1.0/6.0,
        p1 = 0.,
        p2 = 4.0/6.0,
        p3 = 1.0/6.0,
        p4 = 0.,
        p5 = 0.;

void f(double *x, double *k,double h)
{
    k[0] = h*(x[1]);
    k[1] = -h*(x[0] - 1*(1-x[0]*x[0])*x[1]);
}

void coef_calc(double h, double* x, double* k0,double* k1,double* k2,double* k3,double* k4,double* k5, int n)
{
    double a[n];

    f(x, k0, h);
    for(int i = 0; i < n; i++) {
        a[i] = x[i]+a10*k0[i];
    }    
    f(a,k1,h);
    for(int i = 0; i < n; i++) {
        a[i] = x[i]+a20*k0[i]+a21*k1[i];
    }    
    f(a,k2,h);
    for(int i = 0; i < n; i++) {
        a[i] = x[i]+a30*k0[i]+a31*k1[i]+a32*k2[i];
    }    
    f(a,k3,h);
    for(int i = 0; i < n; i++) {
        a[i] = x[i]+a40*k0[i]+a41*k1[i]+a42*k2[i]+a43*k3[i];
    }    
    f(a,k4,h);
    for(int i = 0; i < n; i++) {
        a[i] = x[i]+a50*k0[i]+a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i];
    }    
    f(a,k5,h);
}
double error_calc(double* k0, double* k1, double* k2, double* k3, double* k4, double* k5, int n)
{
    double e[n];
    double err=0;
    for(int i = 0; i < 2; i++) {
        e[i] = r0*k0[i]+r1*k1[i]+r2*k2[i]+r3*k3[i]+r4*k4[i]+r5*k5[i];
        err=err+e[i]*e[i];
    }  
    return sqrt(err);
}

void Runge_Kutta (double* x_start, double* x_finish, int n, int steps)
{
    double* k0 = new double[n];
    double* k1 = new double[n];
    double* k2 = new double[n];
    double* k3 = new double[n];
    double* k4 = new double[n];
    double* k5 = new double[n];
    double error;
    double error_sum = 0;
    double x[n];
    for(int i=0; i<n; i++){
        x[i]=x_start[i];
    }
    double x_prev[n];
    double h = 0.1;
    int count = 0;
    double dist;
    ofstream fout;
    fout.open("data.txt");
    for (int i=0; i<n; i++){
            fout << x[i];
            if (i != n-1) fout << std::endl;
    }
    fout << std::endl;
    while (true) {
        for(int i=0; i<n; i++) {    
            x_prev[i]=x[i];
        }
        coef_calc(h, x, k0, k1, k2, k3, k4, k5, n);
        error = error_calc(k0, k1, k2, k3, k4, k5, n); 
        if (error < EPS) {
            for (int i = 0; i < n; ++i) {            
                x[i] += p0*k1[i]+p2*k3[i]+p3*k4[i];
            }
            for (int i=0; i<n; i++) {
                fout << x[i];
                if (i != n-1) fout << std::endl;
            }
            fout << std::endl;
            ++count;
            error_sum += error;
        }
        h = h*fmin(2.,fmax(0.5, 0.95*pow(EPS/error, 1/6.)));
        dist = sqrt((x[0] - x_finish[0])*(x[0] - x_finish[0]) + (x[1] - x_finish[1])*(x[1] - x_finish[1]));
        if(count > steps) {
            
            break;
            
        }
    }
    cout << error_sum << ' ' << count << endl;
    fout.close();
    delete [] k0;
    delete [] k1;
    delete [] k2;
    delete [] k3;
    delete [] k4;
    delete [] k5;
}


int main (int argc, char** argv)
{
    int n=2;

    double *x_start = new double[2];
    x_start[0] = atof(argv[1]);
    x_start[1] = atof(argv[2]);
    
    double *x_finish = new double[2];
    x_finish[0]=atof(argv[1]);
    x_finish[1]=atof(argv[2]);
    int steps = atoi(argv[3]);
    Runge_Kutta(x_start, x_finish, n, steps);
    
    return 0;
}