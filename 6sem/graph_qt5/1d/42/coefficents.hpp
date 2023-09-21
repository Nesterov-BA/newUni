#ifndef MYCOEFF
#define MYCOEFF

#include <cmath>
#include "valInit.hpp"

using namespace std;

void thomas_algorithm(double* a,
                      double* b,
                      double* c,
                      double* d,
                      double* f, int N) {

  // Create the temporary vectors
  // Note that this is inefficient as it is possible to call
  // this function many times. A better implementation would
  // pass these temporary matrices by non-const reference to
  // save excess allocation and deallocation
  double* c_star = new double[N];
  double* d_star = new double[N];

  // This updates the coefficients in the first row
  // Note that we should be checking for division by zero here
 /* c_star[0] = c[0] / b[0];
  d_star[0] = d[0] / b[0];

  // Create the c_star and d_star coefficients in the forward sweep
  for (int i=1; i<N; i++) {
    double m = 1.0 / (b[i] - a[i-1] * c_star[i-1]);
    if(i<N-1)
        c_star[i] = c[i] * m;
    d_star[i] = (d[i] - a[i-1] * d_star[i-1]) * m;
  }
  c_star[N-1] = 0;
  // This is the reverse sweep, used to update the solution vector f
  for (int i=N-1; i-- > 0; ) {
    f[i] = d_star[i] - c_star[i] * d[i+1];*/
  c[0] /= b[0];
  d[0] /= b[0];
  double m;
  for(int row = 1; row < N; row++)
  {
    b[row] -= c[row-1]*a[row-1];
    d[row] -= d[row-1]*a[row-1];
    c[row] /= b[row];
    d[row] /= b[row];
    b[row] = 1;
  }
    f[N-1] = d[N-1];
    for(int i = 1; i < N; ++i)
    {
        f[N-1-i] = d[N-1-i] - c[N-1-i]*f[N-i];
    }

}

void coefficents_eval(double* Points, double* values, int size, int func, double* coeff, double step)
{
    double* v = new double[size];
    double* a = new double[size-1];
    double* b = new double[size];
    double* c = new double[size-1];
    double* d = new double[size];

    b[0] = 4;
    c[0] = b[0];
    a[size-2] =b[0];
    b[size-1] = b[0];
    d[0] = values[0]*4 + step*step*d2Fun(func, Points[0]);
    d[size-1] = values[size-1]*4 + step*step*d2Fun(func, Points[size-1]);
    for(int row = 1; row < size-1; row++)
    {
        a[row-1] =  1;
        b[row] = 6;
        c[row] = 1;
        d[row] = 4*(values[row-1] + values[row]);
        v[row] = 0;
    }

    thomas_algorithm(a, b, c, d, v, size);

    for (int i = 0; i < size; i++)
    {
        coeff[3*i] = v[i];
        coeff[3*i + 2] = 2*(v[i+1] + v[i] - 2*values[i])/(step*step);
        coeff[3*i + 1] = 2*(values[i] - v[i])/step - step*coeff[3*i +2]/2;
    }
}

double approx_eval(double point, double* coeff, int intervalNum, double* points)
{
    double res = coeff[3*intervalNum] + coeff[3*intervalNum+1]*(point - points[intervalNum]) + coeff[3*intervalNum+2]*pow(point - points[intervalNum],2);
    return res;
}


#endif
