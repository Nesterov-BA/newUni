#ifndef chebyshev_hpp
#define chebyshev_hpp
#include <cmath>
#include <cstdio>

void points(double *cx, double *cy, int nx, int ny, double a, double b, double c, double d);
void Fill_F(double *F, int nx, int ny, double *cx, double *cy, double (*f)(double, double));
void Fill_d2F(double *d2FX, double *d2FY, int nx, int ny,  double *cx, double *cy, double (*d2f)(double, double, int));
void evalGamma(double* gamma, double stepX, double stepY, double* F, int nX, int nY);
void dCoeff(double* d, double* values, int size, double step, double d2ValBeg, double d2ValEnd);
void evalFx(double* F, double* Fx, int nX, int nY , double step, double* d2Fx);
void evalFy(double* F, double* Fx, int nX, int nY , double step, double* d2FY);
void evalFxy(double* F, double* Fx, int nX, int nY , double step, double* d2Fx);
void evalFij(double* F, double* Fx, double* Fy, double* Fxy, double* tensorF, int nX, int nY);
void basis(double x, double y, double xi, double yj, double* arr);
void interpolatetedVal(double* gamma, int xInt, int yInt, double* intVal, double x, double y, double ksiX, double ksiY, int nX, int nY);

#endif /* chebyshev_hpp */
