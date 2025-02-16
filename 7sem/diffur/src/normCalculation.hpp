#include <math.h>
#include <cmath>

extern double alpha;

double* solveQuadratic(double a, double b, double c);
double logNormCalc(double x1, double y1, double x2, double y2);
double regNormCalc(double x1, double y1, double x2, double y2);
double max4(double x1, double x2, double x3, double x4);
inline double* calculateEigen2by2(double* matrix);
