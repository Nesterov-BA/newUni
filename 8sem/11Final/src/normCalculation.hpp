#include <math.h>
#include <cmath>
#include <vector>

double logNorm(double time);
double polynom(double t, double lambda);
double dPolynom(double t, double lambda);
double deltaJacNorm(std::vector<double> (*func)(double, double), double p1, double p2, double globErr);
