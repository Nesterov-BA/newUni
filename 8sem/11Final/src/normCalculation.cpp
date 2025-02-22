#include "normCalculation.hpp"
#include <fstream>

double logNorm(double time)
{
    double root = (1 + sqrt(1 + 4*(time-1)*(time-1)))/2;
    return root/2;
}

double polynom(double t, double lambda)
{
    return pow(lambda, 4) + 2*t*pow(lambda, 3) + (-t*t + 2*t-3)*lambda*lambda + (-2*t*t*t + 4*t*t-2*t)*lambda + t*t - 2*t + 1;
}
double dPolynom(double t, double lambda)
{
    return 4*pow(lambda, 3) + 6*t*pow(lambda, 2) + 2*(-t*t + 2*t-3)*lambda + (-2*t*t*t + 4*t*t-2*t);
}

double deltaJacNorm(std::vector<double> (*func)(double, double), double p1, double p2, double globErr)
{
    double norm = 0;
    double secDeriv = 0;
    double eps = 1.e-9;
    secDeriv += fabs(func(p1 + eps, p2)[0] - func(p1,p2)[0] + func(p1 - eps,p2)[0])/(2*eps);
    secDeriv += fabs(func(p1 + eps, p2)[1] - func(p1,p2)[1] + func(p1 - eps,p2)[1])/(2*eps);
    secDeriv += fabs(func(p1, p2 + eps)[0] - func(p1,p2)[0] + func(p1 - eps,p2)[0])/(2*eps);
    secDeriv += fabs(func(p1, p2 + eps)[1] - func(p1,p2)[1] + func(p1 - eps,p2)[1])/(2*eps);
    norm = secDeriv + 12*globErr/eps;
    return norm;
}
