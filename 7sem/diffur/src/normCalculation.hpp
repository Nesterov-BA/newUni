#include <math.h>
#include <cmath>

double * solveQuadratic(double a, double b, double c);
double logNormCalc(double x1, double y1, double x2, double y2, double step);
double regNormCalc(double x1, double y1, double x2, double y2, double step);

double alpha = 10;


double logNormCalc(double x1, double y1, double x2, double y2, double step)
{
    double logNorm = 0;
    double* lambda1 = new double[2];
    double* lambda2 = new double[2];
    double b1 = -alpha*(1-x1*x1)/2;
    double b2 = -alpha*(1-x2*x2)/2;
    double c1 = alpha*x1*y1;
    c1 = -c1*c1;
    double c2 = alpha*x2*y2;
    c2 = -c2*c2;
    lambda1 = solveQuadratic(1, b1, c1);
    lambda2 = solveQuadratic(1, b2, c2);
    logNorm = step*std::max(std::max(fabs(lambda1[0]), fabs(lambda1[1])),std::max(fabs(lambda2[0]), fabs(lambda2[1])));
    return logNorm;

}

double regNormCalc(double x1, double y1, double x2, double y2, double step)
{
    double norm = 0;
    double a1 = 1; 
    double b1 = alpha*(1-x1*x1);
    double c1 = alpha*(1-x1*x1)*alpha*(1-x1*x1) + (2*alpha*x1*y1 + 1)*(2*alpha*x1*y1 + 1);
    double a2 = 1; 
    double b2 = alpha*(1-x2*x2);
    double c2 = alpha*(1-x2*x2)*alpha*(1-x2*x2) + (2*alpha*x2*y2 + 1)*(2*alpha*x2*y2 + 1);
    double D1 = (a1+c1)*(a1+c1) + 4*b1*b1;
    double D2 = (a2+c2)*(a2+c2) + 4*b2*b2;
    double lambda11 = (a1+c1 + sqrt(D1))/2;
    double lambda21 = (a1+c1 - sqrt(D1))/2;
    double lambda12 = (a2+c2 + sqrt(D2))/2;
    double lambda22 = (a2+c2 - sqrt(D2))/2;

    norm = std::max(std::max(fabs(lambda11), fabs(lambda12)), std::max(fabs(lambda21), fabs(lambda22)));
    norm = step*sqrt(norm);
    return norm;

}

double* solveQuadratic(double a, double b, double c)
{
    double* x = new double[2];
    double D = b*b - 4*a*c;
    if (D > 0)
    {
        x[0] = (-b + sqrt(D))/(2*a);
        x[1] = (-b - sqrt(D))/(2*a);
    }
    else if (D == 0)
    {
        x[0] = -b/(2*a);
        x[1] = -b/(2*a);
    }
    return x;
}