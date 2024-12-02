#include "normCalculation.hpp"

double* logMatrix = new double[4];
double* regMatrix = new double[4];

double logNormCalc(double x1, double y1, double x2, double y2)
{
    double logNorm = 0;
    double* lambda1 = new double[2];
    double* lambda2 = new double[2];
    logMatrix[0] = 0;
    logMatrix[1] = -alpha*x1*y1;
    logMatrix[2] = logMatrix[1];
    logMatrix[3] = alpha*(1-x1*x1);
    lambda1 = calculateEigen2by2(logMatrix);
    logMatrix[0] = 0;
    logMatrix[1] = -alpha*x2*y2;
    logMatrix[2] = logMatrix[1];
    logMatrix[3] = alpha*(1-x2*x2);
    lambda2 = calculateEigen2by2(logMatrix);
    logNorm = max4(lambda1[0], lambda1[1], lambda2[0], lambda2[1]);
    return logNorm;
}

double regNormCalc(double x1, double y1, double x2, double y2)
{
    double regNorm = 0;
    double* lambda1 = new double[2];
    double* lambda2 = new double[2];
    regMatrix[0] = 0.5;
    regMatrix[1] = alpha*(1-x1*x1);
    regMatrix[2] = regMatrix[1];
    regMatrix[3] = alpha*(1-x1*x1)*alpha*(1-x1*x1) + (-2*alpha*x1*y1-1)*(-2*alpha*x1*y1-1);
    lambda1 = calculateEigen2by2(regMatrix);
    regMatrix[0] = 0.5;
    regMatrix[1] = alpha*(1-x2*x2);
    regMatrix[2] = regMatrix[1];
    regMatrix[3] = alpha*(1-x2*x2)*alpha*(1-x2*x2) + (-2*alpha*x2*y2-1)*(-2*alpha*x2*y2-1);
    lambda2 = calculateEigen2by2(regMatrix);
    regNorm = max4(lambda1[0], lambda1[1], lambda2[0], lambda2[1]);
    return sqrt(regNorm);

}
//solve ax^2+bx+c=0
double* solveQuadratic(double a, double b, double c)
{
    double* x = new double[2];
    double D = b*b - 4*a*c;
    x[0] = (-b + sqrt(D))/(2*a);
    x[1] = (-b - sqrt(D))/(2*a);
    // double bigRoot = std::max(x[0], x[1]);
    // x[0] = c/bigRoot;
    // x[1] = -b-x[0];
    return x;
}
//eigenvalues of a 2x2 matrix
double* calculateEigen2by2(double* matrix)
{
    double* eigenvector = new double[2];
    double a = 1;
    double b = -(matrix[0] + matrix[3]);
    double c = matrix[0] * matrix[3] - matrix[1] * matrix[2];
    eigenvector = solveQuadratic(a, b, c);
    return eigenvector;
}

double max4(double x1, double x2, double x3, double x4)
{
    return std::fmax(std::fmax(fabs(x1), fabs(x2)), std::fmax(fabs(x3), fabs(x4)));
}
