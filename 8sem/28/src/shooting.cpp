#include "shooting.hpp"
#include "rungeKutta.hpp"
#include "gaussian.hpp"
#include <cmath>
#include <cstdio>
#include <iterator>
#include <vector>


void reverseJacobiMatrix(double p1, double x2, std::vector<double> (*func)(double, double), double** matrix)
{
    double eps = 1.e-9;
    std::vector<double> res(2);
    std::vector<double> resX(2);
    std::vector<double> resY(2);
    res = func(p1, x2);
    resX = func(p1 - eps, x2);
    resY = func(p1, x2 - eps);
    resX[0] -= res[0];
    resX[0] /= eps;
    matrix[0][0] = -resX[0];
    resX[1] -= res[1];
    resX[1] /= eps;
    matrix[0][1] = -resX[1];
    resY[0] -= res[0];
    resY[0] /= eps;
    matrix[1][0] = -resY[0];
    resY[1] -= res[1];
    resY[1] /= eps;
    matrix[1][1] = -resY[1];
//    printf("At %lf, %lf values are %lf, %lf\nThe reverse jacobi matrix is:\n", p1, x2, res[0], res[1]);
 //   printMatrix(matrix, 2);
}
void jacobiMatrix(double p1, double x2, std::vector<double> (*func)(double, double), double** matrix)
{
    double eps = 1.e-9;
    std::vector<double> res(2);
    std::vector<double> resX(2);
    std::vector<double> resY(2);
    res = func(p1, x2);
    resX = func(p1 + eps, x2);
    resY = func(p1, x2 + eps);
    resX[0] -= res[0];
    resX[0] /= eps;
    matrix[0][0] = resX[0];
    resX[1] -= res[1];
    resX[1] /= eps;
    matrix[0][1] = resX[1];
    resY[0] -= res[0];
    resY[0] /= eps;
    matrix[1][0] = resY[0];
    resY[1] -= res[1];
    resY[1] /= eps;
    matrix[1][1] = resY[1];
  //  printf("At %lf, %lf values are %lf, %lf\nThe jacobi matrix is:\n", p1, x2, res[0], res[1]);
   // printMatrix(matrix, 2);
}

void probe(double* start, std::vector<double> (*function)(double, double))
{
    double alphaStep = 0.001;
    alpha = 0;
    printf("Alpha: 0, start at: %lf, %lf\n", start[0], start[1]);
    while (alpha < 5.1)
    {
        alpha += alphaStep;
        printf("Alpha: %lf, start at: %lf, %lf\n", alpha, start[0], start[1]);
        if(findMinimum(start, function) < 0)
            break;
        if(fabs(alpha - 0.1) < 1.e-11)
            printf("Alpha: 0.1, start at: %lf, %lf\n", start[0], start[1]);
        if(fabs(alpha - 1) < 1.e-11)
            printf("Alpha: 1.0, start at: %lf, %lf\n", start[0], start[1]);
    }
    printf("Alpha: 5.1, start at: %lf, %lf\n", start[0], start[1]);
}

int findMinimum(double* start, std::vector<double> (*function)(double, double))
{
    double correction[2];
    double** jacobi = new double*[2];
    jacobi[0] = new double[2];
    jacobi[1] = new double[2];
    double coefficent = -1;
    std::vector<double> errorVector(2);
    double currNorm;
    int count = 0;
    while(true)
    {
        count++;
        coefficent = 1;
        errorVector = function(start[0], start[1]);
        //printf("Error vector: %lf, %lf\n", errorVector[0], errorVector[1]);
        if(l2Norm(errorVector, 2) < 1.e-7)
        {
            printf("Small Error!\n");
            printf("res = %.e, %.e\n", errorVector[0], errorVector[1]);
            return 0;
            break;
        }
        correction[0] = errorVector[0];
        correction[1] = errorVector[1];
        currNorm = l2Norm(errorVector, 2);
      //  reverseJacobiMatrix(start[0], start[1], function, jacobi);
        jacobiMatrix(start[0], start[1], function, jacobi);
        gauss(jacobi, correction);
        std::vector<double> res = function(start[0] + correction[0]*coefficent, start[1] + correction[1]*coefficent);
    //    printf("Correction vector: %lf, %lf\n", correction[0], correction[1]);
        int count = 0;
        while(l2Norm(function(start[0] + correction[0]*coefficent, start[1] + correction[1]*coefficent), 2) > currNorm)
        {
            coefficent /= 2;
            count++;
            res = function(start[0] + correction[0]*coefficent, start[1] + correction[1]*coefficent);
     //       printf("Coeff: 2^(-%d), Corrected = (%lf, %lf)\n", count, res[0], res[1]);
            if(fabs(coefficent) < 1.e-9)
            {
                printf("Small Coefficent\n");
                printf("Res = (%lf, %lf)\n", errorVector[0], errorVector[1]);
                return -1;
            }
        }

        start[0] += correction[0]*coefficent;
        start[1] += correction[1]*coefficent;
        if(l2Norm(function(start[0], start[1]), 2) < 1.e-7)
        {
            printf("Small Error!\n");
            printf("res = %.e, %.e\n", function(start[0], start[1])[0], function(start[0], start[1])[1]);
            return 0;
            break;
        }
        if(count > 150)
        {
            printf("Big Count!\n");
            printf("Res = (%lf, %lf)\n", res[0], res[1]);
            return -1;
            break;
        }
    }
}

std::vector<double> error(double p2, double x1, function* functions)
{
    double start[4] = {0, p2, x1, 0};
    double end[4];
    std::vector<double> res(2);
    double finish = 1;
    solutionUpToTime(start, end, functions, finish);
    res[0] = end[1];
    res[1] = end[2];
    return res;
}
