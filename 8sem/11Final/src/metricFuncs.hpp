#include "math.h"
#include <vector>


double l1Norm(double* vector, int size);
double lInfNorm(double* vector, int size);
double l2Norm(std::vector<double> vector, int size);
double l1Metric(double* vector1, double vector2, int size);
void normalize(std::vector<double> vector, int size);

void printMatrix(double** matrix, int size);
