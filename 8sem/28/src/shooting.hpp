#include <vector>
#include "mainhead.hpp"
void reverseJacobiMatrix(double p1, double x2, std::vector<double> (*func)(double, double), double** matrix);
void jacobiMatrix(double p1, double p2, std::vector<double> (*func)(double, double), double** matrix);
void probe(double* start, std::vector<double> (*function)(double, double));
int findMinimum(double* start, std::vector<double> (*function)(double, double));
std::vector<double> error(double p2, double x1, function* functions);
