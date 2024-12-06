#include "math.h"
#include <iostream>

void gradient(double func(double, double), double* start, double* gradient);
void newton2d(double func(double, double), double* start);
double derivative(double func(double, double), double* start, double* direction);
void gradient_descent(double func(double, double), double* start);
