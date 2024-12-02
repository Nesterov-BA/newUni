#include <bits/types/FILE.h>
#include <fstream>
#include <vector>
#include <iostream>
#include "metricFuncs.hpp"
// #include "normCalculation.hpp"

extern double tolerance;
extern double alpha;
extern double finish;
typedef double (*function)(double, double, double, double);

using namespace std;
void solutionUpToTime(double* start, double* end, function* functions, double finish);
void solutionReverse(double* start, double* end, function* functions, double finish);
void solutionReverseNoErrors(double* start, double* end, function* functions, double time);
void Runge_Kutta4ClassicSimple(double* start, double* end, function* functions, double step);
void Runge_Kutta4StepVariedSimple(double* start, double* end, function* functions, double step);

void shooting(double* start, double alpha, function* functions, double finish);

void probe(double* start, std::vector<double> (*function)(double, double));
int findMinimum(double* start, std::vector<double> (*function)(double, double));
void jacobiMatrix(double p1, double p2, std::vector<double> (*func)(double, double), double** matrix);
std::vector<double> error(double p1, double p2, function* functions, double alpha);
std::vector<double> revError(double p1, double p2, function* functions, double alpha);
