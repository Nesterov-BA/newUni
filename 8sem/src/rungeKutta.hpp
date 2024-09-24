#include <bits/types/FILE.h>
#include <fstream>
#include <iostream>
#include "metricFuncs.hpp"
// #include "normCalculation.hpp"

inline double tolerance;
typedef double (*function)(double, double, double, double);

using namespace std;
void solutionUpToTime(double* start, double* end, function* functions, double finish);
void Runge_Kutta4ClassicSimple(double* start, double* end, function* functions, double step);
void Runge_Kutta4StepVariedSimple(double* start, double* end, function* functions, double step);

void shooting(double* start, double alpha, function* functions, double finish);

double error(double p1, double p2, function* functions);






