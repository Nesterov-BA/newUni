#include "metricFuncs.hpp"
#include "mainhead.hpp"
// #include "normCalculation.hpp"

void solutionUpToTime(double* start, double* end, function* functions, double finish);
void solutionUpToTime(double* start, double* end, function* functions, double finish, double* globalError);

void Runge_Kutta4ClassicSimple(double* start, double* end, function* functions, double step, double time);
void Runge_Kutta4StepVariedSimple(double* start, double* end, function* functions, double* step, double time, double* globalError);


