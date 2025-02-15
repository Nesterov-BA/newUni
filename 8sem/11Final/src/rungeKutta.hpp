#include "metricFuncs.hpp"
#include "mainhead.hpp"
// #include "normCalculation.hpp"

void solutionUpToTime(double* start, double* end, function* functions, double finish);
void solutionUpToTime(double* start, double* end, function* functions, double finish, string filename, double* integral, double* globalError);
void dorPri5(double* start, double* end, double* err, function* functions, double step, double time);
void dorPri5Varied(double* start, double* end, function* functions, double* step, double time);

void solve(double* start, std::vector<double> (*function)(double, double));
void Runge_Kutta4ClassicSimple(double* start, double* end, function* functions, double step, double time);
void Runge_Kutta4StepVariedSimple(double* start, double* end, function* functions, double* step, double time);

void shooting(double* start, double alpha, function* functions, double finish);

std::vector<double> probe(double* start, std::vector<double> (*function)(double, double));
int findMinimum(double* start, std::vector<double> (*function)(double, double));
void jacobiMatrix(double p1, double p2, std::vector<double> (*func)(double, double), double** matrix);
std::vector<double> error(double p1, double p2, function* functions);
std::vector<double> revError(double p1, double p2, function* functions, double alpha);

void jacobiCheck(function* functions, string filename);
