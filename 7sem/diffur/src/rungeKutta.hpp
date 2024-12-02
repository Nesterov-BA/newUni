#include <bits/types/FILE.h>
#include <fstream>
#include <iostream>
#include "normCalculation.hpp"

using namespace std;
void clear();
void solutionUpToTime(double xStart, double yStart, double f(double, double), double g(double, double), double finish, double* cycleTimeLess, double* cycleTimeMore,double* xEnd, double* yEnd);
void findCycle(double xStart, double yStart, double f(double, double), double g(double, double), double* xEnd, double* yEnd);
void Runge_Kutta4ClassicSimple(double startX, double startY, double step, double f(double, double), double g(double, double), double* endX, double* endY);
void Runge_Kutta4StepVariedSimple(double startX, double startY, double f(double, double), double g(double, double), double* step, double* endX, double* endY, double* errorSum, double* globalError, double* globalErrorRegular);



void fasterFindCycle(double xStart, double yStart, double f(double, double), double g(double, double), double* xEnd, double* yEnd);
int fasterSolutionUpToTime(double xStart, double yStart, double f(double, double), double g(double, double), double timeEnd, double* cycleTimeLess, double* cycleTimeMore, double* xEnd, double* yEnd);

