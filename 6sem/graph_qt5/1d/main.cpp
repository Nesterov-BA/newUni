using namespace std;

#include <iostream>
#include "valInit.hpp"
#include "coefficents.hpp"

using namespace std;
void initValues(double* Points, int size, int func, double* values);
void coefficents_eval(double* Points, double* values, int size, int func, double* coeff);


int main(int argc, char** argv){
	double left, right, step;
	int pointsCount, funcNum; 
	
	if (argc != 5){
		cout << "INVALID NUMBER OF ARGUMENTS";
		return -1;
	}
	
	left = strtod(argv[1], NULL);
	right = strtod(argv[2], NULL);
	
	if (!(right-left>0)){
		cout << "INVALID INTERVAL BORDER";
		return -2;
	}
	
	pointsCount = atoi(argv[3]);
	
	if (pointsCount < 2) {
		cout << "INVALID NUMBER OF INTERPOLATION POINTS";
		return -3;		
	}

	funcNum = atoi(argv[4]);
	step = (right-left)/(pointsCount-1);
	
	double* points = new double[pointsCount];
	double* values = new double[pointsCount]; 
	points[0] = left;
	for (int i = 1; i < pointsCount; ++i) {
		points[i] = points[i-1] + step; 
	}

	initValues(points, pointsCount, funcNum, values);

	double* coefficents = new double[4*pointsCount];
	
	coefficents_eval(points, values, pointsCount, funcNum, coefficents);
	
	for (int i = 0; i < pointsCount; i++)
	{
		cout << points[i] << " " << values[i] << "\n";
	}
	
	return 0;

}
