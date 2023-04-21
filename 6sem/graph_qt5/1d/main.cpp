using namespace std;

#include <iostream>
#include "valInit.hpp"

using namespace std;
void initValues(double* Points, int size, int func, double* values);

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

	double* d = new double[pointsCount];
	double* coefficents = new double[4*pointsCount];

	for (int i = 1; i < pointsCount-1; i++)
	{
		d[i] = ((points[i+1] - points[i])*diffFun2(funcNum, points[i-1], points[i]) - (points[i] - points[i-1])*diffFun2(funcNum, points[i], points[i+1]))/(points[i+1]-points[i-1]);
	}
	d[0] = (3*diffFun2(funcNum, points[0], points[1]) - d2Fun(funcNum, points[0])*(points[1]-points[0])/2 - d[1])/2;
	d[pointsCount-1] = (3*diffFun2(funcNum, points[pointsCount-2], points[pointsCount-1]) + d2Fun(funcNum, points[pointsCount-1])*(points[pointsCount-1]-points[pointsCount-2])/2 - d[pointsCount-2])/2;
	// перекинуть вычисление d и коэффицентов в отдельный файл

	
	for (int i = 0; i < pointsCount; i++)
	{
		cout << points[i] << " " << values[i] << "\n";
	}
	
	return 0;

}
