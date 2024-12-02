#ifndef MYINIT
#define MYINIT



#include <cmath>
double fun(int k, double x);
double dFun(int k, double x);
double diffFun2(int k, double x1, double x2);
double diffFun3(int k, double x1, double x2, double x3);

double fun(int k, double x)
{
	switch (k)
	{
	case 0:
		return 1;
		break;
	case 1:
        return x;
		break;
	case 2:
        return x*x;
		break;
	case 3:
		return pow(x, 3);
		break;
	case 4:
		return pow(x, 4);
		break;
	case 5:
		return exp(x);
		break;
	case 6:
		return 1 / (25 * pow(x, 2) + 1);
		break;
	default:
        //std::cout << "INVALID FUNCTION NUMBER";
		return -1;
		break;
	}
}

double d2Fun(int k, double x)
{
	switch (k)
	{
	case 0:
		return 0;
		break;
	case 1:
		return 0;
		break;
	case 2:
		return 2 ;
		break;
	case 3:
		return 6*x;
		break;
	case 4:
		return 12 * pow(x, 2);
		break;
	case 5:
		return exp(x);
		break;
	case 6:
		return (5000*pow(x,2))/pow(25*pow(x,2)+1, 3) - 50/pow(25*pow(x,2)+1, 2);
		break;
	default:
        //std::cout << "INVALID FUNCTION NUMBER";
		return -1;
		break;
	}
}

double diffFun2(int k, double x1, double x2)
{
	double diff = (fun(k, x2) - fun(k, x1)) / (x2 - x1);
	return diff;
}

double diffFun3(int k, double x1, double x2, double x3)
{
	double diff = (diffFun2(k, x2, x3) - diffFun2(k, x1, x2)) / (x3 - x1);
	return diff;
}

void initValues(double *Points, int size, int func, double *values)
{
	for (int i = 0; i < size; ++i)
	{
		values[i] = fun(func, Points[i]);
	}
}

#endif
