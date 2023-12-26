#include "fun.h"



double egenValue ( const int m, const int N, const double h)
{
	return 4./h/h*sin(M_PI*m/(2.*N))*sin(M_PI*m/(2.*N));
}

int egenVector( double *y, const int m, const int N)
{
	for( int k=0; k<=N; k++)
	{
		y[k] = sin( M_PI*m*k/N);
	}

	return 0;
}

double prod( const double *u, const double *v, const int N, const double h)
{
	double s;

	for( int k= 0; k<=N; k++)
	{
		s+=u[k]*v[k]*h;
	}

	return s;
}

double measure( double *A, const double *y, const double lambda, const int N, const double h)
{
	int k;
	double Measure;

	A[0] = 0;
	for( k = 1; k<=N-1; k++)
	{
		A[k] = (y[k+1] -2*y[k] + y[k-1] )/h/h + lambda*y[k];
	}

	A[N] = 0;

	Measure = prod(A,A,N,h);
	Measure = sqrt(Measure);
	Measure=Measure/lambda;
	return Measure;

}
