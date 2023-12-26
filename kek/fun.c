#include "fun.h"



double egenValue ( const int m, const int N, const double h)
{
	return 4./h/h*sin(M_PI*m/(2.*N - 1)*sin(M_PI*m/(2.*N - 1)));
}

int egenVector( double *y, const int m, const int N)


{
	for( int k=1; k<=N-1; k++)
	{
		y[k] = sin( M_PI*m*2*k/(2.*N-1));
	}
	y[N] = -y[N-1];
	y[0] = 0;

	return 0;
}

double prod( const double *u, const double *v, const int N, const double h)
{
	double s = 0;

	for( int k= 0; k<N; k++)
	{
		s+=u[k]*v[k]*h;
		//printf("s %e\n", s);
	}
	

	return s;
}

double measure( double *A, const double *y, const double lambda, const int N, const double h)
{
	int k;
	double Measure = 0;

	A[0] = lambda*y[0];
	for( k = 1; k<=N-1; k++)
	{
		A[k] = (y[k+1] -2*y[k] + y[k-1] )/h/h + lambda*y[k];
	}

	A[N] =  -( - 2*y[N] + y[N-1])/h/h + lambda*y[N];

	Measure = sqrt(prod(A,A,N,h))/lambda;
	
	return Measure;

}
