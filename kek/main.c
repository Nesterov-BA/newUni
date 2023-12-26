#include "fun.h"

int main(void)
{
	int N; 
	int i;
	int m1, m2;
	int m1_max, m2_max;
	double h;
	double mera, mera_max;
	int m, m_max;
	double *y1, *y2;
	double s, smax;
	double *A;
	double *lambda;

	
	for(i=2; i<=10; i++)
	{
		N=pow(2, i);

		y1=(double *)malloc((N + 1) * sizeof(double));
		if(!y1)
		{
			fprintf(stderr, "[X] Memory allocation error...\n");
			return -1;
		}

		y2=(double *)malloc((N + 1) * sizeof(double));
		if(!y2)
		{
			fprintf(stderr, "[X] Memory allocation error...\n");
			free(y1);
			return -1;
		}


		A=(double *)malloc((N + 1) * sizeof(double));
		if(!A)
		{
			fprintf(stderr, "[X] Memory allocation error...\n");
			free(y1);
			free(y2);
			return -1;
		}

		lambda=(double *)malloc(N * sizeof(double));
		if(!lambda)
		{
			fprintf(stderr, "[X] Memmory allocation error\n");
			free(y1);
			free(y2);
			free(A);
			return -1;
		}

		h=1./(double)N-0.5;

		smax=0.0;
		m1_max=0;
		m2_max=0;
	
		for(m1 = 1; m1 < N; m1++)
		{
			for(m2=m1+1; m2 < N; m2++)
			{
				egenVector(y1, m1, N);
				egenVector(y2, m2, N);
				s = prod(y1, y2, N, h);
				if(s > smax)
				{
					m1_max = m1;
					m2_max = m2;
					smax=s;
				}

			}
		}
	
		
		printf("Size 2^%d\n", i);
		printf("max(y(m) , y(n))h: %e \n", smax);
		printf("NumVector: %d and %d\n\n", m1_max, m2_max);
		mera_max=0.0;
		m_max=0;

		for(m = 1; m < N; m++)
		{
			lambda[m] = egenValue(m, N, h);
			egenVector(y1, m, N);
			mera = measure(A, y1, lambda[m], N, h);
			
			if( mera > mera_max)
			{
				mera_max = mera;
				m_max=m;
			}
		}
	
		printf("max||Ay(m) - Lambda(m)*y(m) || / Lambda(m): %e \n", mera_max);
		printf("Max m: %d\n\n\n", m_max);
	
	
		free(y1);
		free(y2);
		free(A);
		free(lambda);
	}
	printf("egenVector :: y[k] = sin( M_PI*m*2*k/(2.*N-1))\n\n");
	printf("egenValue :: Lambda[k] = 4./h/h*sin(M_PI*m/(2.*N - 1)*sin(M_PI*m/(2.*N - 1)))\n\n");
	return 0;
		

}

