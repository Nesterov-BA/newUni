#ifndef QRREFLECT
#define QRREFLECT

#include "functions.hpp"
#include <complex>
#include <cstdio>


void refmatr(double* buff, double* A, double* first, double* second, int n, bool* problemsolved, double* added);

void qrref(double* A, double* buff, double* eigenval, double precision, int n);

void refmatr(double* buff, double* A, double* first, double* second, int n, bool* problemsolved, double* added)
{
	double* refvec = new double[n];
	double* buffvec = new double[n]; 
	double* buff2 = new double[n*n];
	double norm, scalar1, scalar2, d, check, checkid;
	checkid = 0;
	for (int col = 0; col < n; ++col)
	{
		d = n-col;
		//printf("(");	
		for (int i = 0; i < n; ++i)
		{	
			if (i==col || i == col+1)
			{
				refvec[i] = A[col+n*i];
				//printf("%f, ", refvec[i]);
			} else {
				refvec[i] = 0;
				//printf("%f, ", refvec[i]);
			}		
		}

		//printf(")\n \n");
		norm = vecmod(refvec, n);
		refvec[col] -= norm;
		normalize(refvec, n);
		first[col] = refvec[col];
		if(col < n-1)
			second[col] = refvec[col+1];
		if(col == n-1)
			second[col] = 0;
		//printf("%f, %f\n", refvec[col], refvec[col+1]);
		/*printf("x= (");
		for (int i = 0; i < n-1; ++i)
		{
			printf("%f, ", refvec[i]);
		}*/
		//printf("%f), %f \n\n", refvec[n], vecmod(refvec, n));
		for (int i = col; i < n; ++i)
		{
			if(col < n-1)
				scalar1 = 2*A[i+col*n]*refvec[col] + 2*A[i+(col+1)*n]*refvec[col+1];
			if(col == n-1)
				scalar1 = 2*A[i+col*n]*refvec[col];
			A[i+n*col] -= refvec[col]*scalar1;
			if(col < n-1)
				A[i+n*(col+1)] -= refvec[col+1]*scalar1;
			
				
			
			/*for (int j = col; j < n; ++j)
			{
				buffvec[j] = buff[i+n*j];
			}
			scalar2 = 2*scalar_product(buffvec, refvec, n);
		
			for (int j = col; j < n; ++j)
			{
				buff[i+n*j] -= refvec[j]*scalar2; 	
			}	*/
		}



		//check = checkunit(buff, n);
		
		//printf("%f\n", check);
	}
	print_matrix(n, A);
	if (!*problemsolved)
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				if (i!=j)
				{
					checkid += A[i*n+j]*A[i*n+j];
				}
			}
		}
	}
	if (checkid > 0.000001 && !*problemsolved)
	{
		*problemsolved = true;
	}

	for (int row = 0; row < n; ++row)
	{
		for (int i = row; i < n; ++i)
		{
			//printf("%f, %f\n", first[row], second[row]);
			if(row < n-1)
				scalar1 = 2*A[row+i*n]*first[row] + 2*A[row+1+i*n]*second[row];
			if(row == n-1)
				scalar1 = 2*A[row+i*n]*first[row];
			A[row+n*i] -= first[row]*scalar1;
			if(row < n-1)
				A[row+1+n*i] -= second[row]*scalar1;
		}
	}
	print_matrix(n, A);
	if(!*problemsolved)
	{
		for (int i = 0; i < n; ++i)
		{
			A[i*n+i] += 1;
		}
		*added += 1;
	}
	/*for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			buff2[j+i*n] = buff[j+i*n];
		}
	}	
	sqmatr_mult(buff2, A, n, 0);
	printf("Real result:\n");
	print_matrix(n, buff2);
*/}

void qrref(double* A, double* buff, double* eigenval, double precision, int n)
{
	bool solved  = false;
	bool* solvedadr = &solved;
	double diff, d, added;
	double* addedadr = &added;
	printf(" precision = %f\n", precision);
	double* first = new double[n];
	double* second = new double[n];
	//added = 0;
	diff = 10;
	int k = 0;
	refmatr(buff, A, first, second, n, solvedadr, addedadr);
	
	//sqmatr_mult(A, buff, n, 0);
	for (int i = 0; i < n; ++i)
	{
		eigenval[i] = A[i+i*n];
	}
	while(diff>precision)
	{
		id(buff, n);
		refmatr(buff, A, first, second, n, solvedadr, addedadr);
		//sqmatr_mult(A, buff, n, 0);
		diff = 0;
		for (int i = 0; i < n; ++i)
		{	
			d = eigenval[i] - A[i+i*n];
			//sprintf("%f, %f, d= %f\n", eigenval[i], A[i+i*n], d);
			eigenval[i] = A[i+i*n];
			diff += d*d;
			//printf("%f\n", diff);
		}
		diff = sqrt(diff);
		k++;
	}
	if(added > 0.001)
	{
		for (int i = 0; i < n; ++i)
		{
			A[i*n+i] -= added;
		}
	}

	printf("%d\n", k);
}



#endif