#ifndef MYALGO
#define MYALGO
#include "functions.hpp"

void swaprows(double* matr, int size, int i, int j);
void rowtimesx(double* matr, int size, int i, double x);
void rowplusrowtimesx(double* matr, int size, int i, int j, double x);

void swaprows(double* matr, int size, int i, int j){
	double buffer;
	for(int k = 0; k<size; k++){
		buffer = matr[i*size+k];
		matr[i*size+k] = matr[j*size+k];
		matr[j*size+k] = buffer;
	}
}

void rowtimesx(double* matr, int size, int i, double x){
	for(int k = 0; k<size; k++){
		matr[i*size+k] *= x;
	}
}

void rowplusrowtimesx(double* matr, int size, int i, int j, double x){
	for(int k = 0; k<size; k++){
		matr[i*size+k] += x*matr[j*size+k];
	}	
} 



int my_algo(int n, double* A, double* ans){
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			if (i==j) {
				ans[j+n*i] = 1;
			} else {
				ans[j+n*i] = 0;
			}

		}
	}
	for(int i=0; i<n; i++){
		if (dabs(A[i*n +i]) > 0) {
			rowtimesx(ans, n, i, 1/A[i*n+i]);
			rowtimesx(A, n, i, 1/A[i*n+i]);
		} else {
			std::cout << "incompatible matrix\n";
			return -1;
		}
		for (int j = 0; j < n; j++){
			if (j!=i){
				rowplusrowtimesx(ans, n, j, i, -A[j*n+i]);
				rowplusrowtimesx(A, n, j, i, -A[j*n+i]);
			}	
		}
	};
	return 0;
}
#endif