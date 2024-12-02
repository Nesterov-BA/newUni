#include <stdio.h>
#include <iostream>
#include <time.h>
#include "getmatrix.hpp"
#include "myalgo.hpp"
#include "printmatrix.hpp"
#include "mistakematrix.hpp"
#include "functions.hpp"
#include "matrspin.hpp"
#include "qrreflect.hpp"

void matrix_init(int n,  double * initmat); 
int get_matrix(int n, const char* filename, double* A);
void print_matrix(int n, const double* A); 
int my_algo(int n, double* A, double* ans);
double func(int k, int n, int i, int j); 
void calc_mistake(int n, double* A, double* ans, double* temp); 

int main(int argc, char **argv)
{
    std::cerr << argc << "\n";
    if ((argc != 2) && (argc != 3)) {return -1;} 
    int n = atoi(argv[1]); 

    std::cerr << n << "\n";

    int res = 0;
    double *A = new double[n * n];
    double *temp = new double[n * n];
    if (argc == 2) {
        matrix_init(n, A);
        matrix_init(n, temp);
    } else {
        res = get_matrix(n, argv[2], A);
        res = get_matrix(n, argv[2], temp);
        if (res != 0)  {return -1;}
    }

    //print_matrix(n, A);
    
    double* ans = new double[n * n];
    double* eigenval = new double[n];
    double precision = 0.000001;
    id(ans, n);

    clock_t tStart = clock();
    triang(A, n);
    print_matrix(n, A);
    qrref(A, ans, eigenval, precision, n);
    clock_t tEnd = clock();


    printf("Time taken: %.2fs\n", (double)(tEnd - tStart)/CLOCKS_PER_SEC);

    print_matrix(n, A);
    print_matrix(n, ans);

   
    double * temp2 = new double[n * n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            temp2[i*n + j] = 0;
        }
    }

    calc_mistake(n, temp, ans, temp2);
    



    delete[] A;
    delete[] ans;
    delete[] temp;
    delete[] temp2;
 
    return 0;
}