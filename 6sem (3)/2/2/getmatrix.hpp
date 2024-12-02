#ifndef MYMATRX
#define MYMATRX

#include "funcmatrix.hpp"


void matrix_init(int n, double* initmat)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            initmat[i*n + j] = func(n, i, j);
        }
    }
}

int get_matrix(int n, const char* filename, double* A)
{
    FILE *in = fopen(filename, "r");
    if(in == NULL) {
        std::cerr << "BAD_FILE\n";
        return -1; 
    }
    int flag = 1;
    int counter = 0;
    double curr;

    while((counter < n * n) && (flag = fscanf(in, "%lf", &curr) == 1))
    {
        A[counter] = curr;
        counter++;
    }

    if (counter != n * n) {
        std::cerr << counter << ": counter\n";
        std::cerr << "BAD NUMBER OF ELEM\n";
        return -2; 
    }

    if (flag == 0) {
        std::cerr << counter << ": counter\n";
        std::cerr << flag << ": flag\n";
        std::cerr << "DATA FORMAT\n";
        return -3; 
    }
    

    return 0;
}


#endif
