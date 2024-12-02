#ifndef MYMISTAKE
#define MYMISTAKE



#include <cstdio>
void calc_mistake(double* res, int size)
{
    double result = 0;
    for (int i = 0; i < size; i++)
        res[i*size + i]--;
    for(int i = 0; i < size; i++)
    {
        for(int j = 0; j < size; j++)
        {
            result += dabs(res[i*size + j]);
        }
    }
    printf("Невязка = %f\n", result);
}

void mpiMatrMult(double* A, double* B, double* res, int size, int rank, int rowSize)
{
    int i, j, k;
    int start, end;
    
    if(rank < size - 1)
    {
        start = rank*rowSize;
        end = start + rowSize;
    
    } else {
        start = rank*rowSize;
        end = size;
    
    }
    printf("start = %d, end = %d, rank = %d\n", start, end, rank);
    printf("A [] = %f\n", A[start*size]);
    for (i = start; i < end; i++)
    {
        for (j = 0; j < size; j++)
        {
            for (k = 0; k < size; k++)
            {
                res[i*size + j] += A[i*size + k] * B[k*size + j];
            }

            printf("res = %f, i, j = %d, %d\n", res[i*size+j], i, j);
        }
    }
}


#endif