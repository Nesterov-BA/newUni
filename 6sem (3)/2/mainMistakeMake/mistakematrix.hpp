#ifndef MYMISTAKE
#define MYMISTAKE



void calc_mistake(double* res, int size)
{
    double result;
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
    printf("start = %d, end = %d\n", start, end);
    for (i = start; i < end; i++)
    {
        for (j = 0; j < size; j++)
        {
            for (k = 0; k < size; k++)
            {
                res[i*size + j] += A[i*size + k] * B[k*size + j];
            }
        }
    }
}


#endif