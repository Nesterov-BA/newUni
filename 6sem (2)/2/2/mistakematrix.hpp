#ifndef MYMISTAKE
#define MYMISTAKE



void calc_mistake(int n, double* A, double* ans, double* temp)
{
    double EPS = 1e-15;
    int i;
    int j;
    int k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
            {
                temp[i*n + j] += A[i*n + k] * ans[k*n + j];
            }
        }
    }
    for (i = 0; i < n; i++)
    {
        ans[i*n + i]--;
    }
    double max = 0.0;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            max += dabs(ans[i*n+j]);
        }
    }
    printf("НЕВЯЗКА = %f\n", max);
}


#endif