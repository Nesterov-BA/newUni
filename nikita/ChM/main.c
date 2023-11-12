#include "fun.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void)
{
    int N;
    int i;
    int m1, m2;
    int m1_max;
    int m2_max;
    double mera;
    double mera_max;
    double h;
    int m;
    int m_max;
    double *y1;
    double *y2;
    double s;
    double s_max;
    double *A;
    double *Lambda;

    for(i = 2; i<= 10; i++)
    {
        N = pow(2, i);
        y1 = (double *)malloc((N+1) * sizeof(double));
        if(!y1)
        {
            fprintf(stderr, "memory allocation error\n");
            return -1;
        }

         y2 = (double *)malloc((N+1) * sizeof(double));
        if(!y2)
        {
            fprintf(stderr, "memory allocation error\n");
            free (y1);
            return -1;
        }

         A = (double *)malloc((N+1) * sizeof(double));
        if(!A)
        {
            fprintf(stderr, "memory allocation error\n");
            free(y1);
            free(y2);
            return -1;
        }

         Lambda = (double *)malloc((N+1) * sizeof(double));
        if(!Lambda)
        {
            fprintf(stderr, "memory allocation error\n");
            free(y1);
            free(y2);
            free(A);
            return -1;
        }
    

    h = 1./((double)N - 0.5);

    s_max = 0.0;
    m1_max = 0.0;
    m2_max = 0.0;

    for (m1 = 1; m1 < N; m1++)
    {
        for(m2 = m1 + 1; m2 < N; m2++)
        {
            Eigen_Vector(y1, m1, N);
            Eigen_Vector(y2, m2, N);
            s = Scalar_Prod(y1, y2, N, h);
            if (s > s_max)
            {
                m1_max = m1;
                m2_max = m2;
                s_max = s;
            }
        }
    }

    printf("Step 2 : %d\n", i);
    printf("Measure Ortog: %e Number: %d, and %d \n", s_max, m1_max, m2_max);

    mera_max = 0.0;
    m_max = 0;

    for (m = 1; m < N ; m++)
    {
        Lambda[m] = Eigen_Value(m, N, h);
        Eigen_Vector(y1, m, N);
        mera = Measure(A, y1, Lambda[m], N, h);
        if (mera > mera_max)
        {
            mera_max = mera;
            m_max = m;
        }
    }    

    printf ("Measure2: %e Number : %d\n\n", mera_max, m_max);

    free(y1);
    free(y2);
    free(A);
    free(Lambda);
    return 0;
}
}