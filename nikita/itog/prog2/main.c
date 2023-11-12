#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "matrix.h"
#include "task.h"
#include "oldtask.h"

//#define INPUT_FILE_NAME "data.txt"

// 0       1           2    3
// ./a.out MATRIX_SIZE MODE FILENAME
// ./a.out MATRIX_SIZE MODE

int main(int argc, const char * argv[])
{
	int i;
	int j;
	int n;
    int s;
	int iter;
	int rezult;
	int inputMode;
	double eps = 1e-15;
	double inv1;
	double inv2;
    
	double* a = NULL;
	double* values = NULL;
	
        //double* x;
        //double* y;

	double* diag;
	double* side;
	double* sideup;

	FILE* fin = NULL;
	clock_t t;

        /*
	printf("Input mode : 1 - from file \"%s\".\n", INPUT_FILE_NAME);
	printf("             2 - from formula.\n");
	printf("-> ");
	scanf("%d", &inputMode);
        */

        if((argc != 3) && (argc != 4))
        {
            printf("[X] Wrong input! Usage: ./a.out matrix_size input_mode [filename]\n");
            return -1;
        }

        if(! (n = atoi(argv[1])))
        {
            printf("[X] Wrong input! Matrix size should be a number!\n");
            return -2;
        }
        if (n < 1)
        {
            printf("[X] Wrong input! Matrix size should be > 0\n");
            return -3;
        }

        if (! (inputMode = atoi(argv[2])))
           {
            printf("[X] Wrong input! Input mode should be a number!\n");
            return -4;
        }
        if ((inputMode != 1) && (inputMode != 2))
        {
            printf("[X] Wrong input! Input mode should be equal to 1 or 2!\n");
            return -5;
        }

        if ((inputMode == 1) && (argc != 4))
        {
            printf("[X] Wrong input! If you choose input_mode 1, you should pass 4 arguments!\n");
            return -6;
        }
        if ((inputMode == 2) && (argc != 3))
        {
            printf("[X] Wrong input! If you choose input_mode 2, you should pass 3 arguments!\n");
            return -6;
        }

	// INPUT MODE
	{switch (inputMode)
	{
	case 1:
                fin = fopen(argv[3], "r");
		if (!fin)
		{
			printf("[X] Open file!\n");
			return -1;
		}
                /*if (fscanf(fin, "%d", &n) != 1)
		{
			printf("[X] Read from file!\n");
			fclose(fin);
			return -2;
                }*/
		break;
	case 2:
                //printf("[>] Enter N --> ");
                //scanf("%d", &n);
                printf("[>] Enter S --> ");
                scanf("%d", &s);

                if((s < 1) || (s > 4))
                {
                    printf("[X] Wrong s!\n");
                    return -10;
                }

                break;
	default:
		printf("[X] Wrong mode!\n");
		return -3;
	}

	if (n < 1)
	{
		printf("[X] Wrong N\n");
		if (inputMode == 1)
			fclose(fin);
		return -4;
	}}

	// MEMORY
  	{a = (double*)malloc(n * n * sizeof(double));
	values = (double*)malloc(n * sizeof(double));

	diag = (double*) calloc (n, sizeof (double));
  	side = (double*) calloc (n, sizeof (double));
  	sideup = (double*) calloc (n, sizeof (double));
  	
  	if (!(a && values && diag && side && sideup))
	{
		printf("Not enough memory!\n");

		if (a)
			free(a);
		if (values)
			free(values);
		if (diag)
			free(diag);
		if (side)
			free(side);
		if (sideup)
			free(sideup);

		if (inputMode == 1)
			fclose(fin);

		return -5;
	}

	rezult = InputMatrix(n, s, a, inputMode, fin);
    

	if (inputMode == 1)
		fclose(fin);

	if (rezult == -1)
	{
		printf("[X] Read from file!\n");

		free(a);
		free(values);
		free(diag);
		free(side);
        free(sideup);
		return -6;
	}
	}

	// PRINT
	{printf("\n[A]:\n");
	PrintMatrix(n, a);
	printf("\n");}

	
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < i; j++)
		{
			if (fabs(a[i*n + j] - a[j*n + i]) > 1e-40)
			{
				printf("[X] Matrix Not Symmetrical!\n");
                free(a);
		free(values);
		free(diag);
		free(side);
        free(sideup);
				return -7;
			}
		}
	}
	// TO TRIDIAG
	{
	
	tridiag(a, n);

	printf("\n[Tridiagonal  A]:\n");
	PrintMatrix(n, a);
	printf("\n");}

	// RESIDUALS
	{inv1 = 0.0;
	inv2 = 0.0;
	for (i = 0; i < n; ++i)
	{
		inv1 -= a[i * n + i];
		for (j = 0; j < n; ++j)
			inv2 -= a[i * n + j] * a[j * n + i];
	}}

	// MEMORY OPTIMIZATION
	{optimize(a, diag, side, n);
	free(a);
	}

	// PRINT
	{printf("\n[Optimized A]:\n");
	PrintVector(n, diag);
	PrintVector(n-1, side);
	printf("\nOR\n");
	PrintTriDMatrix2(n, diag, side);
	printf("\n");
	}

	// TASK
	{t = clock();
	
	FindValues(n, diag, side, sideup, values, eps, &iter);
	
	
	t = clock() - t;

	for (i = 0; i < n; ++i)
	{
		inv1 += values[i];
		inv2 += values[i] * values[i];
	}}
	

	// PRINT ANSWER
	{printf("\n[>] VALUES: ");
	PrintVector(n, values);
	printf("\n");

	printf("[>] TIME:  %5.5fs.\n", (double)t / CLOCKS_PER_SEC);
	printf("[>] ITERS: %d\n\n", iter);

	// Невязки
	printf("E(xi) - E(ai):          %g\n", inv1);
	printf("E(xi^2) - E(aij * aji): %g\n", inv2);}

	// FREE
	{free(values);
	free(diag);
	free(side);
	free(sideup);
	}

	return 0;
}
