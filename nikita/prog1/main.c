#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "matrix.h"
#include "task.h"

//#define INPUT_FILE_NAME "data.dat"

// 0       1           2    3
// ./a.out MATRIX_SIZE MODE FILENAME
// ./a.out MATRIX_SIZE MODE


int main(int argc, const char* argv[])
{
	int n;
        int s = 1;
	int rezult;
	int inputMode = 0;
	int* index = NULL;
	double* a = NULL;
	double* x = NULL;
	FILE* fin = NULL;
	clock_t t;
	double norm;

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
    
    switch (inputMode)
	{
	case 1:
		fin = fopen(argv[3], "r");

		if (!fin)
		{
			printf("Cann't open file!\n");

			return -1;
		}

		/*if (fscanf(fin, "%d", &n) != 1)
		{
			printf("Error in reading from file!\n");

			fclose(fin);

			return -2;
		}*/   

		break;
	case 2:
		//printf("Please, enter N: ");
		//scanf("%d", &n);
                printf("Please, enter S: ");
                scanf ("%d", &s);
                       

                if((s < 1) || (s > 5))
                {
                    printf("Incorrect mode!\n");
                    return -10;
                }

		 break;
	default:
		printf("[X] Wrong mode!\n");
		return -3;
	}


	if (n < 1)
	{
		printf("Incorrect N!\n");

		if (inputMode == 1)
			fclose(fin);

		return -4;
	}

	index = (int*)malloc(n * sizeof(double));
	a = (double*)malloc(n * n * sizeof(double));
	x = (double*)malloc(n * n * sizeof(double));

	if (!(index && a && x))
	{
		printf("Not enough memory!\n");

		if (index)
			free(index);
		if (a)
			free(a);
		if (x)
			free(x);

		if (inputMode == 1)
			fclose(fin);

		return -5;
	}

	rezult = InputMatrix(n, s, a, inputMode, fin);

	if (inputMode == 1)
		fclose(fin);

	if (rezult == -1)
	{
		printf("Error in reading from file!\n");

		free(index);
		free(a);
		free(x);

		return -6;
	}
	norm = matrNorm(a, n);
	printf("\nMatrix A:\n");
	PrintMatrix(n, a);
	printf("\n");
	PrintNorm(n, a);
	printf(", funNorm = %f\n", norm);
	

	printf("Calculating...\n");
        PrintVector(n, a);
        //PrintVector(n, x);

	t = clock();
	rezult = InvertMatrix(n, a, x, index, norm, inputMode);
	t = clock() - t;

        PrintVector(n, a);
        PrintVector(n, x);

	switch (rezult)
	{
	case -1:
		printf("Can't invert - matrix is deteriorated.\n");
		break;
	case 0:
		printf("\nMatrix A^{-1}:\n");
		PrintMatrix(n, x);
		printf("\n");

		printf("Inversion time\t\t= %.2f sec.\n\n", (double)t / CLOCKS_PER_SEC);

		if (inputMode == 1)
		{
			fin = fopen(argv[3], "r");
                        //fscanf(fin, "%d", &n);
		}

               
                PrintVector(n, a);
                PrintVector(n, x);
                printf("\n");

                InputMatrix(n, s, a, inputMode, fin);

                if (inputMode == 1)
			fclose(fin);

                
                PrintVector(n, a);
                PrintVector(n, x);
                printf("\n");

		printf("Solution ||A * A^{-1} - b||\t= %e\n", SolutionError(n, a, x));
             
		break;
	default:
		printf("Unknown error!\n");

		break;
	}

	free(index);
	free(a);
	free(x);

	return 0;
}
