#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "help.h"
#include "task.h"

//#define INPUT_FILE "data.dat"

typedef struct
{
	int n;
	double *a;
	double *x;
	int *index;
	int my_rank;
	int total_threads;
	double Norm;
	int *Err;
	int S;
} ARGS;

long int thread_time = 0;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
void *Inversion(void *p_arg);
void *Inversion(void *p_arg)
{
	ARGS *arg = (ARGS*)p_arg;
	long int t1;

	t1 = get_time();

	if(InvMatrix(arg->n, arg->a, arg->x, arg->index, arg->my_rank, arg->total_threads, arg->Norm, arg -> Err, arg -> S))
		return NULL;
	t1 = get_time() - t1;

	pthread_mutex_lock(&mutex);
	thread_time += t1;
	pthread_mutex_unlock(&mutex);

	return NULL;
}

// 0       1         2           3    4
// ./a.out n_threads matrix_size mode filename

int main(int argc, char **argv)
{
	int i;
	int n;
    int s;
	double *a;
	double *x;
	int *index;
	int mode;
	int *err = malloc(sizeof(int));
    int rezult; 
	double norm;
	FILE *input;
	long int t_full;
	int total_threads;
	pthread_t *threads;
	ARGS *args;

	if (argc > 1)
		total_threads = atoi(argv[1]);
	else
	{
		printf("Incorrect usage!\n");
		return -1;
	}
if (total_threads < 1)
    {
        printf("[X] Wrong input! total_threads should be > 0\n");
        return -7; 
    }
	
    if((argc != 4) && (argc != 5))
    {
        printf("[X] Wrong input! Usage: ./a.out n_thread matrix_size input_mode [filename]\n");
        return -1;
    }
    
    if(! (n = atoi(argv[2])))
    {
        printf("[X] Wrong input! Matrix size should be a number!\n");
        return -2;
    }
    if (n < 1)
    {
        printf("[X] Wrong input! Matrix size should be > 0\n");
        return -3; 
    }
    
    if (! (mode = atoi(argv[3])))
       {
        printf("[X] Wrong input! Input mode should be a number!\n");
        return -4;
    }    
    if ((mode != 1) && (mode != 2))
    {
        printf("[X] Wrong input! Input mode should be equal to 1 or 2!\n");
        return -5;
    }
    
    if ((mode == 1) && (argc != 5))
    {
        printf("[X] Wrong input! If you choose input_mode 1, you should pass 5 arguments!\n");
        return -6;
    }
    if ((mode == 2) && (argc != 4))
    {
        printf("[X] Wrong input! If you choose input_mode 2, you should pass 4 arguments!\n");
        return -6;
    }
	
	//mode = 2; /* Change this for inputing from file */

	switch (mode)
	{
	case 1:
		if ((input = fopen(argv[4], "r")) == NULL)
		{
			printf("Can't open file \"%s\"!\n", argv[4]);
			return -2;
		}
		//else 
        //    fscanf(input, "%d", &n);
		break;
	case 2:
		input = NULL;
/*
		printf("Matrix size: ");
		scanf("%d", &n);
*/
		if (argc == 4){
			//n = atoi(argv[2]);
            printf("Please, enter S: ");
            scanf ("%d", &s);

            if((s < 1) || (s > 5))
            {
                printf("Incorrect mode!\n");
                return -10;
            }

            //s = atoi(argv[3]);
        }
		else
		{
			printf("Incorrect usage!\n");
			return -3;
		}

		if (n <= 0)
		{
			printf("Incorrect n!\n");
			return -4;
		}
		break;
	default:
		printf("Unknown mode.\n");
		return -5;
	}

	a = (double*)malloc(n * n * sizeof(double));
	x = (double*)malloc(n * n * sizeof(double));
	index = (int*)malloc(n * sizeof(int));
	threads = (pthread_t*)malloc(total_threads * sizeof(pthread_t));
	args = (ARGS*)malloc(total_threads * sizeof(ARGS));
	// printf("Problem not here1\n");
	*err = 0;
	// printf("Problem not here2\n");
	if (!(a && x && index && threads && args))
	{
		printf("Not enough memory!\n");

		if (a) free(a);
		if (x) free(x);
		if (index) free(index);
		if (threads) free(threads);
		if (args) free(args);

		return -4;
	}

	rezult = InputMatrix(n, s, a, mode, input);
    
    if (mode == 1)
		fclose(input);
    
    if (rezult == -1)
	{
		printf("Error in reading from file!\n");

		free(index);
		free(a);
		free(x);
        free(threads);
        free(args);
        
		return -6;
	}
    
	printf("Matrix A:\n");
	OutputMatrix(n, a);

	printf("\nCalculating...\n");
	norm = matrNorm(a, n);
	for (i = 0; i < total_threads; i++)
	{
		args[i].n = n;
		args[i].a = a;
		args[i].x = x;
		args[i].index = index;
		args[i].my_rank = i;
		args[i].total_threads = total_threads;
		args[i].Norm = norm;
		args[i].Err = err;
		args[i].S = mode;
	}

	t_full = get_full_time();

	for (i = 0; i < total_threads; i++)
		if (pthread_create(threads + i, 0, Inversion, args + i))
		{
			printf("Cannot create thread %d!\n", i);

			if (a) free(a);
			if (x) free(x);
			if (index) free(index);
			if (threads) free(threads);
			if (args) free(args);

			return -7;
		}

	for (i = 0; i < total_threads; i++)
		if (pthread_join(threads[i], 0))
		{
			printf("Cannot wait thread %d!\n", i);

			if (a) free(a);
			if (x) free(x);
			if (index) free(index);
			if (threads) free(threads);
			if (args) free(args);

			return -8;
		}

	t_full = get_full_time() - t_full;

	if (t_full == 0)
		t_full = 1;
	

	printf("err = %d\n", *err);
	if(!*err)
	{
		printf("\nMatrix A^{-1}:\n");
		OutputMatrix(n, x);
		printf("\n\nInversion time = %ld\nTotal threads time = %ld"\
		" (%ld%%)\nPer thread = %ld\n",
		t_full, thread_time, thread_time * 100/t_full,
		thread_time/total_threads);

		if (mode == 1)
		{
			input = fopen(argv[4], "r");
			//fscanf(input, "%d", &n);
		}
		InputMatrix(n, s, a, mode, input);
        if (mode == 1)
        {
            fclose(input);
        }

		printf("\n||A*A^{-1} - E|| = %e\n", TestMatrix(n, a, x));

	} else {
		printf("Матрица вырожденная.\n");
	}


	free(a);
	free(x);
	free(index);
	free(threads);
	free(args);
	free(err);
	return 0;
}
