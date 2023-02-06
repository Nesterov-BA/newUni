#ifndef MYALGO
#define MYALGO
#include "get_time.hpp"

typedef struct _ARGS
{
  double *matrix;       
  double *result; 
  int n;                /* размер матрицы и векторов */
  int thread_num;       /* номер задачи */
  int total_threads;    /* всего задач */
} ARGS;

static long int threads_total_time = 0;
/* Объект типа mutex для синхронизации доступа к
   threads_total_time */
static pthread_mutex_t threads_total_time_mutex 
  = PTHREAD_MUTEX_INITIALIZER;

void synchronize(int total_threads);
void swaprows(double* matr, int size, int i, int j);
void rowtimesx(double* matr, int size, int i, double x);
void rowplusrowtimesx(double* matr, int size, int i, int j, double x);
double dabs(double x);
void rowsum(double* A, double* ans, int n, int thread_num, int total_threads);
static void * rowsum_thr (void *pa);

void synchronize(int total_threads)
{
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
	static int threads_in = 0;
	static int threads_out = 0;

	pthread_mutex_lock(&mutex);

	threads_in++;
	if (threads_in >= total_threads)
	{
		threads_out = 0;
		pthread_cond_broadcast(&condvar_in);
	} else
		while (threads_in < total_threads)
			pthread_cond_wait(&condvar_in,&mutex);

	threads_out++;
	if (threads_out >= total_threads)
	{
		threads_in = 0;
		pthread_cond_broadcast(&condvar_out);
	} else
		while (threads_out < total_threads)
			pthread_cond_wait(&condvar_out,&mutex);

	pthread_mutex_unlock(&mutex);
}

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

double dabs (double x){
	if (x<0) {
		return -x;
	} else {
		return x;
	}
}

void rowsum(double* A, double* ans, int n, int thread_num, int total_threads)
{
	int first_row, last_row;
	//printf("hey\n");
	for (int i = 0; i < n; ++i)
	{
		if (thread_num == 0)
		{
			rowtimesx(ans, n, i, 1/A[i*n+i]);
			rowtimesx(A, n, i, 1/A[i*n+i]);
		}
		synchronize(total_threads);
		first_row = n * thread_num;
		first_row /= total_threads;
		last_row = n * (thread_num + 1);
		last_row = last_row / total_threads - 1;
		/*if (i == 0)
		{
			printf("thread_num: %d, first_row: %d, last_row: %d\n", thread_num, first_row, last_row);
		}*/
		for (int j = first_row; j <= last_row; ++j)
			{
				if (j!=i){
					rowplusrowtimesx(ans, n, j, i, -A[j*n+i]);
					rowplusrowtimesx(A, n, j, i, -A[j*n+i]);
				}
		}
		synchronize(total_threads);
	}
}

static void * rowsum_thr (void *pa)
{
  ARGS *pargs = (ARGS*)pa;
  long int t;
  t = get_time (); 

     rowsum(pargs->matrix, pargs->result, pargs->n, pargs->thread_num, pargs->total_threads);
  t = get_time () - t;  /* время конца работы */

  /* Суммируем времена работы */
  /* "захватить" mutex для работы с threads_total_time */
  pthread_mutex_lock (&threads_total_time_mutex);
  threads_total_time += t;
  /* "освободить" mutex */
  pthread_mutex_unlock (&threads_total_time_mutex);
  return 0;
}




#endif