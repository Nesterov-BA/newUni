#include <stdio.h>
#include <iostream>
#include <time.h>
#include <pthread.h>
#include "getmatrix.hpp"
#include "myalgo.hpp"
#include "printmatrix.hpp"
#include "mistakematrix.hpp"
#include "get_time.hpp"

typedef struct _ARGS
{
  double *matrix; 
  double* Temp;
  double* Temp2; 
  int n;
  double* Norm;                
  int thread_num;       /* номер задачи */
  int total_threads;    /* всего задач */
} ARGS;

static long int threads_total_time = 0;
/* Объект типа mutex для синхронизации доступа к
   threads_total_time */
static pthread_mutex_t threads_total_time_mutex 
  = PTHREAD_MUTEX_INITIALIZER;

static void * mistmatr (void *pa);
void matrix_init(int n,  double * initmat); 
int get_matrix(int n, const char* filename, double* A);
void print_matrix(int n, const double* A); 
int my_algo(int n, double* A, double* ans);
double func(int k, int n, int i, int j); 
void calc_mistake(double* A, double* ans, double* res, int n, double* norm, int thread_num, int total_threads);

void calc_mistake(double* A, double* ans, double* res, int n, double* norm, int thread_num, int total_threads)
{
    //double EPS = 1e-15;
    int first_row, last_row;
    first_row = n * thread_num;
    first_row /= total_threads;
    last_row = n * (thread_num + 1);
    last_row = last_row / total_threads - 1;
    
            printf("thread_num: %d, first_row: %d, last_row: %d\n", thread_num, first_row, last_row);
        
        for (int j = first_row; j <= last_row; ++j)
            {
                for (int i = 0; i < n; ++i)
                {
                    for (int k = 0; k < n; ++k)
                    {
                        res[j*n+i] += A[j*n+k]*ans[k*n+i];
                    }
                }
            }
    


}

static void * mistmatr (void *pa)
{
  ARGS *pargs = (ARGS*)pa;
  long int t;
  t = get_time (); 

     calc_mistake(pargs->matrix, pargs->Temp, pargs->Temp2, pargs->n, pargs->Norm, pargs->thread_num, pargs->total_threads);
  t = get_time () - t;  /* время конца работы */

  /* Суммируем времена работы */
  /* "захватить" mutex для работы с threads_total_time */
  pthread_mutex_lock (&threads_total_time_mutex);
  threads_total_time += t;
  /* "освободить" mutex */
  pthread_mutex_unlock (&threads_total_time_mutex);
  return 0;
}

int main(int argc, char **argv)
{
    std::cerr << argc << "\n";
    if ((argc != 3) && (argc != 4)) {return -1;} 
    int nthreads = atoi(argv[1]);
    int n = atoi(argv[2]); 

    std::cerr << n << "\n";

    int res = 0;
    double norm = 0;
    double *A = new double[n * n];
    double *temp = new double[n * n];
    pthread_t * threads;
     ARGS * args;
    if (argc == 3) {
        matrix_init(n, A);
        matrix_init(n, temp);
    } else {
        res = get_matrix(n, argv[3], A);
        res = get_matrix(n, argv[3], temp);
        if (res != 0)  {return -1;}
    }

    //print_matrix(n, A);
    
    double* ans = new double[n * n];

    //clock_t tStart = clock();
    res = my_algo(n, A, ans);
    //clock_t tEnd = clock();

    if (res != 0)  {return -1;}
    printf("FIn\n");
    

    //print_matrix(n, ans);

   
    double * temp2 = new double[n * n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            temp2[i*n + j] = 0;
        }
    }

       threads = (pthread_t*) malloc (nthreads * sizeof (pthread_t));
    args = (ARGS*) malloc (nthreads * sizeof (ARGS));
    for (int i = 0; i < nthreads; i++)
    {
      args[i].matrix = temp;
      args[i].Temp = ans;
      args[i].Temp2 = temp2;
      args[i].n = n;
      args[i].Norm = &norm;
      args[i].thread_num = i;
      args[i].total_threads = nthreads;
    }

    printf("Not segfault here\n");
    clock_t tStart = clock();

            for (int thr = 0; thr < nthreads; ++thr)
            {
                printf("%d\n", thr);
                if (pthread_create (threads + thr, 0,
                          mistmatr,
                          args + thr))
                {
                    fprintf (stderr, "cannot create thread #%d!\n",
                   thr);
                    return -1;
                }
            }
  /* Ожидаем окончания задач */
            for (int thr = 0; thr < nthreads; thr++)
            {
            if (pthread_join (threads[thr], 0))
                fprintf (stderr, "cannot wait thread #%d!\n", thr);
            }
    printf("thread time: %ld\n", threads_total_time);

    clock_t tEnd = clock();
    printf("Time taken: %.2fs\n", (double)(tEnd - tStart)/CLOCKS_PER_SEC);
    //print_matrix(n,temp2);
    for (int i = 0; i < n; i++)
    {
        temp2[i*n + i]--;

    }
    //print_matrix(n,temp);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            norm += dabs(temp2[i*n+j]);
        }
    }
    //calc_mistake(n, temp, ans, temp2);
    printf("Nevyazka = %f\n", norm);



    delete[] A;
    delete[] ans;
    delete[] temp;
    delete[] temp2;
 
    return 0;
}