#include <stdio.h>
#include <iostream>
#include <time.h>
extern "C"{
#include <pthread.h>
}
#include "getmatrix.hpp"
#include "myalgo.hpp"
#include "printmatrix.hpp"
#include "mistakematrix.hpp"

void matrix_init(int n,  double * initmat); 
int get_matrix(int n, const char* filename, double* A);
void print_matrix(int n, const double* A); 
int my_algo(int n, double* A, double* ans);
double func(int k, int n, int i, int j); 
void calc_mistake(int n, double* A, double* ans, double* temp); 

int main(int argc, char **argv)
{
    std::cerr << argc << "\n";
    if ((argc != 3) && (argc != 4)) {return -1;} 
    int nthreads = atoi(argv[1]); 
    int n = atoi(argv[2]); 

    std::cerr << n << "\n";

    int res = 0;
    double *A = new double[n * n];
    double *temp = new double[n * n];
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
    pthread_t * threads;
  /* массив аргументов для созданных задач */
    ARGS * args;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            if (i==j) {
                ans[j+n*i] = 1;
            } else {
                ans[j+n*i] = 0;
            }

        }
    }

    threads = (pthread_t*) malloc (nthreads * sizeof (pthread_t));
    args = (ARGS*) malloc (nthreads * sizeof (ARGS));
    printf("nthreads = %d\n", nthreads);
    for (int i = 0; i < nthreads; i++)
    {
      args[i].matrix = A;
      args[i].result = ans;
      args[i].n = n;
      args[i].thread_num = i;
      args[i].total_threads = nthreads;
    }

    clock_t tStart = clock();

            for (int thr = 0; thr < nthreads; thr++)
            {
                if (pthread_create (threads + thr, 0,
                          rowsum_thr,
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

    if (res != 0)  {return -1;}

    printf("Time taken: %.2fs\n", (double)(tEnd - tStart)/CLOCKS_PER_SEC);

    //print_matrix(n, ans);

   
    double * temp2 = new double[n * n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            temp2[i*n + j] = 0;
        }
    }

    //calc_mistake(n, temp, ans, temp2);
    



    delete[] A;
    delete[] ans;
    delete[] temp;
    delete[] temp2;
 
    return 0;
}