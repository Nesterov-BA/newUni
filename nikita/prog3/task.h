#ifndef __TASK_H_INCLUDED__
#define __TASK_H_INCLUDED__

int InvMatrix(int n, double *a, double *x, int *index, int my_rank, int total_threads, double norm, int* err, int s);

#endif /* __TASK_H_INCLUDED__ */
