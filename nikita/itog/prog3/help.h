#ifndef __HELP_H_INCLUDED__
#define __HELP_H_INCLUDED__



int InputMatrix(int n,int s, double *a, int mode, FILE *input);

void OutputMatrix(int n, double *a);

double TestMatrix(int n, double *a, double *x);

double matrNorm(double* matrix, int size);

long int get_time(void);

long int get_full_time(void);

#endif /* __HELP_H_INCLUDED__ */
