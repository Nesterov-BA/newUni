#ifndef MYPRINTER
#define MYPRINTER

void print_matrix(int n, const double* A)
{
 
    int i = 0;
    int j = 0;

    while (i < n){
        j = 0;
        while (j < n) 
        {
            printf(" %f", A[i*n + j]);
            j++;
        }
        printf("\n");
        i++;
    }
}


#endif
