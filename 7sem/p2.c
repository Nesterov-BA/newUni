#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#define MAX_CNT 8

int Solver(double* u,int N, int M, double T, FILE* plot);

double f_w(double x){
    return 1-16*x*x*(1-x)*(1-x);
}

double integral(double t, double eps){
    int k;
    int N;
    int M;
    double h;
    double tau;
    double intgrl_sum;
    double intgrl_sum_s;
    double tmp;
    double *u;

    intgrl_sum_s = 0;
    k = 16;
    for(int cnt = 0; cnt < MAX_CNT; cnt++, k*=2){

        M = k;
        N = 2*t*k;

        u = (double*)malloc((N+2)*(M+1)*sizeof(double));

        h = 1.0/M;
        tau = t/N;
        Solver(u,(N+1),M,tau*(N+1), NULL);
        intgrl_sum = 0;
        for(int m = 0; m < M-1; m++){
            tmp = ((u[(N+1)*(M+1)+m]-u[(N)*(M+1)+m])/tau)*((u[(N+1)*(M+1)+m]-u[(N)*(M+1)+m])/tau)+((u[(N)*(M+1)+m+1]-u[(N)*(M+1)+m])/h)*((u[(N)*(M+1)+m+1]-u[(N)*(M+1)+m])/h)+2*f_w(h*m)*u[(N)*(M+1)+m]*u[(N)*(M+1)+m];
            tmp += ((u[(N+1)*(M+1)+m+1]-u[(N)*(M+1)+m+1])/tau)*((u[(N+1)*(M+1)+m+1]-u[(N)*(M+1)+m+1])/tau)+((u[(N)*(M+1)+m+1+1]-u[(N)*(M+1)+m+1])/h)*((u[(N)*(M+1)+m+1+1]-u[(N)*(M+1)+m+1])/h)+2*f_w(h*(m+1))*u[(N)*(M+1)+m+1]*u[(N)*(M+1)+m+1];
            intgrl_sum += 0.5*tmp;
        }
        tmp = ((u[(N+1)*(M+1)+(M-1)]-u[(N)*(M+1)+(M-1)])/tau)*((u[(N+1)*(M+1)+(M-1)]-u[(N)*(M+1)+(M-1)])/tau)+((u[(N)*(M+1)+(M-1)+1]-u[(N)*(M+1)+(M-1)])/h)*((u[(N)*(M+1)+(M-1)+1]-u[(N)*(M+1)+(M-1)])/h)+2*f_w(h*(M-1))*u[(N)*(M+1)+(M-1)]*u[(N)*(M+1)+(M-1)];
        tmp +=((u[(N+1)*(M+1)+(M)]-u[(N)*(M+1)+(M)])/tau)*((u[(N+1)*(M+1)+(M)]-u[(N)*(M+1)+(M)])/tau)+((u[(N)*(M+1)+M]-u[(N)*(M+1)+(M-1)])/h)*((u[(N)*(M+1)+M]-u[(N)*(M+1)+(M-1)])/h)+2*f_w(h*(M))*u[(N)*(M+1)+(M)]*u[(N)*(M+1)+(M)];
        intgrl_sum += 0.5*tmp;
        intgrl_sum *= h;

        if(fabs(intgrl_sum_s-intgrl_sum) < eps){
            free(u);
            break;
        }
        intgrl_sum_s = intgrl_sum;
        //printf("k = %d, %lf\n",k,intgrl_sum);
        free(u);
    }
    return intgrl_sum;
}

int Solver(double* u,int N, int M, double T, FILE* plot){      //u[n*M+m] = u(tau*n,h*m)

    double tmp;

    double tau = T/N;
    double h = 1.0/M;

    for(int m = 0; m < M+1; m++){
        u[0*(M+1)+m] = 0;
    }

    for(int m = 1; m < M; m++){
        u[1*(M+1)+m] = 0;
    }

    u[1*(M+1)+0] = 0;
    u[1*(M+1)+M] = 0;

    if(plot != NULL){
        for(int n = 0; n < 2; n++){
            for(int m = 0; m < M+1; m++){
                fprintf(plot,"%10.7lf %10.7lf %10.7lf\n",m*h,n*tau,u[n*(M+1)+m]);
            }
        }
    }

    //tau/h < 1;  M*T/N < 1;

    for(int n = 1; n+1 < N+1; n++){
        for(int m = 1; m < M; m++){
            u[(n+1)*(M+1)+m] = 2*u[n*(M+1)+m]-u[(n-1)*(M+1)+m]+tau*tau*((u[n*(M+1)+m+1]-2*u[n*(M+1)+m]+u[n*(M+1)+m-1])/(h*h)-2*f_w(m*h)*u[n*(M+1)+m]);
        }

       tmp = u[(n-1)*(M+1)+0]+2*tau*((u[n*(M+1)+0+1]-u[n*(M+1)+0])/h)-h*(tau*2*f_w(h*0)*u[n*(M+1)+0]-((2*u[n*(M+1)+0]-u[(n-1)*(M+1)+0])/tau))+2*tau*sin(tau*n);
        u[(n+1)*(M+1)+0] = tmp/(1+h/tau);

        tmp = u[(n-1)*(M+1)+(M)]-2*tau*((u[n*(M+1)+(M)]-u[n*(M+1)+(M)-1])/h)-h*(tau*2*f_w(h*(M))*u[n*(M+1)+(M)]-((2*u[n*(M+1)+(M)]-u[(n-1)*(M+1)+(M)])/tau));
        u[(n+1)*(M+1)+(M)] = tmp/(1+h/tau);

        if(plot != NULL){
            for(int m = 0; m < M+1; m++){
                fprintf(plot,"%10.7lf %10.7lf %10.7lf\n",m*h,(n+1)*tau,u[(n+1)*(M+1)+m]);
            }
        }
    }

    return 0;

}


int main()
{
    //M*T/N < 1;
    int K = 25;
    double T = 10.0;
    int N = 2*T*K;
    int M = K;
    double t;
    double eps;
    double *u;

    FILE* plot;
    plot = fopen("plot","w");

    u = (double*)malloc((N+1)*(M+1)*sizeof(double));
    Solver(u, N, M, T, plot);

    /*for(int n = 0; n < N+1; n++){
        for(int m = 0; m < M+1; m++){
            printf("%5.3lf ", u[n*(M+1)+m]);
        }
        printf("\n");
    }*/

    free(u);
    fclose(plot);

    t = 1;
    eps = 1e-4;
    printf("f(%3.1lf) = %5.3lf \n",t,integral(t,eps));

    return 0;
}
