#include "help.hpp"

double f0(double x){
    if(x>0)
        return 1;
    return 1;

}

double f1(double x){
        return x;
}

double f2(double x){
        return x*x;
}

double f3(double x){
    return x*x*x;
}

double f4(double x){
    return x*x*x*x;
}

double f5(double x){
    return exp(x);
}

double f6(double x){
    return 1/(25*x*x+1);
}

double d2f0(double x){
    if(x>0)
        return 0;
    return 0;
}

double d2f1(double x ){
    if(x > 0)
        return 1;
    return 1;
}

double d2f2(double x ){
    return 2*x;
}

double d2f3(double x ){
    return 3*x*x;
}

double d2f4(double x ){
    return 4*x*x*x;
}

double d2f5(double x ){
    return exp(x);
}

double d2f6(double x ){
    return (5000*pow(x,2))/pow(25*pow(x,2)+1, 3) - 50/pow(25*pow(x,2)+1, 2);
}


double max_matr(double *M, int size){
    double res=M[0];
    for(int i = 0; i < size; i++)
    {
        if(M[i] > res)
            res = M[i];
    }
    return res;
}

double min_matr(double *M, int size){
    double res=M[0];
    for(int i = 0; i < size; i++)
    {
        if(M[i] < res)
            res = M[i];
    }
    return res;
}

double max4(double a, double b, double c, double d){
    double res=a;
    if(res<b)
        res=b;
    if(res<c)
        res=c;
    if(res<d)
        res=d;
    return res;
}

double min4(double a, double b, double c, double d){
    double res=a;
    if(res>b)
        res=b;
    if(res>c)
        res=c;
    if(res>d)
        res=d;
    return res;
}

double max(double a, double b){
    if(a>b)
        return a;
    return b;
}

double min(double a, double b){
    if(a>b)
        return b;
    return a;
}

double diffFun2(double val1, double val2, double x1, double x2)
{
    double diff = (val2 - val1) / (x2 - x1);
    return diff;
}
