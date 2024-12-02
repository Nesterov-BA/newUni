#include "help.hpp"

double f0(double x, double y){
    if(x>0 && y>0)
        return 1;
    return 1;
}

double f1(double x, double y){
    if(x>0 && y>0)
        return x;
    return x;
}

double f2(double x, double y){
    if(x>0 && y>0)
        return y;
    return y;
}

double f3(double x, double y){
    return x+y;
}

double f4(double x, double y){
    return sqrt(x*x+y*y);
}

double f5(double x, double y){
    return x*x+y*y;
}

double f6(double x, double y){
    return exp(x*x-y*y);
}

double f7(double x, double y){
    return 1/(25*(x*x+y*y)+1);
}

double d2f0(double x, double y, int arg){
    return 0;
}

double d2f1(double x, double y, int arg){
    return 0;
}

double d2f2(double x, double y, int arg){
    return 0;
}

double d2f3(double x, double y, int arg){
    return 0;
}

double d2f4(double x, double y, int arg){ // 0 for der with respect to x
    if(arg == 0)
        return y*y/sqrt(pow(x*x+y*y, 3));
    return x*x/sqrt(pow(x*x+y*y, 3));
}

double d2f5(double x, double y, int arg){
    return 2;
}

double d2f6(double x, double y, int arg){
    if(arg == 0)
        return 2*(2*x*x+1)*exp(x*x-y*y);
    return 2*(2*y*y-1)*exp(x*x-y*y);
}

double d2f7(double x, double y, int arg){
    if(arg == 0)
        return 5000*x*x/pow((25*(x*x+y*y)+1), 3) + 50/pow((25*(x*x+y*y)+1), 2);
    return 5000*y*y/pow((25*(x*x+y*y)+1), 3) + 50/pow((25*(x*x+y*y)+1), 2);
}
double max_matr(double *M, int nx, int ny){
    double res=M[0];
    for(int i=0;i<nx;++i){
        for(int j=0;j<ny;++j){
            if(M[i*ny+j]>res)
                res=M[i*ny+j];
        }
    }
    return res;
}

double min_matr(double *M, int nx, int ny){
    double res=M[0];
    for(int i=0;i<nx;++i){
        for(int j=0;j<ny;++j){
            if(M[i*ny+j]<res)
                res=M[i*ny+j];
        }
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

void transponse(double* matrix, int size)
{
    double buffer;
    for(int row = 0; row < size-1; row++)
    {
        for(int col = row + 1; col < size; col++)
        {
            buffer = matrix[row*size + col];
            matrix[row*size + col] = matrix[col*size + row];
            matrix[col*size + row] = buffer;
        }
    }
}

void sqMatrMult(double* left, double* right, int size, int pos)
{
    double* temp = new double[size*size];
    for(int row = 0; row < size; row++)
        {
            for(int col = 0; col < size; col++)
            {
                temp[row*size+col] = 0;
                for(int k = 0; k < size; k++)
                    temp[row*size + col] += left[row*size + k]*right[k*size + col];
            }
        }
    if(pos)
    {
        for(int i = 0; i < size*size; i++)
            right[i] = temp[i];
    } else {
        for(int i = 0; i < size*size; i++)
            left[i] = temp[i];
    }
    delete [] temp;
}
