#include <iostream>
#include <fstream>
#include <ctime>
#include <math.h>

const int n = 5; //Размерность задачи
const int m = 3; //Количество параметров пристрелки
double alpha;
double eps = 1e-7;
double dr[m][m]; // матрица частных производных
double k[6][n];

void f(double* x, double* res, double*p, double t, double h);
void load(double*x, double*p);
void next_k(double* x, double*p, double h, double t);
void next_x(double* x);
double error();
void Runge_Kutta(double* x, double*p);
void Runge_Kutta_write(double* x, double*p);
void residual(double *p, double *r);
double norm(double *x, int k);
void gauss(double *b);

//Диффур
void f(double* x, double* res, double*p, double t, double h){
    res[0] = (p[0])*h;
    res[1] = (-x[0])*h;
    res[2] = ((1+alpha*pow(t*x[3], 2))*x[1])*h;
    res[3] = x[2]*h;
    res[4] = x[3]*h;
}

//начальные данные для диффура
void load(double*x, double *p){
    x[0] = p[1];
    x[1] = p[2];
    x[2] = 1;
    x[3] = 0;
    x[4] = 0;
}

void next_k(double* x,  double*p, double h, double t){
    //-------------------------
    //Подсчет K1 (k[0])
    f(x, k[0], p, t, h);
    //-------------------------
    //Подсчет K2 (k[1])
    double tmp[n];
    for (int i =0; i<n; i++){
        tmp[i] = x[i]+k[0][i]/2.;
    }
    f(tmp, k[1], p, t+h/2., h);
    //-------------------------
    //Подсчет K3 (k[2])
    for (int i =0; i<n; i++){
        tmp[i] = x[i]+k[0][i]/4.+k[1][i]/4.;
    }
    f(tmp, k[2], p, t+h/2., h);
    //-------------------------
    //Подсчет K4 (k[3])
    for (int i =0; i<n; i++){
        tmp[i] = x[i]-k[1][i]+2.*k[2][i];
    }
    f(tmp, k[3], p, t+h, h);
    //-------------------------
    //Подсчет K5 (k[4])
    for (int i =0; i<n; i++){
        tmp[i] = x[i]+k[0][i]*7./27.+k[1][i]*10./27.+k[3][i]/27.;
    }
    f(tmp, k[4], p, t+2.*h/3., h);
    //-------------------------
    //Подсчет K6 (k[5])
    for(int i =0; i<n; i++){
        tmp[i] = x[i]+k[0][i]*0.0448-k[1][i]*0.2+k[2][i]*0.8736+k[3][i]*0.0864-k[4][i]*0.6048;
    }
    f(tmp, k[5], p, t+h/5., h);
    //-------------------------
}

void next_x(double* x){
    for (int i=0; i<n; i++){
        x[i] = x[i]+ k[0][i]/24.+k[3][i]*5./48.+k[4][i]*27./56.+k[5][i]*125./336.;
    }
}

double error(){
    double err = 0;
    double tmp = 0;
    for (int i=0; i<n; i++){
        tmp = -0.125*k[0][i]-k[2][i]*224./336.-k[3][i]*0.0625+k[4][i]*162./336.+k[5][i]*125./336.;
        err+=tmp*tmp;
    }
    err = sqrt(err);
    return err;
}

void Runge_Kutta(double* x, double *p){
    double t = 0;
    double T = 1;
    double x_tmp;
    double I=0;
    double fac;
    double h_new, h;
    h_new = 0.01;
    while (t<T){
        double err = 1;
        while (err > eps){
            h = h_new;
            next_k(x, p, h, t);
            err = error();
            fac = fmax(0.1, fmin(5, pow(err/eps, 1./5.)));
            h_new = 0.95*h/fac;
        }
        if (t+h>T){
            h = T-t;
            next_k(x, p, h, t);
            err = error();
        }
        t+=h;     
        next_x(x);
    }
}

void Runge_Kutta_write(double* x, double *p){
    double t = 0;
    double T = 1;
    double I=0;
    double x_tmp;
    double h_new, h;
    h_new = 0.01;
    std::ofstream fout("output" +std::to_string(alpha)+".csv");
    fout << "t,p1,p2,y,x,z" << std::endl;
    while (t<T){
        fout << t << "," << x[0] << "," << x[1] << "," << x[2] << "," << x[3] << "," << x[4] << std::endl;
        double err = 1;
        double fac;
        while (err > eps){
            h = h_new;
            next_k(x, p, h, t);
            err = error();
            fac = fmax(0.1, fmin(5, pow(err/eps, 1./5.)));
            h_new = 0.95*h/fac;
        }
        if (t+h>T){
            h = T-t;
            next_k(x, p, h, t);
            err = error();
        }
        t+=h;     
        next_x(x);
    }
    fout << t << "," << x[0] << "," << x[1] << "," << x[2] << "," << x[3] << "," << x[4] << std::endl;
    fout.close();
}

void residual(double *p, double *r){
    double x[n];
    load(x, p);
    Runge_Kutta(x, p);
    r[0] = x[0];
    r[1] = x[2];
    r[2] = x[4]-1;
}

double norm(double *x, int k){
    double sum =0;
    for (int i=0; i<k; i++){
        sum += fabs(x[i]);
    }
    return sum;
}

void gauss(double *b){
    double mx;
    int mxn;
    for (int i=0;i<m;i++){
        mx = fabs(dr[i][i]);
        mxn = i;
        //Находим строку с наибольшим значеием на данной итерации
        for (int j=i+1; j<m; j++){
            if (fabs(dr[j][i]) > mx){
                mx = fabs(dr[j][i]);
                mxn = j;
            }
        }

        for (int j=i; j<m; j++){
            mx = dr[i][j];
            dr[i][j] = dr[mxn][j];
            dr[mxn][j] = mx;
        }
        mx = b[i];
        b[i] = b[mxn];
        b[mxn] = mx;

        mx = dr[i][i];
        dr[i][i] = 1;
        b[i]=b[i]/mx;
        for (int j=i+1; j<m;j++){
            dr[i][j]=dr[i][j]/mx;
        }

        for (int j=0; j<m; j++) {
			if(j != i){
				for (int l=i+1; l<m; l++){
                    dr[j][l] -= dr[j][i]*dr[i][l];
                }
                b[j] -= dr[j][i]*b[i];
				dr[j][i] = 0;
			}
        }
    }
}

double shooting(double *p){
    double d=1e-4;
    double r[m], r_new[m], p_new[m], h[m];
    double r_norm, r_new_norm;
    double gamma = 1.;
    int counter = 0;
    while (true){
        residual(p, r);
        r_norm = norm(r, m);
        if (r_norm < eps){
            return 1;
        }
        for (int i=0; i<m;i++){
            h[i] = r[i];
            for (int j=0; j<m;j++){
                
                p_new[j]=p[j];
            }
            p_new[i]+=d;
            residual(p_new, r_new);
            for (int j=0; j<m;j++){
                dr[j][i]=(r_new[j]-r[j])/d;
            }
        }

        gauss(h);
        counter++;

        while (true){
            for(int i=0;i<m;i++){
                p_new[i] = p[i] - h[i]*gamma;
            }
            residual(p_new, r_new);
            r_new_norm = norm(r_new, m);
            if (r_new_norm<eps){
                for(int i=0;i<m;i++){
                    p[i] = p_new[i];
                }
                return counter;
            }
            else if (r_new_norm < r_norm){
                for(int i=0;i<m;i++){
                    p[i] = p_new[i];
                }
                break;
            }
            else if (gamma>1e-20){
                gamma = gamma/2;
            }
            else {
                throw std::runtime_error("Gamma less then eps");
            }
        }
        if (counter > 10000){
            std::cout << std::endl;
            throw std::runtime_error("Counter too big");
        }
    }
}

int main(){
    int counter;
    double p[m];
    double x[n];
    p[0] = -30;
    p[1] = 30;
    p[2] = 9;
    alpha = -0.1;
    for (int i=0; i<=100; i++){
        alpha += 0.1;
        counter = shooting(p);
        if ((i == 0) || (i == 1) || (i == 10) || (i == 100)){
            std::cout << "Parametrs p (alpha=" << alpha << "):" << std::endl;
            std::cout << "------------------------" << std::endl;
            for(int jm=0;jm<m;jm++){
                std::cout << p[jm] << " ";
            }
            std::cout << std::endl;
            std::cout << "------------------------" << std::endl;
            load(x, p);
            Runge_Kutta_write(x, p);
        }
    }

    return 1;
}