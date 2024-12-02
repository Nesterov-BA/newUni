#include <iostream>
#include <fstream>
#include <ctime>
#include <math.h>

const int n = 5; //Размерность задачи
const int m = 3; //Количество параметров пристрелки
double alpha;
double p[m]; //Вектор параметров пристрелки
double dr[m][m]; // матрица частных производных
double k[6][n];

/*
Вектор p - параметры пристрелки:
____________________
p[0] - lamda_1
p[1] - lambda_2
p[2] - lamda_4
____________________
x[0] - p1   res[0] = p1'
x[1] - p2   res[1] = p2'
x[2] - y    res[2] = y'
x[3] - x    res[3] = x'
*/
void f(double* x, double* res, double t, double h);
void load(double*x);
void residual(double *x, double *r, double eps);
double metric(double* x, double* y, int k);
double metric_abs(double *x, double *y, int k);
double norm(double *x, int k);
void next_k(double* x, double h, double t);
void next_x(double* x);
double error();
void Runge_Kutta(double* x, double eps);
void Runge_Kutta_write(double* x, double eps);
double shooting(double *x, double eps);
void gauss(double *b);



//Диффур
void f(double* x, double* res, double t, double h){
    res[0] = (p[0])*h;
    res[1] = (-x[0])*h;
    res[2] = ((1+alpha*pow(t, 4))*x[1])*h;
    res[3] = x[2]*h;
    res[4] = x[3]*h;
}

//начальные данные для диффура
void load(double*x){
    x[0] = p[1];
    x[1] = p[2];
    x[2] = 1;
    x[3] = 0;
    x[4] = 0;
}

//вычисление невязки
void residual(double *x, double *r, double eps){
    double x_tmp[n];
    for (int i=0; i<n; i++){
        x_tmp[i] = x[i]; 
    }
    Runge_Kutta(x_tmp, eps);
    r[2] = x_tmp[4]-1;
    r[0] = x_tmp[0];
    r[1] = x_tmp[2];
}

double metric(double* x, double* y, int k){
    double res =0;
    for (int i=0; i<k; i++){
        res+= (x[i]-y[i])*(x[i]-y[i]);
    }
    return sqrt(res);
}

double metric_abs(double *x, double *y, int k){
    double sum =0;
    for (int i=0; i<k; i++){
        sum += fabs(y[i]-x[i]);
    }
    return sum;
}

double norm(double *x, int k){
    double sum =0;
    for (int i=0; i<k; i++){
        sum += fabs(x[i]);
    }
    return sum;
}

void next_k(double* x, double h, double t){
    //-------------------------
    //Подсчет K1 (k[0])
    f(x, k[0], t, h);
    //-------------------------
    //Подсчет K2 (k[1])
    double tmp[n];
    for (int i =0; i<n; i++){
        tmp[i] = x[i]+k[0][i]/2.;
    }
    f(tmp, k[1], t+h/2., h);
    //-------------------------
    //Подсчет K3 (k[2])
    for (int i =0; i<n; i++){
        tmp[i] = x[i]+k[0][i]/4.+k[1][i]/4.;
    }
    f(tmp, k[2], t+h/2., h);
    //-------------------------
    //Подсчет K4 (k[3])
    for (int i =0; i<n; i++){
        tmp[i] = x[i]-k[1][i]+2.*k[2][i];
    }
    f(tmp, k[3], t+h, h);
    //-------------------------
    //Подсчет K5 (k[4])
    for (int i =0; i<n; i++){
        tmp[i] = x[i]+k[0][i]*7./27.+k[1][i]*10./27.+k[3][i]/27.;
    }
    f(tmp, k[4], t+2.*h/3., h);
    //-------------------------
    //Подсчет K6 (k[5])
    for(int i =0; i<n; i++){
        tmp[i] = x[i]+k[0][i]*0.0448-k[1][i]*0.2+k[2][i]*0.8736+k[3][i]*0.0864-k[4][i]*0.6048;
    }
    f(tmp, k[5], t+h/5., h);
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

void Runge_Kutta(double* x, double eps){
    load(x);
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
            next_k(x, h, t);
            err = error();
            fac = fmax(0.1, fmin(5, pow(err/eps, 1./5.)));
            h_new = 0.95*h/fac;
        }
        if (t+h>T){
            h = T-t;
            next_k(x, h, t);
            err = error();
        }
        t+=h;     
        next_x(x);;
    }
}

void Runge_Kutta_write(double* x, double eps){
    load(x);
    double t = 0;
    double T = 1;
    double I=0;
    double x_tmp;
    double h_new, h;
    h_new = 0.01;
    std::ofstream fout("output.csv");
    fout << "t,p1,p2,y,x,z" << std::endl;
    while (t<T){
        double err = 1;
        double fac;
        while (err > eps){
            h = h_new;
            next_k(x, h, t);
            err = error();
            fac = fmax(0.1, fmin(5, pow(err/eps, 1./5.)));
            h_new = 0.95*h/fac;
        }
        if (t+h>T){
            h = T-t;
            next_k(x, h, t);
            err = error();
        }
        t+=h;     
        next_x(x);
        fout << t << "," << x[0] << "," << x[1] << "," << x[2] << "," << x[3] << "," << x[4] << std::endl;
    }
    fout.close();
}

double shooting(double *x, double eps){
    std::cout << "Shooting in progress" << std::endl;
    double x_tmp[n];
    double r_old[m];
    double b[m];
    double r[m];
    double p_old[m];
    double d;
    int counter = 0;
    while (true){
        residual(x, r_old, eps);
        double resi = norm(r_old, m);

        //if (norm(r_old, m) < eps){
        //    return 1;
        //}
        for(int j=0; j<m;j++){ 
            p_old[j]=p[j];
            b[j]=p[j];
        }
        for(int i=0;i<m;i++){
            for(int j=0; j<m;j++){ 
                p[j]=p_old[j];
            }
            
            d=1e-7;
            p[i]+=d;
            residual(x, r, eps);

            for(int j=0; j<m;j++){
                dr[j][i] = (r[j]-r_old[j])/d;
            }
        }
        
        //Получили матрицу частных производных dr, пршлое значение параметров p_old
        gauss(b);
        counter++;
        
        double gamma = 1;
        double tmp;
        while (true){
            for(int i=0;i<m;i++){
                p[i] = p_old[i] - b[i]*gamma;
            }
            
            residual(x, r, eps);
            tmp = norm(r, m);
            if (norm(r, m)<eps){
                return counter;
            }
            else if (norm(r, m) < norm(r_old, m)){
                break;
            }
            else if (gamma>1e-20){
                gamma = gamma/2;
            }
            else {
                throw std::runtime_error("Gamma less then eps");
            }

        }
        if ((counter%10 == 0) && (counter !=0)) {
            std::cout << "Counter: " << counter << ", residual: " << norm(r, m) << std::endl;
        }
        if (counter > 10000){
            std::cout << std::endl;
            throw std::runtime_error("Counter too big");
        }
                
    }
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
            dr[i][j] = dr[mxn][i];
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


/*
Вектор p - параметры пристрелки:
____________________
p[0] - lamda_1 - p1'
p[1] - lambda_2 - p1(0)
p[2] - lamda_4 - p2(0)
____________________
x[0] - p1   res[0] = p1'
x[1] - p2   res[1] = p2'
x[2] - y    res[2] = y'
x[3] - x    res[3] = x'
*/
int main(){
    double eps = 1e-7;
    int counter;
    alpha = 0.0001;
    p[0] = -30;
    p[1] = 30;
    p[2] = 9;
    double x[n]; 
    counter = shooting(x, eps);
    std::cout << "Parametrs p:" << std::endl;
    std::cout << "------------------------" << std::endl;
    for(int jm=0;jm<m;jm++){
        std::cout << p[jm] << " ";
    }
    std::cout << std::endl;
    std::cout << "------------------------" << std::endl;
    Runge_Kutta_write(x, eps);




    return 1;
}