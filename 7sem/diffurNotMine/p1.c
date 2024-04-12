#include<stdio.h>
#include<math.h>

void f(double f_val[2], double x, double y[2]);
void f(double f_val[2], double x, double y[2]){
    double a = 10;
    f_val[0] = y[1];
    f_val[1] = a*(1-y[0]*y[0])*y[1]-y[0];
}


int RKstep(double y[2], void (*f)(double[2], double, double[2]), double x0, double y0[2], double h){
    double k1[2], k2[2], k3[2], k4[2], f_val[2], tmp[2];
    double x;
    y[0] = y0[0];
    y[1] = y0[1];
    x = x0;
    f(f_val,x,y);
    for(int i = 0; i < 2; i++){
        k1[i] = h*f_val[i];
    }
    for(int i = 0; i < 2; i++){
        tmp[i] = y[i]+k1[i]/2;
    }

    f(f_val,x+h/2,tmp);
    for(int i = 0; i < 2; i++){
        k2[i] = h*f_val[i];
    }
    for(int i = 0; i < 2; i++){
        tmp[i] = y[i]+k2[i]/2;
    }

    f(f_val,x+h/2,tmp);
    for(int i = 0; i < 2; i++){
        k3[i] = h*f_val[i];
    }
    for(int i = 0; i < 2; i++){
        tmp[i] = y[i]+k3[i];
    }

    f(f_val,x+h,tmp);
    for(int i = 0; i < 2; i++){
        k4[i] = h*f_val[i];
    }

    for(int i = 0; i < 2; i++){
        y[i] = y[i]+(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
    }
}

int RKsolver(double y_x1[2], void (*f)(double[2], double, double[2]), double x0, double y0[2], double x_1, double tol, FILE* plot, double x_T[2]){
    double x, x_, err, h_new, h;
    double y1[2], y[2], y1_[2], y_[2], y_s[2], d[2];
    double h_ = 0.1;// первоначальная длина шага
    double fac = 0.8; // коэфициент для вычисления новой длины шага
    double facmin = 0.001;// минимальная длина шага
    double facmax = 3;// максимальная длина шага
    int nT = 0;
    x = x0;
    y[0] = y0[0];
    y[1] = y0[1];
    if(plot != NULL){
        fprintf(plot,"%lf %lf\n",y[0],y[1]);
    }
    for(; x < x_1;){
        //printf("h_ = %lf x = %lf \n",h_,x);

        y_s[0] = y[0];
        y_s[1] = y[1];

        y_[0] = y[0];
        y_[1] = y[1];
        h = 2*h_;

        // шаг длины 2*h_
        RKstep(y1,f,x,y,h);
        y[0] = y1[0];
        y[1] = y1[1];

        // два шага длины h_
        x_ = x;
        for(int s = 0; s < 2; s++, x_ += h_){
            RKstep(y1_,f,x_,y_,h_);
            y_[0] = y1_[0];
            y_[1] = y1_[1];
        }

        //	printf("h_ = %e\n");
        for(int i = 0; i < 2; i++){
            d[i] = fabs(y_[i]+(y_[i]-y[i])/(16-1));
        }

        err = fabs((y_[0]-y[0])/d[0]) > fabs((y_[1]-y[1])/d[1]) ? fabs((y_[0]-y[0])/d[0]) : fabs((y_[1]-y[1])/d[1]);
        err /= (16-1);
        h_new = fac*pow((tol/err),(double)1/(4+1));
        h_new = facmin > h_new ? facmin : h_new;
        h_new = facmax < h_new ? facmax : h_new;
        h_ = h_*h_new;
        if(err > tol){
            y[0] = y_s[0];
            y[1] = y_s[1];
        }
        else{
            if(!(x_ < x_1)){
                h = x_1-x;
                y[0]=y_s[0];
                y[1]=y_s[1];
                RKstep(y1,f,x,y,h);
                y[0] = y1[0];
                y[1] = y1[1];
                for(int i = 0; i < 2; i++){
                    y_x1[i] = y[i];
                }
                break;

            }
            y[0] = y_[0];
            y[1] = y_[1];
            if(x_T != NULL){
                if(y_s[1]*y[1] < 0 && ++nT == 2){
                    x_T[0] = x;
                    x_T[1] = x_;
                }
            }
            x = x_;
            if(plot != NULL){
                fprintf(plot,"%lf %lf\n",y[0],y[1]);
            }
        }

    }
    return 0;
}

void findCycle(double y_c[2], void (*f)(double[2], double, double[2]), double y0[2], double tol, double tolChord, double tolCycle){
    double x1,x2,x3;
    double x[2];
    double yn0[2], y1[2] , y2[2];
    yn0[0] = y0[0];
    yn0[1] = y0[1];
    for(;;){
        RKsolver(y1,f,0,yn0,200,tol,NULL,x);
        x1 = x[0]; x2 = x[1];
        do{
            RKsolver(y1,f,0,yn0, x1, tol, NULL, NULL);
            RKsolver(y2,f,0,yn0, x2, tol, NULL, NULL);
            x3 = x1 - (x2-x1)*y1[1]/(y2[1]-y1[1]);
            x1 = x2;
            x2 = x3;
        }while(fabs(x2-x1) > tolChord);
        if(fabs(yn0[0]-y2[0]) < tolCycle){
            y_c[0] = y2[0];
            y_c[1] = y2[1];
            break;
        }
        yn0[0] = y2[0];
    }
}


int main(){

    double T = 100;
    double t0 = 0.0;
    double x0 = 1;
    double z0 = 0;


    double tol = 1e-11;
    double tolChord = 1e-11;
    double tolCycle = 1e-9;
    double y0[2];
    double y1[2];
    double y_c[2];

    FILE* f1;

    y0[0] = x0;
    y0[1] = z0;

    f1 = fopen("plot","w");
    RKsolver(y1, f, t0, y0, T, tol,f1,NULL);
    fclose(f1);

    findCycle(y_c, f, y0, tol, tolChord, tolCycle);
    printf("%lf %lf\n", y_c[0], y_c[1]);

    f1 = fopen("plot_cycle","w");
    RKsolver(y1, f, 0, y_c, T, tol, f1, NULL);
    fclose(f1);

    return 0;
}
