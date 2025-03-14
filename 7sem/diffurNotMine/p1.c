#include<stdio.h>
#include<math.h>


void f(double f_val[2], double x, double y[2]);
void f(double f_val[2], double x, double y[2]){
    double a = 0.1;
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
    /*
    if (x0 < 1e-10)
    {
        printf("step = %f", h);
        printf("\n %.10f, %.10f", k1[0]/h, k1[1]/h);
        printf("\n %.10f, %.10f", k2[0]/h, k2[1]/h);
        printf("\n %.10f, %.10f", k3[0]/h, k3[1]/h);
        printf("\n %.10f, %.10f\n", k4[0]/h, k4[1]/h);
        printf("(%.10f, %.10f) -> (%.10f, %.10f)\n", y0[0], y0[1], y[0], y[1]);

    }
    */
}

int RKsolver(double y_x1[2], void (*f)(double[2], double, double[2]), double x0, double y0[2], double x_1, double tol, FILE* plot, double x_T[2]){
    double x, x_, err, h_new, h;
    double y1[2], y[2], y1_[2], y_[2], y_s[2], d[2];
    //y_s - старт шага, который мы не трогаем, на него откатывеася если большая ошибка
    //y1_, y_ - временные значения конца шага
    // y - правильно посчитанное значение на конце шага, его записываем в файл
    // y1 - то, что мы запишем в y, реально ли он нужен не понятно
    double h_ = 0.1;// первоначальная длина шага
    double fac = 0.8; // коэфициент для вычисления новой длины шага
    double facmin = 0.001;// минимальная длина шага
    double facmax = 3;// максимальная длина шага
    int count = 0;
    int nT = 0;
    x = x0;
    y[0] = y0[0];
    y[1] = y0[1];
    if(plot != NULL){
        fprintf(plot,"%lf %lf %f\n",y[0],y[1], x);
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
       /* if(x < 1e-10)
        {
            printf("after double step = (%lf, %lf)\n",y[0], y[1]);
            printf("step = %f\n", h_);
        }
        */
        // в у записан конец шага
        //начало подсчета ошибки на шаге
        // два шага длины h_
        x_ = x;
        for(int s = 0; s < 2; s++, x_ += h_){
            RKstep(y1_,f,x_,y_,h_);
            y_[0] = y1_[0];
            y_[1] = y1_[1];
        }
     //   if(x < 1e-10)
       //     printf("after two steps = (%lf, %lf)\n",y_[0], y_[1]);
        // в у_ запмсан результат за 2 шага

        //	printf("h_ = %e\n");
        for(int i = 0; i < 2; i++){
            d[i] = fabs(y_[i]+(y_[i]-y[i])/(16-1));
        }

        err = fabs((y_[0]-y[0])/d[0]) > fabs((y_[1]-y[1])/d[1]) ? fabs((y_[0]-y[0])/d[0]) : fabs((y_[1]-y[1])/d[1]);
        err /= (16-1);
        h_new = fac*pow((tol/err),(double)1/(4+1));
        h_new = facmin > h_new ? facmin : h_new;
        h_new = facmax < h_new ? facmax : h_new;
        //if(x < 1e-10)
          //  printf("factor = %f\n\n",h_new);
        h_ = h_*h_new;
        //коне подсчета ошибки на шаге
        if(err > tol){
            //заново считаем с новым шагом
            y[0] = y_s[0];
            y[1] = y_s[1];
        }
        else{
            if(!(x_ < x_1))
            {
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
                // если дошли до конца по времени, то выходим из цикла
            }
            y[0] = y_[0];
            y[1] = y_[1];
            // обновляем значение решения, записывая результат за 2 шага
            if(x_T != NULL)
            {
                if(y_s[1]*y[1] < 0)
                    // printf("%f*%f,time = %f, count = %d\n",y_s[1],y[1],x, nT);
                if(y_s[1]*y[1] < 0 && ++nT == 2){
                    x_T[0] = x;
                    x_T[1] = x_;
                    // printf("Cycle time: %f, %f\n", x_T[0], x_T[1]);
                    // х - это время на начале шага
                    // х_ - это время в конце шага
                }
                // если надо записать время цикла , то записываем
            }
            x = x_;
            //обновляем время
            if(plot != NULL)
            {
                fprintf(plot,"%lf %lf %f\n",y[0],y[1], x);
            }
            //записываем в файл текущее значение y
        }

    }
    return 0;
}

void findCycle(double y_c[2], void (*f)(double[2], double, double[2]), double y0[2], double tol, double tolChord, double tolCycle){
    double x1,x2,x3;
    double x[2];
    double yn0[2], y1[2] , y2[2];
    int count = 0;
    // y1, y2 - временные массивы для хранения решений, из них гам надо только 2-е решение
    yn0[0] = y0[0];
    yn0[1] = y0[1];
    // записали, откуда стартуем
    for(;;){
        count++;
        RKsolver(y1,f,0,yn0,200,tol,NULL,x);
        //прогоняем решение на временном интервале 0-200
        x1 = x[0]; x2 = x[1];
        int count2 = 0;
        do{
            count2++;
            RKsolver(y1,f,0,yn0, x1, tol, NULL, NULL);
            RKsolver(y2,f,0,yn0, x2, tol, NULL, NULL);
            x3 = x1 - (x2-x1)*y1[1]/(y2[1]-y1[1]);
            x1 = x2;
            x2 = x3;

        }while(fabs(x2-x1) > tolChord);
        printf("Count = %d, count2 = %d \n", count, count2);
        //printf("Count2 = %d\n", count2);
        // printf("x1 = %f, x2 = %f, x2 - x1 = %.20f\n", x1, x2, x2 - x1);
        // printf("yn[0] = %f, y[2] = %f, yn[0] - y[2] = %.20f\n", yn0[0], y2[0], yn0[0] - y2[0]);
        if(fabs(yn0[0]-y2[0]) < tolCycle){
            //старт совпадает с концом, все отлично
            y_c[0] = y2[0];
            y_c[1] = y2[1];
            break;
        }
        yn0[0] = y2[0];
    }
}


int main(){

    double T = 10;
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

    f1 = fopen("plot.csv","w");
    RKsolver(y1, f, t0, y0, T, tol,f1,NULL);
    fclose(f1);

     findCycle(y_c, f, y0, tol, tolChord, tolCycle);
     printf("%lf %lf\n", y_c[0], y_c[1]);

     f1 = fopen("plot_cycle.csv","w");
     RKsolver(y1, f, 0, y_c, T, tol, f1, NULL);
     fclose(f1);

    return 0;
}
