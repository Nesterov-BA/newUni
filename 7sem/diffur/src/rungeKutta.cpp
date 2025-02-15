#include "rungeKutta.hpp"

int maxNumberOfCycles = 5;
bool limitedCycles = true;
double tolerance = 1.e-11;
double maxStep = 3;
double minStep = 0.00001;
double factor = 0.8;
double* c = new double[3];
double* b = new double[3];
double* a = new double[6];
double* d = new double[4];



double a10 = 1./4.;
double a20 = -189./800.;
double a21 = 729./800.;
double a30 = 214./891.;
double a31 = 1./33.;
double a32 = 650./891.;


double c2 = 1./4.;
double c3 = 27./40.;
double c4 = 1.;

double d1 = 533./2106.;
double d2 = 0;
double d3 = 800./1053.;
double d4 = -1./78.;



void clear()
{
    ofstream file("data.txt");
    file<<"";
    file.close();
}

void solutionUpToTime(double xStart, double yStart, double f(double, double), double g(double, double), double finish, double* cycleTimeLess, double* cycleTimeMore, double* xEnd, double* yEnd)
{
    ofstream file("data.csv");
    ofstream fileErr("errorLog.txt");
    ofstream fileErrReg("errorRegular.txt");
    ofstream numberOfPoints("numberOfPoints.txt");
    int count = 0;
    int numOfPoints = 0;
    double time = 0;
    double step = 0.1;
    double tempX, tempY;
    double errorSum = 0;
    double globalError = 0;
    double globalErrorRegular = 0;
    file << "X,Y,time" << endl;
    file << xStart << "," << yStart << "," << time << endl;
    while(time < finish)
    {
        if(time + 2*step > finish)
        {
            step = (finish - time)/2;
            Runge_Kutta4ClassicSimple(xStart, yStart, step, f, g, &tempX, &tempY);
            Runge_Kutta4ClassicSimple(tempX, tempY, step, f, g, &tempX, &tempY);
            numOfPoints+=2;

        }
        else
        {
            Runge_Kutta4StepVariedSimple(xStart, yStart, f, g, &step, &tempX, &tempY, &errorSum, &globalError, &globalErrorRegular);
            numOfPoints++;
        }

        time += 2*step;

        if (yStart * tempY < 0)
        {
            //printf("%.10f * %.10f, timeL = %f, timeR = %f, count = %d, errorSum = %e\n", yStart, tempY, time-2*step, time, count, errorSum);
        }
        if (yStart * tempY < 0 && ++count == 2)
        {
            *cycleTimeLess = time - 2*step;
            *cycleTimeMore = time;
        }
        fileErr << globalError << " " << logNormCalc(xStart, yStart, tempX, tempY) << " " << time << endl;
        fileErrReg << globalErrorRegular << " " << regNormCalc(xStart, yStart, tempX, tempY) << " " << time << endl;
        xStart = tempX;
        yStart = tempY;
        file << xStart << "," << yStart << "," << time << endl;

    }
    *xEnd = xStart;
    *yEnd = yStart;
    //printf("Tolerance: %e Error: %e, numOfPoints: %d, globalError: %e, globalErrorReg: %e\n", tolerance, errorSum, numOfPoints, globalError, globalErrorRegular);
    numberOfPoints << numOfPoints;


}
void fasterFindCycle(double xStart, double yStart, double f(double, double), double g(double, double), double* xEnd, double* yEnd)
{
    xStart = 1;
    yStart = 0;
    double chordStart, chordEnd, chordTemp;
    double cycleTimeLess, cycleTimeMore;
    double yEnd1, yEnd2;
    double xEnd1, xEnd2;
    double xStartTemp, yStartTemp;
    int count = 0;
    chordEnd = 10;
    xStartTemp = xStart;
    yStartTemp = yStart;

    while(count < 100)
    {
        fasterSolutionUpToTime(xStartTemp, yStartTemp, f, g, 200, &cycleTimeLess, &cycleTimeMore, &xEnd1, &yEnd1);
        chordStart = cycleTimeLess;
        chordEnd = cycleTimeMore;
        int count2 = 0;
        while (fabs(chordEnd - chordStart) > tolerance)
        {
            count2++;
            fasterSolutionUpToTime(xStartTemp, yStartTemp, f, g, chordStart, &cycleTimeLess, &cycleTimeMore, &xEnd1, &yEnd1);
            fasterSolutionUpToTime(xStartTemp, yStartTemp, f, g, chordEnd, &cycleTimeLess, &cycleTimeMore, &xEnd2, &yEnd2);
            //printf("(%.10f, %.10f), (%.10f, %.10f)\n", xEnd1, yEnd1, xEnd2, yEnd2);
            chordTemp = chordStart - (chordEnd-chordStart)*yEnd1/(yEnd2-yEnd1);
            chordStart = chordEnd;
            chordEnd = chordTemp;
            //std::cout << "Cycle time: " << chordStart << " " << chordEnd << ", difference = " << -chordStart + chordEnd << endl;
        }
        printf("count = %d, count2 = %d\n", count, count2);
        //printf("\nCycle end \n");

        if(fabs(xStartTemp - xEnd2) < tolerance*100)
        {
            *xEnd = xEnd2;
            *yEnd = yEnd2;
            //cout << "Cycle found from " << xStartTemp << "," << yStartTemp << " to " << xEnd2 << "," << yEnd2 << endl;
            //std::cout << "Cycle time: " << chordStart << " " << chordEnd << ", difference = " << -chordStart + chordEnd << endl;

            fasterSolutionUpToTime(xStartTemp, yStartTemp, f, g, chordStart, &cycleTimeLess, &cycleTimeMore, &xEnd1, &yEnd1);

            break;
        }
    //    printf("Difference = %.11lf\n", fabs(xStartTemp - xEnd2));
        xStartTemp = xEnd2;
        count++;
    }
    printf("\nCycle end, final count = %d\n", count);

}
void findCycle(double xStart, double yStart, double f(double, double), double g(double, double), double* xEnd, double* yEnd)
{
    xStart = 1;
    yStart = 0;
    double chordStart, chordEnd, chordTemp;
    double cycleTimeLess, cycleTimeMore;
    double yEnd1, yEnd2;
    double xEnd1, xEnd2;
    double xStartTemp, yStartTemp;
    int count = 0;
    chordEnd = 10;
    xStartTemp = xStart;
    yStartTemp = yStart;

    while(count < 100)
    {
        solutionUpToTime(xStartTemp, yStartTemp, f, g, 200, &cycleTimeLess, &cycleTimeMore, &xEnd1, &yEnd1);
        chordStart = cycleTimeLess;
        chordEnd = cycleTimeMore;
        int count2 = 0;
        while (fabs(chordEnd - chordStart) > tolerance)
        {
            count2++;
            solutionUpToTime(xStartTemp, yStartTemp, f, g, chordStart, &cycleTimeLess, &cycleTimeMore, &xEnd1, &yEnd1);
            solutionUpToTime(xStartTemp, yStartTemp, f, g, chordEnd, &cycleTimeLess, &cycleTimeMore, &xEnd2, &yEnd2);
            //printf("(%.10f, %.10f), (%.10f, %.10f)\n", xEnd1, yEnd1, xEnd2, yEnd2);
            chordTemp = chordStart - (chordEnd-chordStart)*yEnd1/(yEnd2-yEnd1);
            chordStart = chordEnd;
            chordEnd = chordTemp;
            //std::cout << "Cycle time: " << chordStart << " " << chordEnd << ", difference = " << -chordStart + chordEnd << endl;
        }
        printf("count = %d, count2 = %d\n", count, count2);
        //printf("\nCycle end \n");

        if(fabs(xStartTemp - xEnd2) < tolerance*100)
        {
            *xEnd = xEnd2;
            *yEnd = yEnd2;
            //cout << "Cycle found from " << xStartTemp << "," << yStartTemp << " to " << xEnd2 << "," << yEnd2 << endl;
            //std::cout << "Cycle time: " << chordStart << " " << chordEnd << ", difference = " << -chordStart + chordEnd << endl;

            solutionUpToTime(xStartTemp, yStartTemp, f, g, chordStart, &cycleTimeLess, &cycleTimeMore, &xEnd1, &yEnd1);

            break;
        }
    //    printf("Difference = %.11lf\n", fabs(xStartTemp - xEnd2));
        xStartTemp = xEnd2;
        count++;
    }
    printf("\nCycle end, final count = %d\n", count);

}


void Runge_Kutta4ClassicSimple(double startX, double startY, double step, double f(double, double), double g(double, double), double* endX, double* endY)
{
    double k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y;

    k1x = f(startX, startY);
    k1y = g(startX, startY);

    k2x = f(startX + step/2*k1x, startY + step/2*k1y);
    k2y = g(startX + step/2*k1x, startY + step/2*k1y);
;
    k3x = f(startX + step/2*k2x, startY + step/2*k2y);
    k3y = g(startX + step/2*k2x, startY + step/2*k2y);

    k4x = f(startX + step*k3x, startY + step*k3y);
    k4y = g(startX + step*k3x, startY + step*k3y);

    *endX = startX + step/6*(k1x + 2*k2x + 2*k3x + k4x);
    *endY = startY + step/6*(k1y + 2*k2y + 2*k3y + k4y);

}

void Runge_Kutta4StepVariedSimple(double startX, double startY, double f(double, double), double g(double, double), double* step, double* endX, double* endY, double* errorSum, double* globalError, double* globalErrorRegular)
{
    double err = 1;
    double tempEndX, tempEndY;
    double stepFac;
    while (err > tolerance)
    {
        Runge_Kutta4ClassicSimple(startX, startY, 2*(*step), f, g, endX, endY);
        tempEndX = *endX;
        tempEndY = *endY;

        Runge_Kutta4ClassicSimple(startX, startY, *step, f, g, endX, endY);
        Runge_Kutta4ClassicSimple(*endX, *endY, *step, f, g, endX, endY);

        double dX = fabs(*endX+(*endX-tempEndX)/(16-1));
        double dY = fabs(*endY+(*endY-tempEndY)/(16-1));

        err = fabs((*endX-tempEndX)/dX) > fabs((*endY-tempEndY)/dY) ? fabs((*endX-tempEndX)/dX) : fabs((*endY-tempEndY)/dY);
        err /= (16-1);
        stepFac = factor*pow((tolerance/err),(double)1/(4+1));
        stepFac = minStep > stepFac ? minStep : stepFac;
        stepFac = maxStep < stepFac ? maxStep : stepFac;
        *step = *step*stepFac;
    }

    *globalError *= exp(*step*logNormCalc(startX, startY, *endX, *endY));
    *globalError += err;
    *globalErrorRegular *= exp(*step*regNormCalc(startX, startY, *endX, *endY));
    *globalErrorRegular += err;
    *errorSum += err;

}

void checkCycle(double f(double, double), double g(double, double))
{
    double xStart = 1;
    double yStart = 0;
    double xEnd, yEnd, cycleTimeLess, cycleTimeMore;
    double xEnd1, yEnd1, xEnd2, yEnd2;
    solutionUpToTime(xStart, yStart, f, g, 50, &cycleTimeLess, &cycleTimeMore, &xEnd, &yEnd);
    //printf("\ncycle1: %.10f, %.10f\n", cycleTimeLess, cycleTimeMore);
    double tStart = cycleTimeLess;
    double tEnd = cycleTimeMore;
    solutionUpToTime(xStart, yStart, f, g, tStart, &cycleTimeLess, &cycleTimeMore, &xEnd1, &yEnd1);
    solutionUpToTime(xStart, yStart, f, g, tEnd, &cycleTimeLess, &cycleTimeMore, &xEnd2, &yEnd2);
    //printf("ends: %.10f, %.10f\n", yEnd1, yEnd2);
}

int fasterSolutionUpToTime(double xStart, double yStart, double f(double, double), double g(double, double), double timeEnd, double* cycleTimeLess, double* cycleTimeMore, double* xEnd, double* yEnd){
    FILE* plot = fopen("plot.txt", "w");
    double time, tempTime, err, stepFactor, step;
    double xValue, yValue, xValueTemp, yValueTemp;
    double xValueTemp2Steps, yValueTemp2Steps, xStartTemp, yStartTemp;
    double d[2];
    //tempStartValues - старт шага, который мы не трогаем, на него откатывеася если большая ошибка
    //tempValues, tempValues2Steps - временные значения конца шага
    // values - правильно посчитанное значение на конце шага, его записываем в файл
    // y1 - то, что мы запишем в values, реально ли он нужен не понятно
    double stepIni = 0.1;// первоначальная длина шага
    double fac = 0.8; // коэфициент для вычисления новой длины шага
    double facmin = 0.001;// минимальная длина шага
    double facmax = 3;// максимальная длина шага
    int count = 0;
    int nT = 0;
    time = 0;
    xValue = xStart;
    yValue = yStart;
    if(plot != NULL){
        fprintf(plot,"%lf %lf %f\n",xValue,yValue, time);
    }
    while(time < timeEnd){
        //printf("stepIni = %lf x = %lf \n",stepIni,x);
        //printf("time = %.e\n",time);
        printf("stepIni = %.e x = %lf \n",stepIni,xValue);
        xStartTemp = xValue;
        yStartTemp = yValue;

        xValueTemp2Steps = xValue;
        yValueTemp2Steps = yValue;
        step = 2*stepIni;

        // шаг длины 2*stepIni
        Runge_Kutta4ClassicSimple(xStartTemp, yStartTemp, step, f, g, &xValue, &yValue);
       /* if(x < 1e-10)
        {
            printf("after double step = (%lf, %lf)\n",xValue, yValue);
            printf("step = %f\n", stepIni);
        }
        */
        // в у записан конец шага
        //начало подсчета ошибки на шаге
        // два шага длины stepIni
        tempTime = time;
        for(int s = 0; s < 2; s++, tempTime += stepIni){
            Runge_Kutta4ClassicSimple(xStartTemp, yStartTemp, stepIni, f, g, &xValueTemp2Steps, &yValueTemp2Steps);
        }
     //   if(x < 1e-10)
       //     printf("after two steps = (%lf, %lf)\n",xValueTemp2Steps, yValueTemp2Steps);
        // в у_ запмсан результат за 2 шага

        //	printf("stepIni = %e\n");
        
        d[0] = fabs(xValueTemp2Steps+(xValueTemp2Steps-xValue)/(16-1));
        d[1] = fabs(yValueTemp2Steps+(yValueTemp2Steps-yValue)/(16-1));
        

        err = fabs((xValueTemp2Steps-xValue)/d[0]) > fabs((yValueTemp2Steps-yValue)/d[1]) ? fabs((xValueTemp2Steps-xValue)/d[0]) : fabs((yValueTemp2Steps-yValue)/d[1]);
        err /= (16-1);
        stepFactor = fac*pow((tolerance/err),(double)1/(4+1));
        stepFactor = facmin > stepFactor ? facmin : stepFactor;
        stepFactor = facmax < stepFactor ? facmax : stepFactor;
        //if(x < 1e-10)
          //  printf("factor = %f\n\n",stepFactor);
        stepIni = stepIni*stepFactor;
        //коне подсчета ошибки на шаге
        printf("err = %e\n",err);
        if(err > tolerance){
            //заново считаем с новым шагом
            xValue = xStartTemp;
            yValue = yStartTemp;
        }
        else{
            printf("error is small\n");
            printf("err = %e\n",err);
            if(!(tempTime < timeEnd))
            {
                step = timeEnd-time;
                Runge_Kutta4ClassicSimple(xStartTemp, yStartTemp, step, f, g, xEnd, yEnd);
                break;
                // если дошли до конца по времени, то выходим из цикла
            }
            xValue = xValueTemp2Steps;
            yValue = yValueTemp2Steps;
            // обновляем значение решения, записывая результат за 2 шага
            if(yStartTemp*yValue < 0)
                 printf("%f*%f,time = %f, count = %d\n",yStartTemp,yValue,time, nT);
            if(yStartTemp*yValue < 0 && ++nT == 2){
                *cycleTimeLess = time;
                *cycleTimeMore = tempTime;
                // printf("Cycle time: %f, %f\n", cycleTimes[0], cycleTimes[1]);
                // х - это время на начале шага
                // х_ - это время в конце шага
            }
            // если надо записать время цикла , то записываем
            time = tempTime;
            //обновляем время
            if(plot != NULL)
            {
                fprintf(plot,"%lf %lf %f\n",xValue,yValue, time);
            }
            //записываем в файл текущее значение values
        }

    }
    return 0;
}
