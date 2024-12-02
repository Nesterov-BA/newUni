#ifndef MYMISTAKE
#define MYMISTAKE



void calc_mistake(double* res, int size)
{
    double result = 0;
    for (int i = 0; i < size; i++)
        res[i*size + i]--;
    for(int i = 0; i < size; i++)
    {
        for(int j = 0; j < size; j++)
        {
            result += dabs(res[i*size + j]);
        }
    }
    printf("Невязка = %f\n", result);
}


#endif