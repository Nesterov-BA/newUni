double scalarMultiplication(double* vec1, double* vec2, int size);


double scalarMultiplication(double* vec1, double* vec2, int size)
{
    double result = 0;
    
    for (int i = 0; i < size; i++)
    {
        result += vec1[i] * vec2[i];
    }

    return result;
}

