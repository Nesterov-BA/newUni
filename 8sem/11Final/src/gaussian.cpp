#include <cmath>
#include <cstdio>
#include <iostream>
#include "gaussian.hpp"

void gaussian(double** matrix, double* vec)
{
    double maxVal;
    int maxIndex;
    maxVal = fabs(matrix[1][0]) > fabs(matrix[0][0]) ? matrix[1][0] : matrix[0][0];
    maxIndex = fabs(matrix[1][0]) > fabs(matrix[0][0]) ? 1 : 0;
    matrix[maxIndex][0] /= maxVal;
    matrix[maxIndex][1] /= maxVal;
    vec[maxIndex] /= maxVal;
    matrix[1 - maxIndex][1] -= matrix[1-maxIndex][0];
    vec[1-maxIndex] -= matrix[1-maxIndex][0]*vec[maxIndex];
    vec[1 - maxIndex] /= matrix[1-maxIndex][1];
    vec[maxIndex] -= matrix[maxIndex][1]*vec[1 - maxIndex];
    if (maxIndex == 1)
    {
        double temp = vec[1];
        vec[1] = vec[0];
        vec[0] = temp;
    }
}
void gauss(double** matrix, double *vector)
{
    int size = 2;
    double max = 0;
    int maxRow;
    double temp, coefficent;
    //printf("Initial matrix:\n");
    //printMatrix(matrix, size);
    //printf("Vector: %lf, %lf, %lf\n", vector[0], vector[1], vector[2]);
    // Приведение к верхнетреугольному виду
    for(int col = 0; col < size; col++)
    {
        max = 0;
        maxRow = col;
      //  printf("Col = %d\n", col);
        for(int row = col; row < size; row++)
        {
            if(fabs(matrix[row][col]) > max)
            {
                max = matrix[row][col];
                maxRow = row;
            }
        }
     //   printf("Max row = %d\n", maxRow);
        coefficent = 1/matrix[maxRow][col];
        swapRows(matrix, maxRow, col, size);
        rowTimesX(matrix, col, coefficent, size);
       // printMatrix(matrix, size);
        temp = vector[maxRow];
        vector[maxRow] = vector[col];
        vector[col] = temp*coefficent;
    //    printf("Vector: %lf, %lf, %lf\n", vector[0], vector[1], vector[2]);

        for(int i = col + 1; i < size; i++)
        {
            vector[i] -= vector[col]*matrix[i][col];
            rowPlusRowTimesX(matrix, size, i, col, -matrix[i][col]);
        }
     //   printMatrix(matrix, size);
     //   printf("Vector: %lf, %lf, %lf\n", vector[0], vector[1], vector[2]);
    }
    for(int col = size - 1; col > 0; col--)
    {
        for(int i = 0; i < col; i++)
        {
            vector[i] -= vector[col]*matrix[i][col];
            matrix[i][col] = 0;
        }
    }
   // printMatrix(matrix, size);
   // printf("Vector: %lf, %lf, %lf\n", vector[0], vector[1], vector[2]);
}

void swapRows(double** matrix, int i1, int i2, int size)
{
    double temp;
    for(int i = 0; i < size; i++)
    {
        temp = matrix[i1][i];
        matrix[i1][i] = matrix[i2][i];
        matrix[i2][i] = temp;
    }
}

void rowTimesX(double** matrix, int row, double x, int size)
{
    for(int i = 0; i < size; i++)
    {
        matrix[row][i] *= x;
    }
}

void rowPlusRowTimesX(double** matrix, int size, int row1, int row2, double x)
{
    for(int i = 0; i < size; i++)
        matrix[row1][i] += matrix[row2][i]*x;
}
