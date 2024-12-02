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
void gauss(double** dr, double *b){
    double mx;
    int m = 2;
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
