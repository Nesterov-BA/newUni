#include <iostream>
#include <fstream>
#include <ctime>
#include <math.h>

const int m = 3; //Количество параметров пристрелки
double alpha;
double eps = 1e-12;
double dr[m][m];

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
        std::cout << "___________________" << std::endl;
        for (int i =0; i<m;i++){
            for (int j=0; j<m; j++){
                std::cout << dr[i][j] << " ";
            }
            std::cout << " " << b[i] << std::endl;
        }
        std::cout << "___________________" << std::endl;

        mx = dr[i][i];
        dr[i][i] = 1;
        b[i]=b[i]/mx;
        for (int j=i+1; j<m;j++){
            dr[i][j]=dr[i][j]/mx;
        }

        std::cout << "___________________" << std::endl;
        for (int i =0; i<m;i++){
            for (int j=0; j<m; j++){
                std::cout << dr[i][j] << " ";
            }
            std::cout << " " << b[i] << std::endl;
        }
        std::cout << "___________________" << std::endl;

        for (int j=0; j<m; j++) {
			if(j != i){
				for (int l=i+1; l<m; l++){
                    dr[j][l] -= dr[j][i]*dr[i][l];
                }
                b[j] -= dr[j][i]*b[i];
				dr[j][i] = 0;
			}
        }
        std::cout << "___________________" << std::endl;
        for (int i =0; i<m;i++){
            for (int j=0; j<m; j++){
                std::cout << dr[i][j] << " ";
            }
            std::cout << " " << b[i] << std::endl;
        }
        std::cout << "___________________" << std::endl;
    }
}

int main(){
    double b[m];
    b[0] = 2;
    b[1] = 2;
    b[2] = 4;
    for (int i =0; i<m;i++){
        for (int j=0; j<m; j++){
            dr[i][j]=i*m+j+1;
        }
    }
    dr[2][2] = 10;
    std::cout << "___________________" << std::endl;
    for (int i =0; i<m;i++){
        for (int j=0; j<m; j++){
            std::cout << dr[i][j] << " ";
        }
        std::cout << " " << b[i] << std::endl;
    }
    std::cout << "___________________" << std::endl;
    gauss(b);
    std::cout << "___________________" << std::endl;
    for (int i =0; i<m;i++){
        for (int j=0; j<m; j++){
            std::cout << dr[i][j] << " ";
        }
        std::cout << " " << b[i] << std::endl;
    }
    std::cout << "___________________" << std::endl;


    return 1;
}