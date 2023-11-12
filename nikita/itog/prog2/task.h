#ifndef __TASK_H_INCLUDED__
#define __TASK_H_INCLUDED__


// ЗАГОЛОВКИ ===================================

//void FindValues(int n, double* a, double* values, double eps, int* iterOut);
void FindValues(int n, double* diag, double* side, double* sideup, double* values, double eps, int* iterOut);
//void TriDiag(double* a, int n);

//void Shift(int n, double* a, int k, double s);
void Shift(double* diag, int k, double s);

//double CalcNorm(int n, double* a);
double CalcNorm(int n, double* diag, double* sign);

//void OldFindValues(int n, double* a, double* values, double eps, int* iterOut);
//void QR_without_arrays(int n, double* a, int k);
//void QR_without_arrays(int n, double* diag, double* sign, double* signup, int k);
void QR_without_arrays(double* diag, double* sign, double* signup, int k);
void Old_almost_triangular(int n, double* a);

void assign(double* diag, double* side, int i, int j, double value);
double get_value(double* diag, double* side, int i, int j);
void PrintTriDMatrix2(int n, double* diag, double* side);
void PrintVector2(int n, double* x);

void SaveMatrix(int n, double* a);
void PrintMatrix2(int n, double* a);

//void left_multiply_by_T_k_kpp(double* a, int n, int len, int k, double cos, double sin, double r);
//void left_multiply_by_T_k_kpp(double* diag, double* side, double* sideup, int n, int len, int k, double cos, double sin, double r);
void left_multiply_by_T_k_kpp(double* diag, double* side, double* sideup,  int k, double cos, double sin);
//void right_multiply_by_T_k_kpp(double* a, int n, int len, int k, double cos, double sin, double r);
//void right_multiply_by_T_k_kpp(double* diag, double* side, double* sideup, int n, int len, int k, double cos, double sin, double r);
void right_multiply_by_T_k_kpp(double* diag, double* side, double* sideup, int k, double cos, double sin);

void tridiag(double* a, int n);
void left_multiplication(double* a, int n, int k);
void right_multiplication(double* a, int n, int k);
void cleaning(double* a, int n, int k, double norm);
// =============================================


#endif /* not __TASK_H_INCLUDED__ */
