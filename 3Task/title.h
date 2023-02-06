#ifndef TITLE_H
#define TITLE_H

#include <random>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <thread>
#include <mutex>
#include <ctime>
#include <cmath>

inline constexpr int ThrNum { 8 };
bool is_double(const std::string s);
double f(double i, double j, double n);

class Vector;

class Matrix
{
private:
    double **a;
    int n;
public:
    Matrix();
    Matrix(const Matrix& other);
    Matrix(Matrix&& other);
    Matrix(int N, double E);
    Matrix(int k, Matrix &u);
    Matrix( Matrix &u, int k);
    ~Matrix();
    void Delete();
    const int Size();
    void ReadFromFile(std::string file, int N);
    void Gen(int N);
    void Print();
    void ToFile(std::string file);
    double Trace();
    double Len();
    double TraceS();
    Vector diag();
    Matrix Inverse();
    double* operator [] (int const &index);
    const double* operator [] (int const &index) const;
    Matrix& operator=(const Matrix& right);
    Matrix& operator=(Matrix&& other);
    Matrix operator+( Matrix const &other );
    Matrix operator-( Matrix const &other );
    Matrix operator*( Matrix const &other );
    Vector GetSpRow(int N, int M, bool fl);
    double Norm();
    double PrNorm();
    Matrix T();
};

class Vector
{
private:
    double *a;
    int n;
public:
    Vector();
    Vector(const Vector& other);
    Vector(Vector&& other);
    Vector(int N);
    Vector(int N, int m);
    Vector(double *arr, int N);
    ~Vector();
    double Norm();
    double Sum();
    void Delete();
    Vector& operator=(const Vector& right);
    Vector& operator=(Vector&& other);
    double& operator [] (int const &index);
    const double& operator [] (int const &index) const;
    const int Size() const;
    Vector operator+( Vector const &other );
    Vector operator-( Vector const &other );
    Vector operator*( double const &num ); 
    Vector operator/( double const &num );
    Matrix operator*( Vector const &other );
    void Print();
    
};

Vector operator*( double const left, Vector const right);
Matrix Triangle(Matrix &a);
void Gauss(Matrix &tmp, Matrix &inv, int k, int r);
void GaussInv(Matrix &tmp, Matrix &inv, int k, int r);
double NormThr(Matrix &a, double &max, int r);
Vector QR(Matrix a, double norm, double eps);


#endif