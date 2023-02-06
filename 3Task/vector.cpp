#include "title.h"

Vector::Vector(){
    a = NULL;
    n = 0;
}

Vector::Vector(const Vector& other){
    n = other.n;
    a = new double [n];
    for (int i = 0; i < n; i++) {
        a[i] = other[i];
    }
}

Vector::Vector(Vector&& other)
   : a(nullptr)
   , n(0)
{
    a = other.a;
    n = other.n;
    other.a = nullptr;
    other.n = 0;
}

Vector::Vector(int N){
    n = N;
    a = new double [n];
    for(int i=0; i<n; i++){
        a[i]=0;
    }
}

Vector::Vector(int N, int m){
    n = N;
    a = new double [n];
    for(int i=0; i<n; i++){
        a[i]=0;
    }
    a[m] = 1;
}

Vector::Vector(double *arr, int N){
    n = N;
    a = new double [n];
    for(int i=0; i<n; i++){
        a[i]=arr[i];
    }
}

Vector::~Vector(){
    delete[] a;
    n = 0;
}

double Vector::Norm(){
    double sum = 0;
    for(int i=0; i<n; i++){
        sum+= (a[i]*a[i]);
    }
    return sqrt(sum);
}

double Vector::Sum(){
    double sum = 0;
    for(int i=0; i<n; i++){
        sum+= a[i];
    }
    return sum;
}

void Vector::Delete(){
    delete[] a;
    n = 0;
}

Vector& Vector::operator=(const Vector& right) {
    if (this == &right) {
        return *this;
    }
    this->Delete();
    n = right.n;
    a = new double [n];
    for (int i = 0; i < n; i++) {
        a[i] = right[i];
    }
    return *this;
}

Vector& Vector::operator=(Vector&& other){
    if (this != &other){
        this->Delete();
        n = other.n;
        a = other.a;
        other.a = nullptr;
        other.n = 0;
    }
    return *this;
}

double & Vector::operator [] (int const &index){
    if (!((index < n) && (index >= 0))) throw "Err0";
    return (a[index]);
}

const double& Vector::operator [] (int const &index) const{
    if (!((index < n) && (index >= 0))) throw "Err1";
    return a[index];
}

const int Vector::Size() const{
    return n;
}

Vector Vector::operator+( Vector const &other ){
    if (n != other.n) {
        throw "Err: different vec sizes";
    }
    Vector res(n);
    for (int i=0; i<n; i++){
        res[i] = a[i] + other[i];
    }
    return res;
}

Vector Vector::operator-( Vector const &other ){
    if (n != other.n) {
        throw "Err: different vec sizes";
    }
    Vector res(n);
    for (int i=0; i<n; i++){
        res[i] = a[i] - other[i];
    }
    return res;
}

Vector Vector::operator*( double const &num ){
    Vector res(n);
    for (int i=0; i<n; i++){
        res[i] = a[i]*num;
    }
    return res;
}

Vector operator*( double const left, Vector const right){
    Vector res(right.Size());
    for (int i=0; i<right.Size(); i++){
        res[i] = right[i]*left;
    }
    return res;
}

Vector Vector::operator/( double const &num ){
    Vector res(n);
    if (num == 0) {
        for (int i=0; i<n; i++){
            res[i] = 0;
        }
    }
    else{
    for (int i=0; i<n; i++){
        res[i] = a[i]/num;
    }
    }
    return res;
}

Matrix Vector::operator*( Vector const &other ){
    if (n != other.n) {
        throw "Err: different vec sizes";
    }
    Matrix mat(n, 0);
    for (int i=0; i<n; i++){
        for (int j =0; j<n; j++){
            mat[i][j] = a[i]*other[j] ;
        }
    }
    return mat;
}

void Vector::Print(){
    for(int i=0; i<n; i++){
        std::cout << a[i] << " ";
    }
    std::cout << std::endl; 
}
