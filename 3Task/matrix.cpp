#include "title.h"

bool is_double(const std::string s) {
    double d;
    return((std::istringstream(s) >> d >> std::ws).eof());
}

double f(double i, double j, double n){
    if (i == j) {
        if (i == 0) {
            return -1;
        }
        else if (i == n-1) return -((n-1)/n);
        else return -2;
    }
    else if ((j == i+1) || (j == i-1)) return 1;
    else return 0;
}

Matrix::Matrix(){
    a = NULL;
    n = 0;
}

Matrix::Matrix( const Matrix& other){
    n = other.n;
    a = new double *[n];
    for (int i = 0; i < n; i++) {
        a[i] = new double[n];
        for (int j = 0; j < n; j++) {
            a[i][j] = other[i][j];
        }
    }
}

Matrix::Matrix(Matrix&& other)
   : a(nullptr)
   , n(0)
{
    a = other.a;
    n = other.n;
    other.a = nullptr;
    other.n = 0;
}

Matrix::Matrix(int N, double E = 0) {
    n = N;
    a = new double *[n];
    for (int i = 0; i < n; i++) {
        a[i] = new double[n];
        for (int j = 0; j < n; j++) {
            a[i][j] = (i == j) * E;
        }
    }
}

Matrix::~Matrix(){
    this->Delete();
}

void Matrix::Delete(){
    if (a != nullptr) {
        if (*a != nullptr) {
            for (int i = 0; i < n; i++) {
                delete[] a[i];
            }
        }
        n = 0;
        delete[] a;
    }
}

const int Matrix::Size(){
    return n;
}

double* Matrix::operator [] (int const &index){
    if (!((index < n) && (index >= 0))) throw "Err0";
    return a[index];
}

const double* Matrix::operator [] (int const &index) const{
    if (!((index < n) && (index >= 0))) throw "Err1";
    return a[index];
}

Matrix& Matrix::operator=(const Matrix& right) {
    if (this == &right) {
        return *this;
    }
    this->Delete();
    n = right.n;
    a = new double *[n];
    for (int i = 0; i < n; i++) {
        a[i] = new double[n];
        for (int j = 0; j < n; j++) {
            a[i][j] = right[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator=(Matrix&& other){
    if (this != &other){
        this->Delete();
        n = other.n;
        a = other.a;
        other.a = nullptr;
        other.n = 0;
    }
    return *this;
}

Matrix Matrix::operator+( Matrix const &other ){
    if (n < other.n) {
        throw "Err: different sizes";
    }
    Matrix res(n);
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            if ((i<other.n) && (j<other.n)) res[i][j] = a[i][j] + other[i][j];
            else res[i][j]=a[i][j];
        }
    }
    return res;
}

Matrix Matrix::operator-( Matrix const &other ){
    if (n < other.n) {
        throw "Err: different sizes";
    }
    Matrix res(n);
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            if ((i<other.n) && (j<other.n)) res[i][j] = a[i][j] - other[i][j];
            else res[i][j]=a[i][j];
        }
    }
    return res;
    return res;
}

Matrix Matrix::operator*( Matrix const &other ){
    int max = 0;
    int min = 0;
    double **uk = nullptr;
    if (n >= other.n) {
        max = n;
        min = other.n;
        uk = this->a;
    }
    else {
        max = other.n;
        min = n;
        uk = other.a;
    }
    Matrix res(max);
    for(int i=0; i<max; i++){
		for(int j=0; j<max; j++){
            if ((i<min) && (j<min)) {
                for(int k=0; k<min; k++){
                    res[i][j]+=a[i][k]*other[k][j];
                }
            }
            else res[i][j] = uk[i][j];
		}
	}
    uk = nullptr;
	return res;
}

std::mutex g_lock;

//Матричная норма ||A||1 подчиненная векторной норме ||x||=sum(|xi|)
double Matrix::Norm(){
    double max=0;
    std::thread thr[ThrNum];
    for (int num = 0; num < ThrNum; num++){
        thr[num] = std::thread(NormThr, std::ref(*this), std::ref(max), num);
    }
    for (int num = 0; num < ThrNum; num++){
        thr[num].join();
    }
    return max;
}

double Matrix::PrNorm(){
    double max=0;
    double sum;
    for (int j=0; j<n; j++){
        sum = 0;
        for (int i=0; i<n; i++){
            sum+=fabs(a[i][j]);
        }
        if (sum > max) max = sum;
    }
    return max;
}

double NormThr(Matrix &a, double &max, int r){
    double sum;
    int n = a.Size();
    while (r<n){
        sum = 0;
        for (int j=0; j<n; j++){
            sum+=fabs(a[r][j]);
        }
        g_lock.lock();
        if (sum > max) max = sum;
        g_lock.unlock();
        r+=ThrNum;
    }
}

void Matrix::ReadFromFile(std::string file, int N){
    this->Delete();
    n = N;
    a = new double*[n];
    for (int i = 0; i<n; i++) a[i] = new double[n];
    std::ifstream in(file);
	if (in.is_open()) {
		std::string str;
		std::string tmp;
		int i = 0;
		while (std::getline(in, str)){
            //std::cout << "Str" << i << ": " << str << std::endl;
			if (i == n) throw "Wrong format: too much data in row";
			std::istringstream ss(str);
            int j = 0;
			while (ss >> tmp)
			{
                //std::cout << "tmp" << j << "=" << tmp << " ";
				if (j == n) throw "Wrong format: too much data in line";
				if (is_double(tmp)) {
					a[i][j] = std::stod(tmp);
				}
				else throw "Wrong format: not double";
				j++;
			}
            //std::cout << std::endl;
			if (j != n) {
				throw "Wrong format: not enough data in line";
			}
			i++;
			
		}
        std::cout << std::endl;
		if (i != n) throw "Wrong format: not enough data in row";
	}
	else throw "Can\'t open file";
}

void Matrix::Print(){
    std::cout << std::fixed;
    std::cout.precision(3);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << std::setw(7) << a[i][j];
        }
        std::cout << std::endl;
    }
}

void Matrix::ToFile(std::string file){
    std::ofstream fout(file);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fout << std::setw(4) << a[i][j];
        }
        fout << std::endl;
    }
    fout.close();
}

void Matrix::Gen(int N){
    this->Delete();
    n = N;
    a = new double *[n];
    for (int i=0; i<n; i++) {
        a[i] = new double[n];
        for (int j = 0; j < n; j++) {
            a[i][j] = f(i, j, n);
        }
    }
}

//Метода Гаусса без выбора главного элемента
Matrix Matrix::Inverse(){
    Matrix inv(n, 1);
    Matrix tmp;
    tmp = *this;
    double dtmp;
    for (int i = 0; i < n; i++) {
        if (tmp[i][i] == 0) throw "Error: diagonal element eqals zero";
        dtmp = tmp[i][i];
        for (int j = 0; j < n; j++) {
            if (tmp[i][j] != 0) tmp[i][j] = tmp[i][j]/dtmp;
            if (inv[i][j] != 0) inv[i][j] = inv[i][j]/dtmp;
        }
        for (int k = i+1; k < n; k++) {
            dtmp = tmp[k][i];
            if (dtmp != 0) {
                for (int j = 0; j < n; j++){
                    tmp[k][j] = tmp[k][j] - tmp[i][j]*dtmp;
                    inv[k][j] = inv[k][j] - inv[i][j]*dtmp;
                }
            }
        }
    } 
    for (int i=n-1; i > 0; i--){
        for (int j=i-1; j>=0; j--){
            dtmp = tmp[j][i];
            if (dtmp != 0) {
                for (int k=n-1; k>=0; k--) {
                    tmp[j][k] = tmp[j][k] - tmp[i][k]*dtmp;
                    inv[j][k] = inv[j][k] - inv[i][k]*dtmp; 
                }
            }
        }
    }
    return inv;
}

void Gauss(Matrix &tmp, Matrix &inv, int key, int r){
    int n = tmp.Size();
    int num = 2;
    int i= key+r;
    if (i == key) i+=num;
    while(i < n){
        double dtmp = tmp[i][key];
        if (dtmp != 0) {
            for (int j = 0; j < n; j++){
                tmp[i][j] = tmp[i][j] - tmp[key][j]*dtmp;
                inv[i][j] = inv[i][j] - inv[key][j]*dtmp;
            }
        }
        i+=num;
    }
}

void GaussInv(Matrix &tmp, Matrix &inv, int key, int r){
    int n = tmp.Size();
    int num = 2;
    int i = key-r;
    if (i == key) i-=num;
    while(i >= 0){
        double dtmp = tmp[i][key];
        if (dtmp != 0) {
            for (int j = 0; j < n; j++){
                tmp[i][j] = tmp[i][j] - tmp[key][j]*dtmp;
                inv[i][j] = inv[i][j] - inv[key][j]*dtmp;
            }
        }
        i-=num;
    }
}

Matrix Matrix::T(){
    Matrix mat(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mat[i][j] = a[j][i];
        }
    }
    return mat;
}

Vector Matrix::GetSpRow(int N, int M, bool fl = 1){
    if (M >= n) throw "Err: M>=n";
    Vector vec(n-M);
    for (int i=0; i<n-M; i++) {
        vec[i] = a[i+M*fl][N];
    }
    return vec;
}

Matrix::Matrix(int k, Matrix &u){
    n = u.n+k;
    a = new double *[n];
    for (int i = 0; i < n; i++) {
        a[i] = new double[n];
        for (int j = 0; j < n; j++) {
            if ((i<k) || (j<k)) a[i][j] = (i == j);
            else a[i][j] = u[i-k][j-k];
        }
    }
}

Matrix::Matrix(Matrix &u, int k){
    n = u.n+k;
    a = new double *[n];
    for (int i = 0; i < n; i++) {
        a[i] = new double[n];
        for (int j = 0; j < n; j++) {
            if ((i<u.n) && (j<u.n)) a[i][j] = u[i][j];
            else a[i][j] = (i == j);
        }
    }
}

double Matrix::Trace(){
    double sum = 0;
    for(int i=0; i<n; i++){
        sum+=a[i][i];
    }
    return sum;
}

double Matrix::Len(){
    double sum = 0;
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            sum+=a[i][j]*a[i][j];
        }
    }
    return sqrt(sum);
}

double Matrix::TraceS(){
    double sum = 0;
    for(int i=0; i<n; i++){
        sum+=a[i][i]*a[i][i];
    }
    return sqrt(sum);
}

Vector Matrix::diag(){
    Vector res(n);
    for (int i=0; i<n; i++){
        res[i] = a[i][i];
    }
    return res;
}

Matrix Triangle(Matrix &a){
    int n = a.Size();
    Matrix Q(n, 1);
    for (int i=0; i<n-2; i++){
        Vector tmp = a.GetSpRow(i, i+1);
        //std::cout << "tmp: ";
        //tmp.Print();
        Vector e(n-i-1, 0);
        Vector x = tmp-(e*tmp.Norm());
        //std::cout << "x = ";
        //x.Print();
        if (x.Norm() != 0) x = x/(x.Norm());
        Matrix E(n-i-1, 1);
        Matrix u = E - x*(x*2);
        Matrix u1(i+1, u);
        //std::cout << "Matrix u1: " << std::endl;
        //u1.Print();
        a = (u1*a)*u1;
        Q = Q*u1;
        //std::cout << "Matrix a: " << std::endl;
        //a.Print();
    }
    return Q;
}

Vector QR(Matrix a, double norm, double eps= 1e-12){
    int n = a.Size();
    for (int k=0; k<n-2; k++){
        double pr =0;
        while ((fabs(a[n-k-1][n-k-2]) > eps) && (fabs(fabs(a[n-k-1][n-k-2])-fabs(pr)) > eps)) {
            pr = a[n-k-1][n-k-2];
            Matrix Q(n-k, 1);
            double sk = 1;
            if (fabs(a[n-k-1][n-k-1]) > eps)  sk =a[n-k-1][n-k-1];
            Matrix Ek(n-k, sk);
            a = a - Ek; 
            for (int i=0; i<n-k-1; i++){
                Vector x(n-k-i);
                x[0] = a[i][i]-sqrt(a[i][i]*a[i][i]+a[i+1][i]*a[i+1][i]);
                x[1] = a[i+1][i];
                double t = sqrt(x[0]*x[0]+x[1]*x[1]);
                if (t!=0){
                    x[0] = x[0]/t;
                    x[1] = x[1]/t;
                }
                Matrix E(n-k-i, 1);
                Matrix u = E - x*(x*2);
                Matrix u1(i, u);
                //std::cout << "Matrix u1: " << std::endl;
                //u1.Print();
                //std::cout << "Matrix a: " << std::endl;
                //a.Print();
                a = u1*a;
                //std::cout << "Matrix a (after): " << std::endl;
                //a.Print();
                Q = Q*u1;
            }
            a = a*Q+Ek;
        }
        /*
        std::cout << "fabs(a[n-k-1][n-k-2]-pr: " << fabs(fabs(a[n-k-1][n-k-2])-fabs(pr)) << std::endl;
        std::cout << "a[1][2]: " << a[n-k-1][n-k-2] << std::endl;
        std::cout << "eps*norm: " << eps*norm << std::endl;
        std::cout << "lamm" << n-k << ": " << a[n-k-1][n-k-1] <<std::endl;
        */
    }
    //std::cout << "Matrix A:" << std::endl;
	//a.Print();
    Vector res = a.diag();
    if (a[0][1] > eps){
        double d = (a[0][0]-a[1][1])*(a[0][0]-a[1][1])+4*(a[0][1]*a[1][0]);
        if (d < 0) throw "Complex root";
        res[0]=(a[0][0]+a[1][1]+sqrt(d))/2;
        res[1]=(a[0][0]+a[1][1]-sqrt(d))/2;
    }
    return res;
}
