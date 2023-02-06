#include "title.h"


int main(void)
{
	std::string str;
	std::string fout;
	std::string f;
	bool fl = 0;
	int n;
	std::getline(std::cin, str);
	std::istringstream ss(str);
	ss >> fout;
	ss >> n;
	if (ss >> f) fl=1;
	try {
		Matrix A;
		if (fl) A.ReadFromFile(f, n);
		else A.Gen(n);
		int start_time =  clock();
		//Matrix B = A.Inverse();
		//B.ToFile(fout);
		//Matrix C = A*B;
		//Matrix E(n, 1);
		//double res = (C-E).Norm();
		double res = A.Norm();
		int end_time = clock();
		std::cout << "ThrNum: " << ThrNum << std::endl;
		std::cout << "Norm A*Inv(A)-E: " << res << std::endl;
		std::cout << "Time (with threads): " << end_time-start_time << std::endl << std::endl;
		start_time = clock();
		res = A.PrNorm();
		end_time = clock();
		std::cout << "Norm A*Inv(A)-E: " << res << std::endl;
		std::cout << "Time (without threads): " << end_time-start_time<< std::endl;;

	} 
	catch (const char* msg){
		std::cerr << msg << std::endl;
		return -1;
	}
	return 0;
}