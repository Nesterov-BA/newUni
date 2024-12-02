#include "testFuncs.hpp"
#include <vector>

std::vector<double> x2y2(double x, double y)
{
    std::vector<double> result = {x*x, y*y};
    return result;
}
std::vector<double> expofx2y2(double x, double y)
{
    std::vector<double> result = {exp(x*x)-2, exp(y*y)-3};
    return result;
}
