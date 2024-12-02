#include "errors.hpp"

std::vector<double> error(double p2, double x1, function* functions)
{
    double start[4] = {p2, 0, 0, x1};
    double end[4];
    std::vector<double> res(2);
//    double finish = M_PI/2;
    solutionUpToTime(start, end, functions, finish);
    res[0] = end[2];
    res[1] = end[3] + M_PI/2;
    return res;
}
