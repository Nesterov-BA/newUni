#include "graph.hpp"
#include <cstdio>

void graph(std::string filename, double func(double, double), double* bounds)
{
    double step = 0.05;
    int width = (bounds[1] - bounds[0])/step;
    int length = (bounds[3] - bounds[2])/step;
    FILE* graph = fopen(filename.c_str(), "w");

    fprintf(graph, "x1,x2,val\n");
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j<length; j++)
        {
            fprintf(graph, "%lf,%lf,%lf\n", bounds[0] + i*step, bounds[2] + j*step, func(bounds[0] + i*step, bounds[2] + j*step));
        }
    }

}
