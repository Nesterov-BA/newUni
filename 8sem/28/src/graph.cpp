#include "graph.hpp"
#include <cstdio>

void graph(std::string filename, double func(double), double* bounds)
{
    int width = 300;
    double stepX = (bounds[1] - bounds[0])/width;
    FILE* graph = fopen(filename.c_str(), "w");

    fprintf(graph, "x1,val\n");
    for (int i = 0; i < width; i++)
    {
        fprintf(graph, "%lf,%lf\n", bounds[0] + i*stepX, func(bounds[0] + i*stepX));
        printf("%lf,%lf\n", bounds[0] + i*stepX, func(bounds[0] + i*stepX));
    }

}
