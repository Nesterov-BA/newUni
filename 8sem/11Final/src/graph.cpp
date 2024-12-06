#include "graph.hpp"
#include <cstdio>

void graph(std::string filename, double func(double, double), double* bounds)
{
    int width = 30;
    int length = 30;
    double stepX = (bounds[1] - bounds[0])/width;
    double stepY = (bounds[3] - bounds[2])/length;
    FILE* graph = fopen(filename.c_str(), "w");

    fprintf(graph, "x1,x2,val\n");
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j<length; j++)
        {
            fprintf(graph, "%lf,%lf,%lf\n", bounds[0] + i*stepX, bounds[2] + j*stepY, func(bounds[0] + i*stepX, bounds[2] + j*stepY));
            printf("%lf,%lf,%lf\n", bounds[0] + i*stepX, bounds[2] + j*stepY, func(bounds[0] + i*stepX, bounds[2] + j*stepY));
        }
    }

}
