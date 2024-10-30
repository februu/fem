#include "../include/quad.h"

double rect1D(double (*func)(double), int N)
{
    double step_size = 2 / (double)N;
    double sum = 0.0;
    for (int i = 0; i < N; i++)
        sum += func(-1 + i * step_size) * step_size;
    return sum;
}

double gauss1D(double (*func)(double), int N)
{
    double sum = 0.0;
    for (int i = 0; i < N; i++)
        sum += func(gaussNodes[N - 1][i]) * gaussWeights[N - 1][i];
    return sum;
}

double gauss2D(double (*func)(double, double), int N)
{
    double sum = 0.0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            sum += func(gaussNodes[N - 1][i], gaussNodes[N - 1][j]) * gaussWeights[N - 1][i] * gaussWeights[N - 1][j];
    return sum;
}
