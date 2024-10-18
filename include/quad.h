#pragma once

#include <cmath>

const double gaussNodes[3][3] = {{0.0}, {-1 / sqrt(3), 1 / sqrt(3)}, {-sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0)}};
const double gaussWeights[3][3] = {{2.0}, {1.0, 1.0}, {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0}};

double rect1D(double (*func)(double), int N);
double gauss1D(double (*func)(double), int N);
double gauss2D(double (*func)(double, double), int N);