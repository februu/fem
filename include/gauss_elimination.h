#pragma once
#include <cmath>
#include <algorithm>

double *gaussElimination(double **a, double *b, int n)
{

    double *x = new double[n];
    for (int i = 0; i < n; i++)
        x[i] = 0.0;

    int *w = new int[n + 1];
    for (int i = 0; i < n + 1; i++)
        w[i] = i;

    // Gauss elimination
    for (int i = 0; i < n - 1; i++)
    {

        // Crout
        int k = i;
        for (int j = i + 1; j < n; j++)
            if (fabs(a[i][w[k]]) < fabs(a[i][w[j]]))
                k = j;
        std::swap(w[k], w[i]);

        // ====

        for (int k = i + 1; k < n; k++)
        {
            double tmp;
            if (w[i] != n)
                tmp = a[k][w[i]] / a[i][w[i]];
            else
                tmp = b[k] / b[i];

            for (int j = i + 1; j <= n; j++)
                if (w[j] != n)
                    a[k][w[j]] -= tmp * a[i][w[j]];
                else
                    b[k] -= tmp * b[i];
        }
    }

    // Solve for xn
    for (int i = n - 1; i >= 0; i--)
    {
        double s = b[i];
        for (int j = n - 1; j >= i + 1; j--)
        {
            if (w[j] != n)
                s -= a[i][w[j]] * x[w[j]];
            else
                s -= b[i] * x[w[j]];
        }

        if (w[i] != n)
            x[w[i]] = s / a[i][w[i]];
        else
            x[w[i]] = s / b[i];
    }

    delete[] w;
    return x;
}