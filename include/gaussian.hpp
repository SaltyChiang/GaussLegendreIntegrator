#pragma once

#include <cmath>

double Gaussian(double x)
{
    return Gaussian(x, 0.0, 1.0);
}

double Gaussian(double x, double mu, double sigma)
{
    double a = (x - mu) / sigma;
    double pi = 3.14159265358979326846;
    return 1 / sqrt(2 * pi) / sigma * std::exp(-0.5 * a * a);
}