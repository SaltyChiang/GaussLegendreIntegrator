#pragma once

#include <cmath>

double mu = 0.0;
double sigma = 1.0;

double Gaussian(double x)
{
    double a = (x - mu) / sigma;
    double pi = 3.14159265358979326846;
    return 1 / sqrt(2 * pi) / sigma * std::exp(-0.5 * a * a);
}

double Gaussian(double x, double mu, double sigma)
{
    ::mu = mu;
    ::sigma = sigma;
    return Gaussian(x);
}