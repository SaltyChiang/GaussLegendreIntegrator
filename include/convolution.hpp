#pragma once

namespace Convolution
{
    bool useFFT;
    double a, b, x;
    double (*function1)(double), (*function2)(double);

    double Convolution(double a, double b, double x, double function1(double), bool useFFT);
    double Convolution(double a, double b, double x, double function1(double), double function2(double), bool useFFT);
    double PreConvolution(double y);
    double DoConvolution();
}