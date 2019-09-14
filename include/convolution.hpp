#ifndef CONVOLUTION_HPP_
#define CONVOLUTION_HPP_

#include <fftw3.h>

static double x;
static double (*function1)(double), (*function2)(double);

class Convolution
{
private:
    double *fIn, *gIn, *cOut;
    fftw_complex *fOut, *gOut, *cIn;
    fftw_plan p1, p2, pp;
    double a, b;
    int num;
    int den;

public:
    Convolution(double, double, double function1(double), bool, int);
    Convolution(double, double, double function1(double), double function2(double), bool, int);
    ~Convolution();
    void PreConvolutionFFT(double, double);
    double DoConvolution(double);
    double DoConvolutionFFT(double);
};

#endif