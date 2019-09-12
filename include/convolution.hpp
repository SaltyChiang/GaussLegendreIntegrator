#pragma once

#include <fftw3.h>

double *f1in, *f2in;
fftw_complex *f1out, *f2out;
fftw_plan p1, p2, pp;

double a, b, x;
double (*function1)(double), (*function2)(double);
int num;

void SetConvolution(double a, double b, double function1(double), bool useFFT, int num);
void SetConvolution(double a, double b, double function1(double), double function2(double), bool useFFT, int num);
void PreConvolutionFFT(double a, double b);
void AfterConvolutionFFT();
double Convolution(double x);
double ConvolutionFFT(double x);