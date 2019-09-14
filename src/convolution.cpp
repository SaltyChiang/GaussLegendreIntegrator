#include <iostream>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include "convolution.hpp"
#include "integral.hpp"
#include "gaussian.hpp"

Convolution::Convolution(double a, double b, double function1(double), bool useFFT, int den) : Convolution(a, b, function1, Gaussian, useFFT, den)
{
}

Convolution::Convolution(double a, double b, double function1(double), double function2(double), bool useFFT, int den)
{
    this->a = a - (b - a) / 2;
    this->b = b + (b - a) / 2;
    ::function1 = function1;
    ::function2 = function2;
    this->den = den;
    this->num = (int)((double)den * (this->b - this->a));
    if (useFFT)
        PreConvolutionFFT(this->a, this->b);
}

double Convolution::DoConvolution(double x)
{
    ::x = x;
    GaussLegendreIntegrator integrator(num);
    return integrator.DoIntegral(a, b, [](double y) -> double { return ::function1(y) * ::function2(::x - y); });
}

double Convolution::DoConvolutionFFT(double x)
{
    int i;
    double xx = a;
    double yy = (b - a) / ((double)num - 1);
    for (i = 0; i < num; i++)
    {
        if (x >= xx && x < xx + yy)
            break;
        xx += yy;
    }
    i = (i + num / 2) % num;
    return cOut[(i + 1) % num] * (x - xx) / yy + cOut[i] * (xx + yy - x) / yy;
}

void Convolution::PreConvolutionFFT(double a, double b)
{
    fIn = (double *)fftw_malloc(sizeof(double) * num);
    fOut = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * num);
    gIn = (double *)fftw_malloc(sizeof(double) * num);
    gOut = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * num);
    cIn = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * num);
    cOut = (double *)fftw_malloc(sizeof(double) * num);
    double xx = a;
    double yy = (b - a) / ((double)num - 1);
    GaussLegendreIntegrator integrator(den);
    for (int i = 0; i < num; i++)
    {
        fIn[i] = integrator.DoIntegral(xx - yy / 2., xx + yy / 2., function1);
        gIn[i] = integrator.DoIntegral(xx - yy / 2., xx + yy / 2., function2);
        xx += yy;
    }
    p1 = fftw_plan_dft_r2c_1d(num, fIn, fOut, FFTW_ESTIMATE);
    p2 = fftw_plan_dft_r2c_1d(num, gIn, gOut, FFTW_ESTIMATE);
    pp = fftw_plan_dft_c2r_1d(num, cIn, cOut, FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_execute(p2);
    for (int i = 0; i < num; i++)
    {
        if (i <= (num / 2 + 1))
        {
            cIn[i][0] = (fOut[i][0] * gOut[i][0] - fOut[i][1] * gOut[i][1]);
            cIn[i][1] = (fOut[i][0] * gOut[i][1] + fOut[i][1] * gOut[i][0]);
        }
        // else
        // {
        //     cIn[i][0] = (fOut[num - i][0] * gOut[num - i][0] - fOut[num - i][1] * gOut[num - i][1]);
        //     cIn[i][1] = (-fOut[num - i][0] * gOut[num - i][1] - fOut[num - i][1] * gOut[num - i][0]);
        // }
    }
    fftw_execute(pp);
    for (int i = 0; i < num; i++) /*OUTPUT*/
    {
        cOut[i] /= num;
        cOut[i] /= yy;
        // printf("%9f,%9f  %9f,%9f  %9f,%9f\n", fOut[i][0], fOut[i][1], gOut[i][0], gOut[i][1], cIn[i][0], cOut[i]);
    }
}

Convolution::~Convolution()
{
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(pp);
    fftw_free(fIn);
    fftw_free(fOut);
    fftw_free(gIn);
    fftw_free(gOut);
    fftw_free(cIn);
    fftw_free(cOut);
}
