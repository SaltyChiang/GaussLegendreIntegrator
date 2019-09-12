#include <cmath>
#include <complex>
#include <fftw3.h>
#include "convolution.hpp"
#include "integral.hpp"
#include "gaussian.hpp"

void SetConvolution(double a, double b, double function1(double), bool useFFT, int num)
{
    SetConvolution(a, b, function1, Gaussian, useFFT, num);
}

void SetConvolution(double a, double b, double function1(double), double function2(double), bool useFFT, int num)
{
    ::a = a - (b - a) / 5;
    ::b = b + (b - a) / 5;
    ::function1 = function1;
    ::function2 = function2;
    ::num = num;
    if (useFFT)
        PreConvolutionFFT(::a, ::b);
}

double Convolution(double x)
{
    ::x = x;
    GaussLegendreIntegrator integrator(1000);
    return integrator.DoIntegral(a, b, [](double y) -> double { return ::function1(y) * ::function2(::x - y); });
}

double ConvolutionFFT(double x)
{
    int i;
    double xx = a;
    double yy = (b - a) / (num - 1);
    for (i = 0; i < num; i++)
    {
        if (x >= xx && x < xx + yy)
            break;
    }
    return f1in[i + 1] * (x - xx) / yy + f1in[i] * (xx + yy - x) / yy;
}

void PreConvolutionFFT(double a, double b)
{
    f1in = (double *)fftw_malloc(sizeof(double) * num);
    f1out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * num);
    f2in = (double *)fftw_malloc(sizeof(double) * num);
    f2out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * num);
    double xx = a;
    double yy = (b - a) / (num - 1);
    for (int i = 0; i < num; i++)
    {
        f1in[i] = function1(xx);
        f2in[i] = function2(xx);
        xx += yy;
    }
    p1 = fftw_plan_dft_r2c_1d(10, f1in, f1out, FFTW_ESTIMATE);
    p2 = fftw_plan_dft_r2c_1d(10, f2in, f2out, FFTW_ESTIMATE);
    pp = fftw_plan_dft_c2r_1d(10, f1out, f1in, FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_execute(p2);
    for (int i = 0; i < num; i++)
    {
        if (i <= (num / 2 + 1))
        {
            f1out[i][0] = f1out[i][0] * f2out[i][0] - f1out[i][1] * f2out[i][1];
            f1out[i][1] = f1out[i][0] * f2out[i][1] + f1out[i][1] * f2out[i][0];
        }
        else
        {
            f1out[i][0] = f1out[num - i][0] * f2out[num - i][0] - f1out[num - i][1] * f2out[num - i][1];
            f1out[i][1] = -f1out[num - i][0] * f2out[num - i][1] - f1out[num - i][1] * f2out[num - i][0];
        }
        f1out[i][0] /= num;
        f1out[i][1] /= num;
    }
    fftw_execute(pp);
}

void AfterConvolutionFFT()
{
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(pp);
    fftw_free(f1in);
    fftw_free(f1out);
    fftw_free(f2in);
    fftw_free(f2out);
}
