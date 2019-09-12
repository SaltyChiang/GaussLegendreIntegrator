#include <cmath>
#include <fftw3.h>
#include "convolution.hpp"
#include "integral.hpp"
#include "gaussian.hpp"

namespace Convolution
{
    double Convolution(double a, double b, double x, double function1(double), bool useFFT)
    {
        Convolution::useFFT = useFFT;
        Convolution::a = a;
        Convolution::b = b;
        Convolution::x = x;
        Convolution::function1 = function1;
        Convolution::function2 = Gaussian;
    }

    double Convolution(double a, double b, double x, double function1(double), double function2(double), bool useFFT)
    {
        Convolution(a, b, x, function1, useFFT);
        Convolution::function2 = function2;
    }

    double PreConvolution(double y)
    {
        return function1(y) * function2(x - y);
    }

    double DoConvolution()
    {
        GaussLegendreIntegrator integrator(1000);
        return integrator.DoIntegral(a, b, PreConvolution);
    }

}
