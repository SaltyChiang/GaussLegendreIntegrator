#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include "integral.hpp"
#include "convolution.hpp"
#include "gaussian.hpp"

int main()
{
	//	GaussLegendreIntegrator integrator(10000);
	// 	long start = clock();
	// #pragma omp parallel for
	// 	for (int i = -1000; i < 1000; i += 50)
	// 	{
	// 		// printf("%18.12e\n", integrator.DoIntegral(i, i + 50, Gauss));
	// 		integrator.DoIntegral(i, i + 50, Gaussian);
	// 	}
	// 	long end = clock();
	// 	printf("%ld\n", end - start);
	Convolution convolution(-10, 10, Gaussian, true, 1000);
	
	std::cout << convolution.DoConvolutionFFT(0.0) << std::endl;
	std::cout << convolution.DoConvolution(0.0) << std::endl;
	return 1;
}