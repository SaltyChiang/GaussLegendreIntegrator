#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

class GaussLegendreIntegrator
{
private:
	double* sampleAbscissa;
	double* sampleWeight;
	double eps;
	int sampleNum;

public :
	GaussLegendreIntegrator(int sampleNum)
	{
		this->sampleNum = sampleNum;
		this->eps = 1e-12;
		sampleAbscissa = new double[sampleNum];
		sampleWeight = new double[sampleNum];
	}

	double DoIntegral(double a, double b, double function(double x))
	{
		// Gauss-Legendre integral, see CalcGaussLegendreSamplingPoints.

		if (sampleNum <= 0 || sampleAbscissa == 0 || sampleWeight == 0)
			return 0;

		const double a0 = (b + a) / 2;
		const double b0 = (b - a) / 2;

		double xx;

		double result = 0.0;
		for (int i = 0; i < sampleNum; i++)
		{
			xx = a0 + b0 * sampleAbscissa[i];
			result += sampleWeight[i] * function(xx);
		}

		return result * b0;
	}

	void CalcGaussLegendreSamplingPoints()
	{
		// The roots of symmetric is the interval, so we only have to find half of them
		const unsigned int m = (sampleNum + 1) / 2;

		double z, pp, p1, p2, p3;

		// Loop over the desired roots
		for (unsigned int i = 0; i < m; i++) {
			z = std::cos(3.14159265358979323846 * (i + 0.75) / (sampleNum + 0.5));

			// Starting with the above approximation to the i-th root, we enter
			// the main loop of refinement by Newton's method
			do {
				p1 = 1.0;
				p2 = 0.0;

				// Loop up the recurrence relation to get the Legendre
				// polynomial evaluated at z
				for (int j = 0; j < sampleNum; j++)
				{
					p3 = p2;
					p2 = p1;
					p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1.0);
				}
				// p1 is now the desired Legendre polynomial. We next compute pp, its
				// derivative, by a standard relation involving also p2, the polynomial
				// of one lower order
				pp = sampleNum * (z * p1 - p2) / (z * z - 1.0);
				// Newton's method
				z -= p1 / pp;

			} while (std::fabs(p1 / pp) > eps);

			// Put root and its symmetric counterpart
			sampleAbscissa[i] = -z;
			sampleAbscissa[sampleNum - i - 1] = z;

			// Compute the weight and put its symmetric counterpart
			sampleWeight[i] = 2.0 / ((1.0 - z * z) * pp * pp);
			sampleWeight[sampleNum - i - 1] = sampleWeight[i];
		}
	}
}; // namespace Integral

double Gauss(double x) 
{
	return std::exp(-x * x / 2);
}

int main()
{
	GaussLegendreIntegrator integrator(10000);
	integrator.CalcGaussLegendreSamplingPoints();
	long start = clock();
#pragma omp parallel for
	for (int i = -1000; i < 1000; i += 50) 
	{
		// printf("%18.12e\n", integrator.DoIntegral(i, i + 50, Gauss));
		integrator.DoIntegral(i, i + 50, Gauss);
	}
	long end = clock();
	printf("%ld\n", end - start);
    return 1;
}