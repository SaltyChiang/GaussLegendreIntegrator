#ifndef INTEGRAL_HPP_
#define INTEGRAL_HPP_

class GaussLegendreIntegrator
{
private:
	double* sampleAbscissa;
	double* sampleWeight;
	double eps;
	int sampleNum;
	void CalcGaussLegendreSamplingPoints();

public :
	GaussLegendreIntegrator(int sampleNum);
	double DoIntegral(double a, double b, double function(double));
};

#endif