#pragma once

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