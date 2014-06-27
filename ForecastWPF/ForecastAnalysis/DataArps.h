#include <vector>


#ifndef DataArps_Class_
#define DataArps_Class_

class DataArps
{
public:
	DataArps(int i_method, std::vector<double> x, std::vector<double> y);

	std::vector<double>    computeForecast(std::vector<double> future);
	double		computeEur(double Qi, double qf);

private:
	int                                           m_method;
	std::vector<double>							m_x;
	std::vector<double>							m_y;
	int                                           m_size;

	double                                        m_A0;
	double                                        m_A1;

	double                                        m_qi;
	double                                        m_ai;
	double                                        m_b;

	bool    precheck();

	std::vector<double>    exponential_rate_time(std::vector<double> future);
	std::vector<double>    exponential_rate_Q(std::vector<double> future);

	std::vector<double>    harmonic_rate_time(std::vector<double> future);
	std::vector<double>    harmonic_rate_Q(std::vector<double> future);

	std::vector<double>    hyperbolic_onestep_rate_time(double b, bool forecast, std::vector<double> future);
	std::vector<double>    hyperbolic_onestep_rate_Q(double b, bool forecast, std::vector<double> future);

	std::vector<double>    hyperbolic_regression_rate_time(std::vector<double> future);
	std::vector<double>    hyperbolic_regression_rate_Q(std::vector<double> future);

	void    swap_row(double *a, double *b, int r1, int r2, int n);
	void    gauss_eliminate(double *a, double *b, double *x, int n);
};


#endif
