#include <vector>


#ifndef DataSEPD_Class_
#define DataSEPD_Class_

// SEPD method: i_method/m_method is not used, but leave this parameter here 
class DataSEPD
{
public:
	DataSEPD(std::vector<double> x, std::vector<double> y);
	std::vector<double>    computeForecast(std::vector<double> future);
	double	ComputeEur();

private:
	std::vector<double>							m_x;
	std::vector<double>							m_y;
	int                                           m_size;


	double                                        m_q0;
	double                                        m_tao;
	double                                        m_n;


	bool    precheck();

	void    SEPD_match();

	void    swap_row(double *a, double *b, int r1, int r2, int n);
	void    gauss_eliminate(double *a, double *b, double *x, int n);

};


#endif
