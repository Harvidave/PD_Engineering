#include <vector>


#ifndef DataDuong_Class_
#define DataDuong_Class_


// Note: Duong only provide Rate vs. time function
//       Duong's time is daily or monthly, 1, 2, .... 10, .... 

class DataDuong
{
public:
	DataDuong(int i_method, std::vector<double> x, std::vector<double> y);
	std::vector<double>    computeForecast(std::vector<double> future);
	double		computeEur();

private:
	int                                           m_method;
	std::vector<double>							m_x;
	std::vector<double>							m_y;
	int                                           m_size;

	double                                        m_a;
	double                                        m_m;
	double                                        m_q1;
	double                                        m_qinf;


	bool    precheck();

	void    choose_q1(int iuse);

	void    Duong_match();


	void    swap_row(double *a, double *b, int r1, int r2, int n);
	void    gauss_eliminate(double *a, double *b, double *x, int n);

};


#endif
