#include "DataArps.h"
#include<math.h>
#include<iostream>

using namespace std;

DataArps::DataArps(int i_method, std::vector<double> x, std::vector<double> y)
	:m_method(i_method),
	m_x(x),
	m_y(y)
{
	if (x.size() == y.size()){
		m_size = x.size();
	}
}


//---------------------------------------------------------------------------
bool
DataArps::precheck()
{
	if (m_size < 2)       // total size of data is less than 2, can not be solved
		return false;

	for (int i = 0; i < m_size; i++)
	{
		if (m_x[i] < 0 || m_y[i] < 0){
			return false;
		}
	}
	return true;
}

//---------------------------------------------------------------------------
// public compute forecast function
std::vector<double>
DataArps::computeForecast(std::vector<double> future)
{
	// simply check the data, (size >=2, all x > 0 , all y > 0)
	if (precheck() == false)
		return std::vector<double>();

	// exponential decline
	if (m_method == 0)
		return exponential_rate_time(future);
	else if (m_method == 1)
		return exponential_rate_Q(future);

	// harmonic decline
	else if (m_method == 2)
		return harmonic_rate_time(future);
	else if (m_method == 3)
		return harmonic_rate_Q(future);

	// hyperbolic decline, n is specified
	else if (m_method == 4)
		return hyperbolic_onestep_rate_time(0.2, true, future);
	else if (m_method == 5)
		return hyperbolic_onestep_rate_Q(0.2, true, future);

	//hyperbolic decline, automatically nonlinear-regression
	else if (m_method == 6)
		return hyperbolic_regression_rate_time(future);
	else if (m_method == 7)
		return hyperbolic_regression_rate_Q(future);

	else
	{
		std::cout << " Do not have this case \n";
		return std::vector<double>();
	}

}

//---------------------------------------------------------------------------
// public compute eur function
double
DataArps::computeEur(double Qi, double qf)
{
	// simply check the data, (size >=2, all x > 0 , all y > 0)
	if (precheck() == false)
		return -1.0;

	// exponential decline
	if (m_method == 0) {
		exponential_rate_time(std::vector<double>());
		return Qi + (m_qi - qf) / m_ai;
	}
	else if (m_method == 1) {
		exponential_rate_Q(std::vector<double>());
		return Qi + (m_qi - qf) / m_ai;
	}

	// harmonic decline
	else if (m_method == 2) {
		harmonic_rate_time(std::vector<double>());
		return Qi + m_qi / m_ai * log(m_ai / qf);
	}
	else if (m_method == 3) {
		harmonic_rate_Q(std::vector<double>());
		return Qi + m_qi / m_ai * log(m_ai / qf);
	}

	// hyperbolic decline, n is specified
	else if (m_method == 4) {
		hyperbolic_onestep_rate_time(0.2, true, std::vector<double>());
		return Qi + pow(m_qi, m_b) * (pow(m_qi, 1 - m_b) - pow(qf, 1 - m_b)) / m_ai / (1 - m_b);
	}
	else if (m_method == 5) {
		hyperbolic_onestep_rate_Q(0.2, true, std::vector<double>());
		return Qi + pow(m_qi, m_b) * (pow(m_qi, 1 - m_b) - pow(qf, 1 - m_b)) / m_ai / (1 - m_b);
	}

	//hyperbolic decline, automatically nonlinear-regression
	else if (m_method == 6) {
		hyperbolic_regression_rate_time(std::vector<double>());
		return Qi + pow(m_qi, m_b) * (pow(m_qi, 1 - m_b) - pow(qf, 1 - m_b)) / m_ai / (1 - m_b);
	}
	else if (m_method == 7) {
		hyperbolic_regression_rate_Q(std::vector<double>());
		return Qi + pow(m_qi, m_b) * (pow(m_qi, 1 - m_b) - pow(qf, 1 - m_b)) / m_ai / (1 - m_b);
	}

	else
	{
		std::cout << " Do not have this case \n";
		return -1.0;
	}

}



//---------------------------------------------------------------------------
// exponential decline, rate vs time; here ai = a = const
// q = qi * exp( - a * t )
// known: q, t 
// unknown: qi, a
// forecast: Qf = Qi + (qi - qf)/a
std::vector<double>
DataArps::exponential_rate_time(std::vector<double> future)
{
	// transform, y = ln(q), x = t; y = ln(qi) + (-a) * x

	// construct matrix
	//     [ a11,   a12 ]   [ ln(qi) ]   =  [c1]
	//     [ a21,   a22 ]   [ (-a)   ]      [c2]
	double a11, a12, a21, a22, c1, c2;
	a11 = a12 = a21 = a22 = c1 = c2 = 0;

	double x, y;
	for (int i = 0; i < m_size; i++)
	{
		x = m_x[i];
		y = log(m_y[i]);

		a11 += 1.0;
		a12 += x;
		a22 += x * x;

		c1 += y;
		c2 += x*y;
	}
	a21 = a12;  // symatrix 


	// solve A0, and A1
	m_A0 = (c1 * a22 - c2 * a12) / (a11 * a22 - a21 * a12);
	m_A1 = (c1 * a21 - c2 * a11) / (a12 * a21 - a11 * a22);

	// get qi, ai
	m_qi = exp(m_A0);
	m_ai = -m_A1;
	m_b = 0;

	// forecast
	std::vector<double> results;
	for (int i = 0; i < future.size(); i++) {
		results.push_back(m_qi * exp(0 - m_ai * future[i]));
	}
	return results;
}



//---------------------------------------------------------------------------
// exponential decline, rate vs Qumulative
// q = qi - a * Q
// known: q, Q 
// unknown: qi, a
// forecast: Qf = Qi + (qi - qf)/a
std::vector<double>
DataArps::exponential_rate_Q(std::vector<double> future)
{
	// construct matrix
	//     [ a11,   a12 ]   [ qi   ]  =  [c1]
	//     [ a21,   a22 ]   [ (-a) ]     [c2]
	double a11, a12, a21, a22, c1, c2;
	a11 = a12 = a21 = a22 = c1 = c2 = 0;

	double x, y;
	for (int i = 0; i < m_size; i++)
	{
		x = m_x[i];
		y = m_y[i];

		a11 += 1.0;
		a12 += x;
		a22 += x * x;

		c1 += y;
		c2 += x*y;
	}
	a21 = a12;  // symatrix 


	// solve A0, and A1
	m_A0 = (c1 * a22 - c2 * a12) / (a11 * a22 - a21 * a12);
	m_A1 = (c1 * a21 - c2 * a11) / (a12 * a21 - a11 * a22);

	// get qi, ai
	m_qi = m_A0;
	m_ai = -m_A1;
	m_b = 0;

	// forecast
	std::vector<double> results;
	for (int i = 0; i < future.size(); i++) {
		results.push_back(m_qi - m_ai * future[i]);
	}
	return results;
}



//---------------------------------------------------------------------------
// harmonic decline, rate vs time
// q = qi / (1 + ai * t)
// known: q, t 
// unknown: qi, ai
// forecast: Qf = Qi + (qi/ai) * Ln(qi/qf)
std::vector<double>
DataArps::harmonic_rate_time(std::vector<double> future)
{
	// transform, y = 1/q , x = t; y = 1/qi + (ai/qi) * x

	// construct matrix
	//     [ a11,   a12 ]   [ 1/qi  ]   =  [c1]
	//     [ a21,   a22 ]   [ ai/qi ]      [c2]
	double a11, a12, a21, a22, c1, c2;
	a11 = a12 = a21 = a22 = c1 = c2 = 0;

	double x, y;
	for (int i = 0; i < m_size; i++)
	{
		x = m_x[i];
		y = 1 / m_y[i];

		a11 += 1.0;
		a12 += x;
		a22 += x * x;

		c1 += y;
		c2 += x*y;
	}
	a21 = a12;  // symatrix 


	// solve A0, and A1
	m_A0 = (c1 * a22 - c2 * a12) / (a11 * a22 - a21 * a12);
	m_A1 = (c1 * a21 - c2 * a11) / (a12 * a21 - a11 * a22);

	// get qi, ai
	m_qi = 1 / m_A0;
	m_ai = m_A1 * m_qi;
	m_b = 1;

	// forecast
	std::vector<double> results;
	for (int i = 0; i < future.size(); i++) {
		results.push_back(m_qi / (1 + m_ai * future[i]));
	}
	return results;
}



//---------------------------------------------------------------------------
// harmonic decline, rate vs Qumulative
// q = qi * exp( -ai/qi * Q)
// known: q, Q 
// unknown: qi, ai
// forecast: Qf = Qi + (qi/ai) * Ln(qi/qf)
std::vector<double>
DataArps::harmonic_rate_Q(std::vector<double> future)
{
	// transform, y = ln(q) , x = Q; y = ln(qi) + (-ai/qi) * x

	// construct matrix
	//     [ a11,   a12 ]   [ ln(qi)  ]   =  [c1]
	//     [ a21,   a22 ]   [ -ai/qi  ]      [c2]
	double a11, a12, a21, a22, c1, c2;
	a11 = a12 = a21 = a22 = c1 = c2 = 0;

	double x, y;
	for (int i = 0; i < m_size; i++)
	{
		x = m_x[i];
		y = log(m_y[i]);

		a11 += 1.0;
		a12 += x;
		a22 += x * x;

		c1 += y;
		c2 += x*y;
	}
	a21 = a12;  // symatrix 


	// solve A0, and A1
	m_A0 = (c1 * a22 - c2 * a12) / (a11 * a22 - a21 * a12);
	m_A1 = (c1 * a21 - c2 * a11) / (a12 * a21 - a11 * a22);

	// get qi, ai
	m_qi = exp(m_A0);
	m_ai = -m_A1 * m_qi;
	m_b = 1;

	// forecast
	std::vector<double> results;
	for (int i = 0; i < future.size(); i++) {
		results.push_back(m_qi * exp(0 - m_ai / m_qi * future[i]));
	}
	return results;
}



//---------------------------------------------------------------------------
// hyperbolic decline, rate vs time, b is given, attention: b != 0
// q = qi / (1 + b * ai * t)^(1/b)
// known: q, t 
// unknown: qi, ai 
// forecast: Qf = Qi + (qi^b/ai/(1-b)) * (qi^(1-b) - qf^(1-b))
std::vector<double>
DataArps::hyperbolic_onestep_rate_time(double b, bool forecast, std::vector<double> future)
{
	if (b == 0)
		return std::vector<double>();
	// transform, y = (1/q)^b , x = t; y = 1/qi^b + b * ai/qi^b * x

	// construct matrix
	//     [ a11,   a12 ]   [ 1/qi^b   ]   =  [c1]
	//     [ a21,   a22 ]   [ ai/qi^b  ]      [c2]
	double a11, a12, a21, a22, c1, c2;
	a11 = a12 = a21 = a22 = c1 = c2 = 0;

	double x, y;
	for (int i = 0; i < m_size; i++)
	{
		x = m_x[i];
		y = pow(1 / m_y[i], b);

		a11 += 1.0;
		a12 += x;
		a22 += x * x;

		c1 += y;
		c2 += x*y;
	}
	a21 = a12;  // symatrix 


	// solve A0, and A1
	m_A0 = (c1 * a22 - c2 * a12) / (a11 * a22 - a21 * a12);
	m_A1 = (c1 * a21 - c2 * a11) / (a12 * a21 - a11 * a22);

	// get qi, ai
	if (m_A0 > 0)
	{
		m_qi = pow(1 / m_A0, 1 / b);
		m_ai = m_A1 / m_A0 / b;
	}
	else
	{
		if (forecast)
			std::cout << " the input exponential value, " << b << ", leads negative value for power function, no match is found \n";
		else
		{
			// just give an initial guess
			m_qi = pow(abs(1 / m_A0), 1 / b);
			m_ai = abs(m_A1 / m_A0) / b;
		}
	}

	// forecast
	if (forecast)
	{
		std::vector<double> results;
		for (int i = 0; i < future.size(); i++) {
			results.push_back(m_qi / pow(1 + b * m_ai * future[i], 1 / b));
		}
		return results;
	}

	return std::vector<double>();
}



//---------------------------------------------------------------------------
// hyperbolic decline, rate vs Qumulative, b is given, b != 1
// q^(1-b) = qi^(1-b) - Q * ai*(1-b)/qi^b
// known: q, Q 
// unknown: qi, ai 
// forecast: Qf = Qi + (qi^b/ai/(1-b)) * (qi^(1-b) - qf^(1-b))
std::vector<double>
DataArps::hyperbolic_onestep_rate_Q(double b, bool forecast, std::vector<double> future)
{
	if (b == 1)
		return std::vector<double>();
	// transform, y = q^(1-b) , x = Q; y = qi^(1-b) - ai*(1-b)/qi^b * x

	// construct matrix
	//     [ a11,   a12 ]   [ qi^(1-b)       ]   =  [c1]
	//     [ a21,   a22 ]   [ ai*(1-b)/qi^b  ]      [c2]
	double a11, a12, a21, a22, c1, c2;
	a11 = a12 = a21 = a22 = c1 = c2 = 0;

	double x, y;
	for (int i = 0; i < m_size; i++)
	{
		x = m_x[i];
		y = pow(m_y[i], 1 - b);

		a11 += 1.0;
		a12 += x;
		a22 += x * x;

		c1 += y;
		c2 += x*y;
	}
	a21 = a12;  // symatrix 


	// solve A0, and A1
	m_A0 = (c1 * a22 - c2 * a12) / (a11 * a22 - a21 * a12);
	m_A1 = (c1 * a21 - c2 * a11) / (a12 * a21 - a11 * a22);

	// get qi, ai
	if (m_A0 > 0)
	{
		m_qi = pow(m_A0, 1 / (1 - b));
		m_ai = -m_A1 *  pow(m_qi, b) / (1 - b);
	}
	else
	{
		if (forecast)
			std::cout << " the input exponential value, " << b << ", leads negative value for power function, no match is found \n";
		else
		{
			// just give an initial guess
			m_qi = pow(abs(m_A0), 1 / (1 - b));
			m_ai = abs(m_A1 *  pow(m_qi, b) / (1 - b));
		}
	}

	// forecast
	if (forecast)
	{
		std::vector<double> results;
		for (int i = 0; i < future.size(); i++) {
			results.push_back(pow(pow(m_qi, 1 - b) - future[i] * m_ai * (1 - b) / pow(m_qi, b), 1 / (1 - b)));
		}
		return results;
	}

	return std::vector<double>();
}



//---------------------------------------------------------------------------
// hyperbolic decline, rate vs time, nonlinear regression: LM method
// q = qi / (1 + b * ai * dt)^(1/b) 
// forecast: Qf = Qi + (qi^b/ai/(1-b)) * (qi^(1-b) - qf^(1-b))
std::vector<double>
DataArps::hyperbolic_regression_rate_time(std::vector<double> future)
{
	// Create space for Jocobian (column), Hessian matrix, and rhs
	vector<double> J1, J2, J3, y;
	J1.resize(m_size);
	J2.resize(m_size);
	J3.resize(m_size);
	y.resize(m_size);

	double Hessian[9], rhs[3], A[9], dx[3], x[3];
	double lambda0, lambda, v = 10;


	double yave_obs = 0;
	for (int i = 0; i < m_size; i++)
		yave_obs += m_y[i];

	yave_obs /= m_size;

	double SStot = 0;
	for (int i = 0; i < m_size; i++)
		SStot += (m_y[i] - yave_obs) * (m_y[i] - yave_obs);


	// Initial guess, b = 0.5
	m_b = 0.5;
	hyperbolic_onestep_rate_time(m_b, false, future);  // after this, get initial guess of m_qi and m_ai 

	// calculate Jocobian J, J1 = de/dqi, J2 = de/dai, J3 = de/db
	double qi, ai, b;
	qi = x[0] = m_qi;
	ai = x[1] = m_ai;
	b = x[2] = m_b;

	double R2, SSres = 0;
	for (int i = 0; i < m_size; i++)
	{
		double fi = qi / pow(1 + b * ai * m_x[i], 1 / b);
		SSres += (fi - m_y[i]) * (fi - m_y[i]);
	}
	R2 = 1 - SSres / SStot;

	int i_iter = 0;

	while (R2 < 0.9999998 && i_iter < 16)  // maximum 16 iterations
	{
		// step1: assign qi, ai, b, and calculate Jocabian
		qi = m_qi;
		ai = m_ai;
		b = m_b;

		for (int i = 0; i < m_size; i++)
		{
			double dt = m_x[i];

			double temp = 1 + b * ai * dt;

			J1[i] = 1 / pow(temp, 1 / b);
			J2[i] = qi / pow(temp, 1 / b + 1) * (-dt);
			J3[i] = qi / pow(temp, 1 / b) * (1 / b * log(temp) - ai * dt / temp) / b;

			y[i] = qi / pow(temp, 1 / b);
		}

		// step2: calculate Hessian and rhs
		for (int i = 0; i < 9; i++)
			Hessian[i] = 0;

		rhs[0] = rhs[1] = rhs[2] = 0;

		for (int i = 0; i < m_size; i++)
		{
			Hessian[0] += J1[i] * J1[i];
			Hessian[1] += J1[i] * J2[i];
			Hessian[2] += J1[i] * J3[i];

			Hessian[3] += J2[i] * J1[i];
			Hessian[4] += J2[i] * J2[i];
			Hessian[5] += J2[i] * J3[i];

			Hessian[6] += J3[i] * J1[i];
			Hessian[7] += J3[i] * J2[i];
			Hessian[8] += J3[i] * J3[i];

			rhs[0] += J1[i] * (m_y[i] - y[i]);
			rhs[1] += J2[i] * (m_y[i] - y[i]);
			rhs[2] += J3[i] * (m_y[i] - y[i]);
		}

		// select lambda = average of diag(Hessian) , and v = 10
		lambda0 = (Hessian[0] + Hessian[4] + Hessian[8]) / 3;

		// first iteration, guess a lambda
		if (i_iter == 0)
			lambda = lambda0;

		A[0] = Hessian[0] + lambda * Hessian[0];
		A[1] = Hessian[1];
		A[2] = Hessian[2];
		A[3] = Hessian[3];
		A[4] = Hessian[4] + lambda * Hessian[4];
		A[5] = Hessian[5];
		A[6] = Hessian[6];
		A[7] = Hessian[7];
		A[8] = Hessian[8] + lambda * Hessian[8];

		// solve delta_paramter, Gaussian elimination
		double rhs0[3];
		rhs0[0] = rhs[0];
		rhs0[1] = rhs[1];
		rhs0[2] = rhs[2];

		gauss_eliminate(A, rhs0, dx, 3);   // After this one, A and rhs0 are changed

		double R2_0, SSres0 = 0;
		for (int i = 0; i < m_size; i++)
		{
			double temp = pow(1 + (x[2] + dx[2]) * (x[1] + dx[1]) * m_x[i], 1 / (x[2] + dx[2]));
			if (temp < 1.0e-10)
				temp = 1.0e-10;

			double fi = (x[0] + dx[0]) / temp;
			SSres0 += (fi - m_y[i]) * (fi - m_y[i]);
		}
		R2_0 = 1 - SSres0 / SStot;

		// check lambda*v for best update
		A[0] = Hessian[0] + lambda * v *  Hessian[0];
		A[1] = Hessian[1];
		A[2] = Hessian[2];
		A[3] = Hessian[3];
		A[4] = Hessian[4] + lambda * v * Hessian[4];
		A[5] = Hessian[5];
		A[6] = Hessian[6];
		A[7] = Hessian[7];
		A[8] = Hessian[8] + lambda * v * Hessian[8];

		double rhs1[3];
		rhs1[0] = rhs[0];
		rhs1[1] = rhs[1];
		rhs1[2] = rhs[2];

		double dx1[3];
		gauss_eliminate(A, rhs1, dx1, 3);

		double R2_1, SSres1 = 0;
		for (int i = 0; i < m_size; i++)
		{
			double temp = pow(1 + (x[2] + dx1[2]) * (x[1] + dx1[1]) * m_x[i], 1 / (x[2] + dx1[2]));
			if (temp < 1.0e-10)
				temp = 1.0e-10;

			double fi = (x[0] + dx1[0]) / temp;
			SSres1 += (fi - m_y[i]) * (fi - m_y[i]);
		}
		R2_1 = 1 - SSres1 / SStot;

		// check lambda/v for best update
		A[0] = Hessian[0] + lambda / v *  Hessian[0];
		A[1] = Hessian[1];
		A[2] = Hessian[2];
		A[3] = Hessian[3];
		A[4] = Hessian[4] + lambda / v * Hessian[4];
		A[5] = Hessian[5];
		A[6] = Hessian[6];
		A[7] = Hessian[7];
		A[8] = Hessian[8] + lambda / v * Hessian[8];

		double rhs2[3], dx2[3];
		rhs2[0] = rhs[0];
		rhs2[1] = rhs[1];
		rhs2[2] = rhs[2];

		gauss_eliminate(A, rhs2, dx2, 3);

		double R2_2, SSres2 = 0;
		for (int i = 0; i < m_size; i++)
		{
			double temp = pow(1 + (x[2] + dx2[2]) * (x[1] + dx2[1]) * m_x[i], 1 / (x[2] + dx2[2]));
			if (temp < 1.0e-10)
				temp = 1.0e-10;

			double fi = (x[0] + dx2[0]) / temp;
			SSres2 += (fi - m_y[i]) * (fi - m_y[i]);
		}
		R2_2 = 1 - SSres2 / SStot;


		// Choose the best value of lambda for next iteration, and update result
		if (R2_0 >= R2_1)
		{
			if (R2_0 >= R2_2)
			{
				R2 = R2_0;
				lambda = lambda;
				m_qi = x[0] + dx[0];
				m_ai = x[1] + dx[1];
				m_b = x[2] + dx[2];
			}
			else
			{
				R2 = R2_2;
				lambda = lambda / v;
				m_qi = x[0] + dx2[0];
				m_ai = x[1] + dx2[1];
				m_b = x[2] + dx2[2];
			}
		}
		else
		{
			if (R2_1 >= R2_2)
			{
				R2 = R2_1;
				lambda = lambda * v;
				m_qi = x[0] + dx1[0];
				m_ai = x[1] + dx1[1];
				m_b = x[2] + dx1[2];
			}
			else
			{
				R2 = R2_2;
				lambda = lambda / v;
				m_qi = x[0] + dx2[0];
				m_ai = x[1] + dx2[1];
				m_b = x[2] + dx2[2];
			}
		}

		x[0] = m_qi;
		x[1] = m_ai;
		x[2] = m_b;
		i_iter++;   // increase iteration
	}  // end of L-M loop


	//forecast
	std::vector<double> results;
	for (int i = 0; i < future.size(); i++) {
		results.push_back(m_qi / pow(1 + m_b * m_ai * future[i], 1 / m_b));
	}
	return results;
}



//--------------------------------------------------------------------------
// hyperbolic decline, rate vs Qumulative, nonlinear regression: LM method
// q^(1-b) = qi^(1-b) - Q * ai * (1-b) / qi^b
// forecast: Qf = Qi + (qi^b/ai/(1-b)) * (qi^(1-b) - qf^(1-b))
std::vector<double>
DataArps::hyperbolic_regression_rate_Q(std::vector<double> future)
{
	// Create space for Jocobian (column), Hessian matrix, and rhs
	vector<double> J1, J2, J3, y;
	J1.resize(m_size);
	J2.resize(m_size);
	J3.resize(m_size);
	y.resize(m_size);

	double Hessian[9], rhs[3], A[9], dx[3], x[3];
	double lambda0, lambda, v = 10;


	double yave_obs = 0;
	for (int i = 0; i < m_size; i++)
		yave_obs += m_y[i];

	yave_obs /= m_size;

	double SStot = 0;
	for (int i = 0; i < m_size; i++)
		SStot += (m_y[i] - yave_obs) * (m_y[i] - yave_obs);


	// Initial guess, b = 0.5
	m_b = 0.5;
	hyperbolic_onestep_rate_Q(m_b, false, future);  // after this, get initial guess of m_qi and m_ai 

	// calculate Jocobian J, J1 = de/dqi, J2 = de/dai, J3 = de/db
	double qi, ai, b;
	qi = x[0] = m_qi;
	ai = x[1] = m_ai;
	b = x[2] = m_b;

	double R2, SSres = 0;
	for (int i = 0; i < m_size; i++)
	{
		double temp, fi;
		temp = pow(qi, (1 - b)) - m_x[i] * ai * (1 - b) / pow(qi, b);
		fi = pow(temp, 1 / (1 - b));

		SSres += (fi - m_y[i]) * (fi - m_y[i]);
	}
	R2 = 1 - SSres / SStot;

	int i_iter = 0;

	while (R2 < 0.99998 && i_iter < 16)  // maximum 16 iterations
	{
		// step1: assign qi, ai, b, and calculate Jocabian
		qi = m_qi;
		ai = m_ai;
		b = m_b;

		for (int i = 0; i < m_size; i++)
		{
			double Q = m_x[i];

			double temp, y1, y2;
			temp = pow(qi, (1 - b)) - Q * ai * (1 - b) / pow(qi, b);
			y1 = pow(temp, 1 / (1 - b));

			temp = pow(qi*(1 + 0.0001), (1 - b)) - Q * ai * (1 - b) / pow(qi*(1 + 0.0001), b);
			y2 = pow(temp, 1 / (1 - b));

			J1[i] = (y2 - y1) / (0.0001 * qi);

			// ----------------
			temp = pow(qi, (1 - b)) - Q * ai *(1 + 0.0001) * (1 - b) / pow(qi, b);
			y2 = pow(temp, 1 / (1 - b));

			J2[i] = (y2 - y1) / (0.0001 * ai);

			// ----------------
			temp = pow(qi, (1 - b*(1 + 0.0001))) - Q * ai * (1 - b*(1 + 0.0001)) / pow(qi, b*(1 + 0.0001));
			y2 = pow(temp, 1 / (1 - b*(1 + 0.0001)));

			J3[i] = (y2 - y1) / (0.0001 * b);

			y[i] = y1;
		}

		// step2: calculate Hessian and rhs
		for (int i = 0; i < 9; i++)
			Hessian[i] = 0;

		rhs[0] = rhs[1] = rhs[2] = 0;

		for (int i = 0; i < m_size; i++)
		{
			Hessian[0] += J1[i] * J1[i];
			Hessian[1] += J1[i] * J2[i];
			Hessian[2] += J1[i] * J3[i];

			Hessian[3] += J2[i] * J1[i];
			Hessian[4] += J2[i] * J2[i];
			Hessian[5] += J2[i] * J3[i];

			Hessian[6] += J3[i] * J1[i];
			Hessian[7] += J3[i] * J2[i];
			Hessian[8] += J3[i] * J3[i];

			rhs[0] += J1[i] * (m_y[i] - y[i]);
			rhs[1] += J2[i] * (m_y[i] - y[i]);
			rhs[2] += J3[i] * (m_y[i] - y[i]);
		}

		// select lambda = average of diag(Hessian) , and v = 10
		lambda0 = (Hessian[0] + Hessian[4] + Hessian[8]) / 3;

		// first iteration, guess a lambda
		if (i_iter == 0)
			lambda = lambda0;

		A[0] = Hessian[0] + lambda * Hessian[0];
		A[1] = Hessian[1];
		A[2] = Hessian[2];
		A[3] = Hessian[3];
		A[4] = Hessian[4] + lambda * Hessian[4];
		A[5] = Hessian[5];
		A[6] = Hessian[6];
		A[7] = Hessian[7];
		A[8] = Hessian[8] + lambda * Hessian[8];

		// solve delta_paramter, Gaussian elimination
		double rhs0[3];
		rhs0[0] = rhs[0];
		rhs0[1] = rhs[1];
		rhs0[2] = rhs[2];

		gauss_eliminate(A, rhs0, dx, 3);   // After this one, A and rhs0 are changed

		double R2_0, SSres0 = 0;
		for (int i = 0; i < m_size; i++)
		{
			qi = x[0] + dx[0];
			ai = x[1] + dx[1];
			b = x[2] + dx[2];
			double temp = pow(qi, (1 - b)) - m_x[i] * ai * (1 - b) / pow(qi, b);
			if (temp < 1.0e-10)
				temp = 1.0e-10;

			double fi = pow(temp, 1 / (1 - b));
			SSres0 += (fi - m_y[i]) * (fi - m_y[i]);
		}
		R2_0 = 1 - SSres0 / SStot;

		// check lambda*v for best update
		A[0] = Hessian[0] + lambda * v *  Hessian[0];
		A[1] = Hessian[1];
		A[2] = Hessian[2];
		A[3] = Hessian[3];
		A[4] = Hessian[4] + lambda * v * Hessian[4];
		A[5] = Hessian[5];
		A[6] = Hessian[6];
		A[7] = Hessian[7];
		A[8] = Hessian[8] + lambda * v * Hessian[8];

		double rhs1[3];
		rhs1[0] = rhs[0];
		rhs1[1] = rhs[1];
		rhs1[2] = rhs[2];

		double dx1[3];
		gauss_eliminate(A, rhs1, dx1, 3);

		double R2_1, SSres1 = 0;
		for (int i = 0; i < m_size; i++)
		{
			qi = x[0] + dx1[0];
			ai = x[1] + dx1[1];
			b = x[2] + dx1[2];
			double temp = pow(qi, (1 - b)) - m_x[i] * ai * (1 - b) / pow(qi, b);
			if (temp < 1.0e-10)
				temp = 1.0e-10;

			double fi = pow(temp, 1 / (1 - b));
			SSres1 += (fi - m_y[i]) * (fi - m_y[i]);
		}
		R2_1 = 1 - SSres1 / SStot;

		// check lambda/v for best update
		A[0] = Hessian[0] + lambda / v *  Hessian[0];
		A[1] = Hessian[1];
		A[2] = Hessian[2];
		A[3] = Hessian[3];
		A[4] = Hessian[4] + lambda / v * Hessian[4];
		A[5] = Hessian[5];
		A[6] = Hessian[6];
		A[7] = Hessian[7];
		A[8] = Hessian[8] + lambda / v * Hessian[8];

		double rhs2[3], dx2[3];
		rhs2[0] = rhs[0];
		rhs2[1] = rhs[1];
		rhs2[2] = rhs[2];

		gauss_eliminate(A, rhs2, dx2, 3);

		double R2_2, SSres2 = 0;
		for (int i = 0; i < m_size; i++)
		{
			qi = x[0] + dx2[0];
			ai = x[1] + dx2[1];
			b = x[2] + dx2[2];
			double temp = pow(qi, (1 - b)) - m_x[i] * ai * (1 - b) / pow(qi, b);
			if (temp < 1.0e-10)
				temp = 1.0e-10;

			double fi = pow(temp, 1 / (1 - b));
			SSres2 += (fi - m_y[i]) * (fi - m_y[i]);
		}
		R2_2 = 1 - SSres2 / SStot;


		// Choose the best value of lambda for next iteration, and update result
		if (R2_0 >= R2_1)
		{
			if (R2_0 >= R2_2)
			{
				R2 = R2_0;
				lambda = lambda;
				m_qi = x[0] + dx[0];
				m_ai = x[1] + dx[1];
				m_b = x[2] + dx[2];
			}
			else
			{
				R2 = R2_2;
				lambda = lambda / v;
				m_qi = x[0] + dx2[0];
				m_ai = x[1] + dx2[1];
				m_b = x[2] + dx2[2];
			}
		}
		else
		{
			if (R2_1 >= R2_2)
			{
				R2 = R2_1;
				lambda = lambda * v;
				m_qi = x[0] + dx1[0];
				m_ai = x[1] + dx1[1];
				m_b = x[2] + dx1[2];
			}
			else
			{
				R2 = R2_2;
				lambda = lambda / v;
				m_qi = x[0] + dx2[0];
				m_ai = x[1] + dx2[1];
				m_b = x[2] + dx2[2];
			}
		}

		x[0] = m_qi;
		x[1] = m_ai;
		x[2] = m_b;
		i_iter++;   // increase iteration
	}  // end of L-M loop


	//forecast
	std::vector<double> results;
	for (int i = 0; i < future.size(); i++) {
		results.push_back(pow(pow(m_qi, 1 - m_b) - future[i] * m_ai * (1 - b) / pow(m_qi, b), 1 / (1 - b)));
	}
	return results;
}


//---------------------------------------------------------------------------
// The following is gaussian elimination 
#define mat_elem(a, y, x, n) (a + ((y) * (n) + (x)))

void
DataArps::swap_row(double *a, double *b, int r1, int r2, int n)
{
	double tmp, *p1, *p2;
	int i;

	if (r1 == r2) return;
	for (i = 0; i < n; i++) {
		p1 = mat_elem(a, r1, i, n);
		p2 = mat_elem(a, r2, i, n);
		tmp = *p1, *p1 = *p2, *p2 = tmp;
	}
	tmp = b[r1], b[r1] = b[r2], b[r2] = tmp;
}

void
DataArps::gauss_eliminate(double *a, double *b, double *x, int n)
{
#define A(y, x) (*mat_elem(a, y, x, n))
	int i, j, col, row, max_row;
	double max, tmp;

	for (row = 0; row < n; row++) {
		max_row = row, max = A(row, row);

		for (i = row + 1; i < n; i++)
			if ((tmp = fabs(A(i, row))) > max)
				max_row = i, max = tmp;

		swap_row(a, b, row, max_row, n);

		for (i = row + 1; i < n; i++) {
			tmp = A(i, row) / A(row, row);
			for (col = row + 1; col < n; col++)
				A(i, col) -= tmp * A(row, col);
			A(i, row) = 0;
			b[i] -= tmp * b[row];
		}
	}
	for (row = n - 1; row >= 0; row--) {
		tmp = b[row];
		for (j = n - 1; j > row; j--)
			tmp -= x[j] * A(row, j);
		x[row] = tmp / A(row, row);
	}
#undef A
}

