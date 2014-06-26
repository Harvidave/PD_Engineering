#include "DataSEPD.h"
#include<math.h>
#include<iostream>

using namespace std;


DataSEPD::DataSEPD(std::vector<double> x, std::vector<double> y)
	:m_x(x),
	m_y(y)
{
	if (x.size() == y.size()){
		m_size = x.size();
	}
}


//---------------------------------------------------------------------------
bool
DataSEPD::precheck()
{
	bool check = true;

	if (m_size < 2)       // total size of data is less than 2, can not be solved
		check = false;

	double x, y;
	for (int i = 0; i < m_size && check; i++)
	{
		x = m_x[i];
		y = m_y[i];

		if (x < 0 || y < 0)
			check = false;
	}

	if (check == false)
	{
		std::cout << " Input value is not sufficient to do decline analysis \n";
	}

	return check;
}




//---------------------------------------------------------------------------
// q = q0 * exp(- (t/tao)^n ) 
std::vector<double>
DataSEPD::compute(std::vector<double> future)
{
	// simply check the data, (size >=2, all x > 0 , all y > 0)
	bool a = precheck();
	if (a == false)
		return std::vector<double>();

	SEPD_match();

	return SEPD_forecast(future);
}



//---------------------------------------------------------------------------
void
DataSEPD::SEPD_match()
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


	// SEPD model match, initial guess: tao = 1.0; n = 1.1; q0 = first data
	m_tao = 1.0;
	m_n = 1.1;
	m_q0 = m_y[0];

	// calculate Jocobian J, J1 = de/dqinf, J2 = de/da, J3 = de/dm
	double q0, tao, n;
	q0 = x[0] = m_q0;
	tao = x[1] = m_tao;
	n = x[2] = m_n;

	double R2, SSres = 0;
	for (int i = 0; i < m_size; i++)
	{
		double dt = m_x[i];

		double fi = m_q0 * exp(-pow(dt / tao, n));

		SSres += (fi - m_y[i]) * (fi - m_y[i]);
	}
	R2 = 1 - SSres / SStot;

	int i_iter = 0;

	while (R2 < 0.9999 && i_iter < 36)  // maximum 16 iterations
	{
		// step1: assign qi, ai, b, and calculate Jocabian
		q0 = m_q0;
		tao = m_tao;
		n = m_n;

		for (int i = 0; i < m_size; i++)
		{
			double dt = m_x[i];

			double temp = exp(-pow(dt / tao, n));

			J1[i] = temp;
			J2[i] = m_q0 * temp * pow(dt / tao, n) * n / tao;
			J3[i] = -m_q0 * temp * pow(dt / tao, n) * log(dt / tao);

			y[i] = m_q0 * temp;
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
			double dt = m_x[i];
			double temp = exp(-pow(dt / (x[1] + dx[1]), x[2] + dx[2]));
			if (temp < 1.0e-10)
				temp = 1.0e-10;

			double fi = (x[0] + dx[0]) * temp;

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
			double dt = m_x[i];
			double temp = exp(-pow(dt / (x[1] + dx1[1]), x[2] + dx1[2]));
			if (temp < 1.0e-10)
				temp = 1.0e-10;

			double fi = (x[0] + dx1[0]) *  temp;

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
			double dt = m_x[i];
			double temp = exp(-pow(dt / (x[1] + dx2[1]), x[2] + dx2[2]));
			if (temp < 1.0e-10)
				temp = 1.0e-10;

			double fi = (x[0] + dx2[0]) * temp;

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
				m_q0 = x[0] + dx[0];
				m_tao = x[1] + dx[1];
				m_n = x[2] + dx[2];
			}
			else
			{
				R2 = R2_2;
				lambda = lambda / v;
				m_q0 = x[0] + dx2[0];
				m_tao = x[1] + dx2[1];
				m_n = x[2] + dx2[2];
			}
		}
		else
		{
			if (R2_1 >= R2_2)
			{
				R2 = R2_1;
				lambda = lambda * v;
				m_q0 = x[0] + dx1[0];
				m_tao = x[1] + dx1[1];
				m_n = x[2] + dx1[2];
			}
			else
			{
				R2 = R2_2;
				lambda = lambda / v;
				m_q0 = x[0] + dx2[0];
				m_tao = x[1] + dx2[1];
				m_n = x[2] + dx2[2];
			}
		}

		x[0] = m_q0;
		x[1] = m_tao;
		x[2] = m_n;
		i_iter++;   // increase iteration
	}  // end of L-M loop

	//
}



//---------------------------------------------------------------------------
std::vector<double>
DataSEPD::SEPD_forecast(std::vector<double> future)
{
	std::vector<double> results;
	for (int i = 0; i < future.size(); i++) {
		results.push_back(m_q0 * exp(0.0 - pow(future[i] / m_tao, m_n)));
	}
	return results;
}


//---------------------------------------------------------------------------
// The following is gaussian elimination 
#define mat_elem(a, y, x, n) (a + ((y) * (n) + (x)))

void
DataSEPD::swap_row(double *a, double *b, int r1, int r2, int n)
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
DataSEPD::gauss_eliminate(double *a, double *b, double *x, int n)
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



