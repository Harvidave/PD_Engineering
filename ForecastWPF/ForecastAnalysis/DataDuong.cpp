#include "DataDuong.h"
#include<math.h>
#include<iostream>

using namespace std;

DataDuong::DataDuong(int i_method, std::vector<double> x, std::vector<double> y)
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
DataDuong::precheck()
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
void
DataDuong::choose_q1(int iuse)
{
	// no data is provided
	if (m_size == 0)
		return;

	if (iuse == 0)          // use day 1's q as q1; assume the first data is day 1??
	{
		m_q1 = m_y[0];
	}
	else if (iuse == 1)    // use qmax as q1
	{
		double qmax;
		qmax = m_y[0];

		for (int i = 1; i < m_size; i++)
		{
			if (qmax < m_y[i])
			{
				qmax = m_y[i];
			}
		}

		m_q1 = qmax;
	}

}


//---------------------------------------------------------------------------
void
DataDuong::compute()
{
	// simply check the data, (size >=2, all x > 0 , all y > 0)
	bool a = precheck();
	if (a == false)
		return;


	// use m_method to control choosing q1
	choose_q1(m_method);

	Duong_match();

	Duong_forecast();

}



//---------------------------------------------------------------------------
void
DataDuong::Duong_match()
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


	// Duong model match, initial guess: a = 1.0; m = 1.1; qinf = 0. Note: m = 1.0 can not be found.
	//m_q1   = 40;
	m_a = 1.0;
	m_m = 1.1;
	m_qinf = 0;

	// calculate Jocobian J, J1 = de/dqinf, J2 = de/da, J3 = de/dm
	double qinf, a, m;
	qinf = x[0] = m_qinf;
	a = x[1] = m_a;
	m = x[2] = m_m;

	double R2, SSres = 0;
	for (int i = 0; i < m_size; i++)
	{
		double dt = m_x[i];
		double temp = (pow(dt, 1 - m) - 1) / (1 - m);
		double fi = m_q1 * pow(dt, -m) * exp(a * temp) + qinf;

		SSres += (fi - m_y[i]) * (fi - m_y[i]);
	}
	R2 = 1 - SSres / SStot;

	int i_iter = 0;

	while (R2 < 0.9999 && i_iter < 16)  // maximum 16 iterations
	{
		// step1: assign qi, ai, b, and calculate Jocabian
		qinf = m_qinf;
		a = m_a;
		m = m_m;

		for (int i = 0; i < m_size; i++)
		{
			double dt = m_x[i];

			double temp = (pow(dt, 1 - m) - 1) / (1 - m);

			J1[i] = 1;
			J2[i] = m_q1 * pow(dt, -m) * exp(a * temp) * temp;
			//J3[i] = m_q1 * (-m) * pow(dt, -m-1) * exp(a * temp) + 
			//        m_q1 * pow(dt, -m) * exp(a * temp) * a * (-1 + pow(dt, (1-m)) - (1-m)*(1-m)*pow(dt,-m) )/ ((1-m)*(1-m))  ;

			y[i] = m_q1 * pow(dt, -m) * exp(a * temp) + qinf;


			double m2 = m + 0.0001;
			double temp2 = (pow(dt, 1 - m2) - 1) / (1 - m2);
			double y2 = m_q1*pow(dt, -m2) * exp(a * temp2) + qinf;

			double dydm = (y2 - y[i]) / 0.0001;
			J3[i] = dydm;
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
			double temp = (pow(dt, 1 - (x[2] + dx[2])) - 1) / (1 - (x[2] + dx[2]));
			if (temp < 1.0e-10)
				temp = 1.0e-10;

			double fi = m_q1 * pow(dt, -(x[2] + dx[2])) * exp((x[1] + dx[1]) * temp) + (x[0] + dx[0]);

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
			double temp = (pow(dt, 1 - (x[2] + dx1[2])) - 1) / (1 - (x[2] + dx1[2]));
			if (temp < 1.0e-10)
				temp = 1.0e-10;

			double fi = m_q1 * pow(dt, -(x[2] + dx1[2])) * exp((x[1] + dx1[1]) * temp) + (x[0] + dx1[0]);

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
			double temp = (pow(dt, 1 - (x[2] + dx2[2])) - 1) / (1 - (x[2] + dx2[2]));
			if (temp < 1.0e-10)
				temp = 1.0e-10;

			double fi = m_q1 * pow(dt, -(x[2] + dx2[2])) * exp((x[1] + dx2[1]) * temp) + (x[0] + dx2[0]);

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
				m_qinf = x[0] + dx[0];
				m_a = x[1] + dx[1];
				m_m = x[2] + dx[2];
			}
			else
			{
				R2 = R2_2;
				lambda = lambda / v;
				m_qinf = x[0] + dx2[0];
				m_a = x[1] + dx2[1];
				m_m = x[2] + dx2[2];
			}
		}
		else
		{
			if (R2_1 >= R2_2)
			{
				R2 = R2_1;
				lambda = lambda * v;
				m_qinf = x[0] + dx1[0];
				m_a = x[1] + dx1[1];
				m_m = x[2] + dx1[2];
			}
			else
			{
				R2 = R2_2;
				lambda = lambda / v;
				m_qinf = x[0] + dx2[0];
				m_a = x[1] + dx2[1];
				m_m = x[2] + dx2[2];
			}
		}

		x[0] = m_qinf;
		x[1] = m_a;
		x[2] = m_m;
		i_iter++;   // increase iteration
	}  // end of L-M loop


	//

}



//---------------------------------------------------------------------------
void
DataDuong::Duong_forecast()
{

}



//---------------------------------------------------------------------------
// The following is gaussian elimination 
#define mat_elem(a, y, x, n) (a + ((y) * (n) + (x)))

void
DataDuong::swap_row(double *a, double *b, int r1, int r2, int n)
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
DataDuong::gauss_eliminate(double *a, double *b, double *x, int n)
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

