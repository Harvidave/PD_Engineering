using System;
using System.Linq;

namespace ProductionDirector.Engineering.Interpolation.Kriging
{
    public class Kriging2D
    {
        private double _nugget;

        private double _range;

        private double _sill;

        private double[] _k;

        private double[] _m;

        private int _n;

        public double[] TargetValue { get; set; }

        public double[] XCoord { get; set; }

        public double[] YCoord { get; set; }

        public double Sigma2 { get; set; }

        public double Alpha { get; set; }

        public double A { get; set; }

        public KrigingModel Model { get; set; }

        public void Train()
        {
            int i, j, k, l;
            // Lag distance/semivariance
            int targetLength = TargetValue.Length;
            int distanceLength = (targetLength*targetLength - targetLength)/2;
            double[][] distance = new double[distanceLength][];
            for (i = 0, k = 0; i < targetLength; i++)
                for (j = 0; j < i; j++, k++)
                {
                    distance[k] = new double[2];
                    distance[k][0] = Math.Pow(
                        Math.Pow(XCoord[i] - XCoord[j], 2) +
                        Math.Pow(YCoord[i] - YCoord[j], 2), 0.5);
                    distance[k][1] = Math.Abs(TargetValue[i] - TargetValue[j]);
                }
            distance = distance.OrderBy(doubles => doubles[0]).ToArray();
            _range = distance[distanceLength - 1][0];

            // Bin lag distance
            int lags = (distanceLength) > 30 ? 30 : distanceLength;
            double tolerance = _range/lags;
            double[] lag = new double[lags];
            double[] semi = new double[lags];
            if (lags < 30)
            {
                for (l = 0; l < lags; l++)
                {
                    lag[l] = distance[l][0];
                    semi[l] = distance[l][1];
                }
            }
            else
            {
                for (i = 0, j = 0, k = 0, l = 0; i < lags && j < (distanceLength); i++, k = 0)
                {
                    while (distance[j][0] <= ((i + 1)*tolerance))
                    {
                        lag[l] += distance[j][0];
                        semi[l] += distance[j][1];
                        j++;
                        k++;
                        if (j >= (distanceLength))
                        {
                            break;
                        }
                    }
                    if (k > 0)
                    {
                        lag[l] /= k;
                        semi[l] /= k;
                        l++;
                    }
                }
                if (l < 2)
                {
                    return;
                } // Error: Not enough points
            }

            // Feature transformation
            int currentL = l;
            _range = lag[currentL - 1] - lag[0];
            double[] x = Enumerable.Repeat(1.0, 2*currentL).ToArray();
            double[] y = new double[currentL];
            for (i = 0; i < currentL; i++)
            {
                switch (Model)
                {
                    case KrigingModel.Gaussian:
                        x[i*2 + 1] = 1.0 - Math.Exp(-(1.0/A)*Math.Pow(lag[i]/_range, 2));
                        break;
                    case KrigingModel.Exponential:
                        x[i*2 + 1] = 1.0 - Math.Exp(-(1.0/A)*lag[i]/_range);
                        break;
                    case KrigingModel.Spherical:
                        x[i*2 + 1] = 1.5*(lag[i]/_range) -
                                     0.5*Math.Pow(lag[i]/_range, 3);
                        break;
                }
                y[i] = semi[i];
            }

            // Least squares
            double[] xt = KrigingMatrixTranspose(x, currentL, 2);
            double[] z = KrigingMatrixMultiply(xt, x, 2, currentL, 2);
            z = KrigingMatrixAdd(z, KrigingMatrixDiag(1/Alpha, 2), 2, 2);
            double[] cloneZ = {z[0]};
            if (KrigingMatrixChol(z, 2))
            {
                KrigingMatrixChol2Inv(z, 2);
            }
            else
            {
                KrigingMatrixSolve(cloneZ, 2);
                z = cloneZ;
            }
            double[] w = KrigingMatrixMultiply(KrigingMatrixMultiply(z, xt, 2, 2, currentL), y, 2, currentL, 1);

            // Variogram parameters
            _nugget = w[0];
            _sill = w[1]*_range + _nugget;
            _n = XCoord.Length;

            // Gram matrix with prior
            currentL = XCoord.Length;
            double[] kk = new double[currentL*currentL];
            for (i = 0; i < currentL; i++)
            {
                for (j = 0; j < i; j++)
                {
                    kk[i*currentL + j] =
                        CalculateProbablity(
                            Math.Pow(Math.Pow(XCoord[i] - XCoord[j], 2) + Math.Pow(YCoord[i] - YCoord[j], 2), 0.5));
                    kk[j*currentL + i] = kk[i*currentL + j];
                }
                kk[i*currentL + i] = CalculateProbablity(0);
            }

            // Inverse penalized Gram matrix projected to target vector
            double[] c = KrigingMatrixAdd(kk, KrigingMatrixDiag(Sigma2, currentL), currentL, currentL);
            double[] cloneC = c.Clone() as double[];
            if (KrigingMatrixChol(c, currentL))
            {
                KrigingMatrixChol2Inv(c, currentL);
            }
            else
            {
                KrigingMatrixSolve(cloneC, currentL);
                c = cloneC;
            }

            // Copy unprojected inverted matrix as K
            if (c != null)
            {
                _k = c.Clone() as double[];
            }
            double[] m = KrigingMatrixMultiply(c, TargetValue, currentL, currentL, 1);
            _m = m;
        }

        public double Predict(double x, double y)
        {
            double[] k = new double[_n];
            for (int i = 0; i < _n; i++)
            {
                k[i] = CalculateProbablity(Math.Pow(Math.Pow(x - XCoord[i], 2) + Math.Pow(y - YCoord[i], 2), 0.5));
            }
            return KrigingMatrixMultiply(k, _m, 1, _n, 1)[0];
        }

	    public double[,] Calculate(double xMin, double xMax, double yMin, double yMax, int dimension,
		    CoordinateSystemOrigin origin)
	    {
		    double[,] result = new double[dimension, dimension];
		    double xRange = xMax - xMin;
		    double yRange = yMax - yMin;
		    if (origin == CoordinateSystemOrigin.LeftUpper)
		    {
			    for (int i = 0; i < dimension; i++)
			    {
				    for (int j = 0; j < dimension; j++)
				    {
					    result[i, j] = Predict(i*xRange/(dimension - 1), j*yRange/(dimension - 1));
				    }
			    }
		    }
		    else if (origin == CoordinateSystemOrigin.LeftLower)
		    {
			    for (int i = 0; i < dimension; i++)
			    {
				    for (int j = dimension - 1; j >= 0; j--)
				    {
					    result[i, j] = Predict(i*xRange/(dimension - 1), j*yRange/(dimension - 1));
				    }
			    }
		    }
		    return result;
	    }

	    public double Variance(double x, double y)
        {
            double[] k = new double[_n];
            for (int i = 0; i < _n; i++)
            {
                k[i] = CalculateProbablity(Math.Pow(Math.Pow(x - XCoord[i], 2) + Math.Pow(y - YCoord[i], 2), 0.5));
            }
            return CalculateProbablity(0) +
                   KrigingMatrixMultiply(KrigingMatrixMultiply(k, _k, 1, _n, _n), k, 1, _n, 1)[0];
        }

        private static void KrigingMatrixSolve(double[] x, int n)
        {
            int m = n;
            double[] b = new double[n*n];
            int[] indxc = new int[n];
            int[] indxr = new int[n];
            int[] ipiv = new int[n];
            int i, j;

            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                {
                    if (i == j)
                    {
                        b[i*n + j] = 1;
                    }
                    else b[i*n + j] = 0;
                }
            for (j = 0; j < n; j++)
            {
                ipiv[j] = 0;
            }
            int l;
            int k;
            double temp;
            for (i = 0; i < n; i++)
            {
                double big = 0;
                int icol = 0;
                int irow = 0;
                for (j = 0; j < n; j++)
                {
                    if (ipiv[j] != 1)
                    {
                        for (k = 0; k < n; k++)
                        {
                            if (ipiv[k] == 0)
                            {
                                if (Math.Abs(x[j*n + k]) >= big)
                                {
                                    big = Math.Abs(x[j*n + k]);
                                    irow = j;
                                    icol = k;
                                }
                            }
                        }
                    }
                }
                ++(ipiv[icol]);

                if (irow != icol)
                {
                    for (l = 0; l < n; l++)
                    {
                        temp = x[irow*n + l];
                        x[irow*n + l] = x[icol*n + l];
                        x[icol*n + l] = temp;
                    }
                    for (l = 0; l < m; l++)
                    {
                        temp = b[irow*n + l];
                        b[irow*n + l] = b[icol*n + l];
                        b[icol*n + l] = temp;
                    }
                }
                indxr[i] = irow;
                indxc[i] = icol;

                if (Math.Abs(x[icol*n + icol]) < Double.Epsilon)
                {
                    return;
                } // Singula

                double pivinv = 1/x[icol*n + icol];
                x[icol*n + icol] = 1;
                for (l = 0; l < n; l++) x[icol*n + l] *= pivinv;
                for (l = 0; l < m; l++) b[icol*n + l] *= pivinv;

                for (int ll = 0; ll < n; ll++)
                {
                    if (ll != icol)
                    {
                        double dum = x[ll*n + icol];
                        x[ll*n + icol] = 0;
                        for (l = 0; l < n; l++) x[ll*n + l] -= x[icol*n + l]*dum;
                        for (l = 0; l < m; l++) b[ll*n + l] -= b[icol*n + l]*dum;
                    }
                }
            }
            for (l = (n - 1); l >= 0; l--)
                if (indxr[l] != indxc[l])
                {
                    for (k = 0; k < n; k++)
                    {
                        temp = x[k*n + indxr[l]];
                        x[k*n + indxr[l]] = x[k*n + indxc[l]];
                        x[k*n + indxc[l]] = temp;
                    }
                }
        }

        private static void KrigingMatrixChol2Inv(double[] x, int n)
        {
            int i, j, k;
            for (i = 0; i < n; i++)
            {
                x[i*n + i] = 1/x[i*n + i];
                for (j = i + 1; j < n; j++)
                {
                    double sum = 0.0;
                    for (k = i; k < j; k++)
                        sum -= x[j*n + k]*x[k*n + i];
                    x[j*n + i] = sum/x[j*n + j];
                }
            }
            for (i = 0; i < n; i++)
                for (j = i + 1; j < n; j++)
                    x[i*n + j] = 0;
            for (i = 0; i < n; i++)
            {
                x[i*n + i] *= x[i*n + i];
                for (k = i + 1; k < n; k++)
                    x[i*n + i] += x[k*n + i]*x[k*n + i];
                for (j = i + 1; j < n; j++)
                    for (k = j; k < n; k++)
                        x[i*n + j] += x[k*n + i]*x[k*n + j];
            }
            for (i = 0; i < n; i++)
                for (j = 0; j < i; j++)
                    x[i*n + j] = x[j*n + i];
        }

        private static bool KrigingMatrixChol(double[] x, int n)
        {
            int i;
            double[] p = new double[n];
            for (i = 0; i < n; i++) p[i] = x[i*n + i];
            for (i = 0; i < n; i++)
            {
                int j;
                for (j = 0; j < i; j++)
                    p[i] -= x[i*n + j]*x[i*n + j];
                if (p[i] <= 0)
                {
                    return false;
                }
                p[i] = Math.Sqrt(p[i]);
                for (j = i + 1; j < n; j++)
                {
                    for (int k = 0; k < i; k++)
                        x[j*n + i] -= x[j*n + k]*x[i*n + k];
                    x[j*n + i] /= p[i];
                }
            }
            for (i = 0; i < n; i++) x[i*n + i] = p[i];
            return true;
        }

        private static double[] KrigingMatrixAdd(double[] x, double[] y, int n, int m)
        {
            double[] z = new double[n*m];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                    z[i*m + j] = x[i*m + j] + y[i*m + j];
            return z;
        }

        private static double[] KrigingMatrixDiag(double c, int n)
        {
            double[] z = new double[n*n];
            for (int i = 0; i < n; i++) z[i*n + i] = c;
            return z;
        }

        private static double[] KrigingMatrixMultiply(double[] x, double[] y, int n, int m, int p)
        {
            double[] z = new double[n*p];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < p; j++)
                {
                    z[i*p + j] = 0;
                    for (int k = 0; k < m; k++)
                        z[i*p + j] += x[i*m + k]*y[k*p + j];
                }
            }
            return z;
        }

        private static double[] KrigingMatrixTranspose(double[] x, int n, int m)
        {
            double[] z = new double[m*n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    z[j*n + i] = x[i*m + j];
                }
            }
            return z;
        }

        private static void KrigingMatrixScale(double[] x, double c, int n, int m)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    x[i*m + j] *= c;
                }
            }
        }

        private double CalculateProbablity(double h)
        {
            return CalculateProbablity(h, _nugget, _range, _sill, A);
        }

        private double CalculateProbablity(double h, double nugget, double range, double sill, double a)
        {
            switch (Model)
            {
                case KrigingModel.Gaussian:
                    return nugget + ((sill - nugget)/range)*(1.0 - Math.Exp(-(1.0/a)*Math.Pow(h/range, 2)));
                case KrigingModel.Exponential:
                    return nugget + ((sill - nugget)/range)*(1.0 - Math.Exp(-(1.0/a)*(h/range)));
                case KrigingModel.Spherical:
                    if (h > range)
                    {
                        return nugget + (sill - nugget)/range;
                    }
                    return nugget + ((sill - nugget)/range)*(1.5*(h/range) - 0.5*Math.Pow(h/range, 3));
            }
            return 0.0;
        }
    }
}
