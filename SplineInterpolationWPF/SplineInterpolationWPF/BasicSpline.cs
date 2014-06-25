using System.Collections.Generic;

namespace SplineInterpolationWPF
{
    public abstract class BasicSpline
    {
        public Cubic[] CalcNaturalCubic(double[] value)
        {
            int num = value.Length - 1;

            double[] gamma = new double[num + 1];
            double[] delta = new double[num + 1];
            double[] d = new double[num + 1];

            int i;
            /*
                 We solve the equation
                [2 1       ] [D[0]]   [3(x[1] - x[0])  ]
                |1 4 1     | |D[1]|   |3(x[2] - x[0])  |
                |  1 4 1   | | .  | = |      .         |
                |    ..... | | .  |   |      .         |
                |     1 4 1| | .  |   |3(x[n] - x[n-2])|
                [       1 2] [D[n]]   [3(x[n] - x[n-1])]
          
                by using row operations to convert the matrix to upper triangular
                and then back sustitution.  The D[i] are the derivatives at the knots.
            */
            gamma[0] = 1.0f / 2.0f;
            for (i = 1; i < num; i++)
            {
                gamma[i] = 1.0f / (4.0f - gamma[i - 1]);
            }
            gamma[num] = 1.0f / (2.0f - gamma[num - 1]);

            double p0 = value[0];
            double p1 = value[0];

            delta[0] = 3.0f * (p1 - p0) * gamma[0];
            for (i = 1; i < num; i++)
            {
                p0 = value[i - 1];
                p1 = value[i + 1];
                delta[i] = (3.0f * (p1 - p0) - delta[i - 1]) * gamma[i];
            }
            p0 = value[num - 1];
            p1 = value[num];

            delta[num] = (3.0f * (p1 - p0) - delta[num - 1]) * gamma[num];

            d[num] = delta[num];
            for (i = num - 1; i >= 0; i--)
            {
                d[i] = delta[i] - gamma[i] * d[i + 1];
            }

            /*
                 now compute the coefficients of the cubics 
            */
            List<Cubic> cubicCollection = new List<Cubic>();

            for (i = 0; i < num; i++)
            {
                p0 = value[i];
                p1 = value[i + 1];

                cubicCollection.Add(new Cubic(
                               p0,
                               d[i],
                               3 * (p1 - p0) - 2 * d[i] - d[i + 1],
                               2 * (p0 - p1) + d[i] + d[i + 1]
                             )
                         );
            }
            return cubicCollection.ToArray();
        }
    }
}
