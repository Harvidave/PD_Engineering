using System;

namespace SplineInterpolationWPF
{
    public class Spline3D : BasicSpline
    {
        private readonly double[] _x;
        private readonly double[] _y;
        private readonly double[] _z;

        private Cubic[] _xCubics;
        private Cubic[] _yCubics;
        private Cubic[] _zCubics;

        public Spline3D(double[] x, double[] y, double[] z)
        {
            _x = x;
            _y = y;
            _z = z;
        }

        public void CalcSpline()
        {
            _xCubics = CalcNaturalCubic(_x);
            _yCubics = CalcNaturalCubic(_y);
            _zCubics = CalcNaturalCubic(_z);
        }

        public Tuple<double, double, double> GetPoint(double position)
        {
            position = position * _xCubics.Length;
            int cubicNum = (int)position;
            double cubicPos = (position - cubicNum);

            return new Tuple<double, double, double>(_xCubics[cubicNum].Eval(cubicPos),
                           _yCubics[cubicNum].Eval(cubicPos),
                           _zCubics[cubicNum].Eval(cubicPos));
        }
    }
}
