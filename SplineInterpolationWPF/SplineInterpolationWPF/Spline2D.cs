//using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Text;
//using System.Threading.Tasks;

//namespace SplineInterpolationWPF
//{
//    public class Spline2D : BasicSpline
//    {
//        private readonly double[] _x;
//        private readonly double[] _y;

//        private Cubic[] _xCubics;
//        private Cubic[] _yCubics;

//        public Spline2D(double[] x, double[] y)
//        {
//            _x = x;
//            _y = y;
//        }

//        public void CalcSpline()
//        {
//            _xCubics = CalcNaturalCubic(_x);
//            _yCubics = CalcNaturalCubic(_y);
//        }

//        public Vector2f GetPoint(double position)
//        {
//            position = position * _xCubics.Length; // extrapolate to the arraysize
//            int cubicNum = (int)position;
//            double cubicPos = (position - cubicNum);

//            return new Vector2f(_xCubics.get(cubicNum).eval(cubicPos),
//                           _yCubics.get(cubicNum).eval(cubicPos));
//        }
//    }
//}
