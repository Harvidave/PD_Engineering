namespace SplineInterpolationWPF
{
    public class Cubic
    {
        private readonly double _a;
        private readonly double _b;
        private readonly double _c;
        private readonly double _d;

        public Cubic(double a, double b, double c, double d)
        {
            _a = a;
            _b = b;
            _c = c;
            _d = d;
        }

        public double Eval(double u)
        {
            return (((_d * u) + _c) * u + _b) * u + _a;
        }
    }
}
