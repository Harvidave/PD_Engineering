using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows;

namespace SplineInterpolationWPF
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        private void ComputeClick(object sender, RoutedEventArgs e)
        {
            double[] xx = {5.0, 2.0, 9.0};
            double[] yy = {7.0, 11.0, 4.0};
            double[] zz = {8.0, 6.0, 3.0};
            Spline3D s = new Spline3D(xx, yy, zz);
            s.CalcSpline();
            var point = s.GetPoint(0.6667);

            return;


            const int r = 1000;
            var known = new Dictionary<double, double>
                { 
                    { 0.0, 0.0 },
                    { 100.0, 0.50 * r },
                    { 300.0, 0.75 * r },
                    { 500.0, 1.00 * r },
                };
            foreach (var pair in known)
            {
                Console.WriteLine("{0:0.000}\t{1:0.000}", pair.Key, pair.Value);
            }

            var scaler = new SplineInterpolator(known);
            var start = known.First().Key;
            var end = known.Last().Key;
            var step = (end - start) / 50;

            for (var x = start; x <= end; x += step)
            {
                if (Math.Abs(x - start) < double.Epsilon)
                {
                    continue;
                }
                var y = scaler.GetValue(x);
                Console.WriteLine("\t\t{0:0.000}\t{1:0.000}", x, y);
            }
        }
    }
}
