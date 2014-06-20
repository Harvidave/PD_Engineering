using System;
using System.Linq;
using System.Windows;
using ProductionDirector.Engineering.Forecast;

namespace ForecastWPF
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
            double[] x = X.Text.Split(',').Select(double.Parse).ToArray();
            double[] y = Y.Text.Split(',').Select(double.Parse).ToArray();
            double[] future = Future.Text.Split(',').Select(double.Parse).ToArray();
            int method = int.Parse(Method.Text);

            double[] result = Arps.ComputeForecast((ArpsMethodEnum) method, x, y, future);
            Result.Text = string.Empty;
            foreach (double d in result)
            {
                Result.Text = string.Concat(Result.Text, d, Environment.NewLine);
            }
        }
    }
}
