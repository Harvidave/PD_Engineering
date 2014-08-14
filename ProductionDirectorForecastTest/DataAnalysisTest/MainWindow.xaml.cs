using System;
using System.Linq;
using System.Windows;
using ProductionDirector.Engineering.Forecast.DataAnalysis;

namespace ProductionDirectorForecastTest
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

        private void ArpsClick(object sender, RoutedEventArgs e)
        {
            // Rate_Time x: 0.0,1.0,2.0,3.0,4.0
            // Rate_Time y: 30.0,28.0,26.0,24.0,20.0
            // Rate_Time future: 5.0,6.0,7.0,8.0,9.0,10.0
            // Rate_Q x: 0.0,10.0,20.0,30.0,40.0
            // Rate_Q y: 100.0,73.3904,52.77319,37.07398,25.35525
            // Rate_Q future: 50.0,60.0,70.0,80.0,90.0,100.0

            double[] x = X.Text.Split(',').Select(double.Parse).ToArray();
            double[] y = Y.Text.Split(',').Select(double.Parse).ToArray();
            double[] future = Future.Text.Split(',').Select(double.Parse).ToArray();
            int method = int.Parse(Method.Text);

            double[] result = Arps.ComputeForecast((ArpsMethodEnum)method, x, y, future);
            Result.Text = string.Empty;
            foreach (double d in result)
            {
                Result.Text = string.Concat(Result.Text, d, Environment.NewLine);
            }

            var computeEur = Arps.ComputeEur((ArpsMethodEnum)method, x, y, 100.0, 0.0001);
            Eur.Text = computeEur.ToString();
        }

        private void DuongClick(object sender, RoutedEventArgs e)
        {
            // x: 1.0,2.0,3.0,4.0,5.0
            // y: 50.0,44.21524,38.9416,35.28946,32.581
            // z: 6.0,7.0,8.0,9.0,10.0,11.0

            double[] x = X.Text.Split(',').Select(double.Parse).ToArray();
            double[] y = Y.Text.Split(',').Select(double.Parse).ToArray();
            double[] future = Future.Text.Split(',').Select(double.Parse).ToArray();
            int method = int.Parse(Method.Text);

            double[] result = Duong.ComputeForecast((DuongMethodEnum)method, x, y, future);
            Result.Text = string.Empty;
            foreach (double d in result)
            {
                Result.Text = string.Concat(Result.Text, d, Environment.NewLine);
            }

            Eur.Text = Duong.ComputeEur((DuongMethodEnum)method, x, y, 20.0).ToString();
        }

        private void SEPDClick(object sender, RoutedEventArgs e)
        {
            // x: 1.0,2.0,3.0,4.0,5.0  --  Cannot start with 0.0
            // y: 100.0,73.3904,52.77319,37.07398,25.35525
            // z: 6.0,7.0,8.0,9.0,10.0,11.0
            double[] x = X.Text.Split(',').Select(double.Parse).ToArray();
            double[] y = Y.Text.Split(',').Select(double.Parse).ToArray();
            double[] future = Future.Text.Split(',').Select(double.Parse).ToArray();

            double[] result = SEPD.ComputeForecast(x, y, future);
            Result.Text = string.Empty;
            foreach (double d in result)
            {
                Result.Text = string.Concat(Result.Text, d, Environment.NewLine);
            }

            Eur.Text = SEPD.ComputeEur(x, y).ToString();
        }
    }
}
