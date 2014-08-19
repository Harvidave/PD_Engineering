using System;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using ProductionDirector.Engineering.Interpolation.Kriging;

namespace ProductionDirectorInterpolationTest
{
	/// <summary>
	/// Interaction logic for MainWindow.xaml
	/// </summary>
	public partial class MainWindow
	{
		private KrigingModel _model;

		public MainWindow()
		{
			InitializeComponent();

			Spherical.IsChecked = true;
		}

		private void CalculateClick(object sender, RoutedEventArgs e)
		{
			double[] targets = Val.Text.Split(',').Select(double.Parse).ToArray();
			double[] x = X.Text.Split(',').Select(double.Parse).ToArray();
			double[] y = Y.Text.Split(',').Select(double.Parse).ToArray();

			var v = new Kriging2D
			{
				TargetValue = targets,
				XCoord = x,
				YCoord = y,
				Model = _model,
				Sigma2 = double.Parse(Sigma2.Text),
				Alpha = double.Parse(Alpha.Text),
				A = double.Parse(A.Text)
			};
			v.Train();

			string result = string.Empty;
			double xRange = x.Max() + 2.0;
			double yRange = y.Max() + 2.0;
			int gridDimension = int.Parse(Dimension.Text);
			double[,] calculate = v.Calculate(0.0, xRange, 0.0, yRange, gridDimension, CoordinateSystemOrigin.LeftUpper);
			for (int i = 0; i < Math.Sqrt(calculate.Length); i++)
			{
				for (int j = 0; j < Math.Sqrt(calculate.Length); j++)
				{

					result = string.Concat(result, calculate[j, i], ", ");
				}
				result = string.Concat(result.Remove(result.Length - 1), Environment.NewLine);
			}
			//for (int i = 0; i < gridDimension; i++)
			//{
			//	for (int j = 0; j < gridDimension; j++)
			//	{
			//		result = string.Concat(result, v.Predict(j * xRange / (gridDimension - 1), i * yRange / (gridDimension - 1)),
			//			", ");
			//	}
			//	result = string.Concat(result.Remove(result.Length - 1), Environment.NewLine);
			//}
			string txtFile = string.Format("Result{0}.txt", DateTime.Now.Ticks);
			File.WriteAllText(txtFile, result);
			Txt.Text = txtFile;
			string csvFile = string.Format("Result{0}.csv", DateTime.Now.Ticks);
			File.WriteAllText(csvFile, result);
			Excel.Text = csvFile;
		}

		private void ResetClick(object sender, RoutedEventArgs e)
		{
		}

		private void GaussianChecked(object sender, RoutedEventArgs e)
		{
			SetModel(sender);
		}

		private void ExponentialChecked(object sender, RoutedEventArgs e)
		{
			SetModel(sender);
		}

		private void SphericalChecked(object sender, RoutedEventArgs e)
		{
			SetModel(sender);
		}

		private void SetModel(object sender)
		{
			var radioButton = sender as RadioButton;
			if (radioButton != null)
			{
				_model = (KrigingModel)Enum.Parse(typeof(KrigingModel), radioButton.Content.ToString());
			}
		}

		private void OpenTxtClick(object sender, RoutedEventArgs e)
		{
			Process.Start(Txt.Text);
		}

		private void OpenExcelClick(object sender, RoutedEventArgs e)
		{
			Process.Start(Excel.Text);
		}
	}
}
