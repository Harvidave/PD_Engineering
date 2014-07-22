How to use it:


            double[] targets = Val.Text.Split(',').Select(double.Parse).ToArray();
            double[] x = X.Text.Split(',').Select(double.Parse).ToArray();
            double[] y = Y.Text.Split(',').Select(double.Parse).ToArray();

            var v = new Variogram
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
            double xRange = x.Max() - x.Min();
            double yRange = y.Max() - y.Min();
            int gridDimension = int.Parse(Dimension.Text);
            for (int i = 0; i < gridDimension; i++)
            {
                for (int j = 0; j < gridDimension; j++)
                {
                    result = string.Concat(result, v.Predict(j * xRange / (gridDimension - 1), i * yRange / (gridDimension - 1)),
                        ", ");
                }
                result = string.Concat(result.Remove(result.Length - 1), Environment.NewLine);
            }
            string txtFile = string.Format("Result{0}.txt", DateTime.Now.Ticks);
            File.WriteAllText(txtFile, result);
            Txt.Text = txtFile;
            string csvFile = string.Format("Result{0}.csv", DateTime.Now.Ticks);
            File.WriteAllText(csvFile, result);
            Excel.Text = csvFile;