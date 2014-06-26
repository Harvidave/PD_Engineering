using System.Collections.Generic;
using System.Windows;

namespace PolygonWPF
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

        //Polygon: 0.0,0.0;10.0,0.0;10.0,10.0;0.0,10.0
        //Point: 20.0,20.0 - No
        //Point: 5.0,5.0 - Yes

        //Polygon: 0.0,0.0;5.0,5.0;5.0,0.0
        //Point: 3.0,3.0 - Yes
        //Point: 5.0,1.0 - Yes
        //Point: 8.0,1.0 - No

        //Polygon: 0.0,0.0;10.0,0.0;10.0,10.0;0.0,10.0
        //Point: -1.0,10.0 - No

        //Polygon: 1.0,4.0;1.0,6.0;2.0,7.0;4.0,7.0;5.0,6.0;2.0,3.0
        //Point: 3.0,6.0 - Yes

        //Polygon: 2.0,-2.0;4.0,0.0;2.0,2.0;0.0,0.0
        //Point: 2.0,0.0 - Yes

        private void CheckClick(object sender, RoutedEventArgs e)
        {
            List<Point> polygonList = new List<Point>();
            foreach (var s in PolygonText.Text.Split(';'))
            {
                var ploygonPoints = s.Split(',');
                polygonList.Add(new Point(double.Parse(ploygonPoints[0]), double.Parse(ploygonPoints[1])));
            }

            var points = PointText.Text.Split(',');
            Point p = new Point(double.Parse(points[0]), double.Parse(points[1]));

            ResultText.Text = PolygonUtility.IsInside(PolygonUtility.GetClockWisedPoints(polygonList), p) ? "Yes" : "No";
        }
    }
}
