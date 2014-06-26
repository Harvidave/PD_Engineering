using System;
using System.Collections.Generic;
using System.Windows;

namespace PolygonWPF
{
    public static class PolygonUtility
    {
        public static Point[] GetClockWisedPoints(List<Point> points)
        {
            List<Point> pointList = new List<Point>(points);
            Point center = GetCenter(pointList);
            pointList.Sort((point1, point2) =>
            {
                if (point1 == point2)
                {
                    return 0;
                }
                if (IsLess(point1, point2, center))
                {
                    return -1;
                }
                return 1;
            });
            return pointList.ToArray();
        }

        private static bool IsLess(Point a, Point b, Point center)
        {
            if (a.X >= 0 && b.X < 0)
            {
                return true;
            }
            if (Math.Abs(a.X) < double.Epsilon && Math.Abs(b.X) < Double.Epsilon)
            {
                return a.Y > b.Y;
            }

            double det = (a.X - center.X) * (b.Y - center.Y) - (b.X - center.X) * (a.Y - center.Y);
            if (det < 0)
            {
                return true;
            }
            if (det > 0)
            {
                return false;
            }

            double d1 = (a.X - center.X) * (a.X - center.X) + (a.Y - center.Y) * (a.Y - center.Y);
            double d2 = (b.X - center.X) * (b.X - center.X) + (b.Y - center.Y) * (b.Y - center.Y);
            return d1 > d2;
        }

        private static Point GetCenter(List<Point> points)
        {
            Point pSum = new Point(0.0, 0.0);
            foreach (Point p in points)
            {
                pSum.X += p.X;
                pSum.Y += p.Y;
            }
            return new Point(pSum.X / points.Count, pSum.Y / points.Count);
        }

        /// <summary>
        /// Given three colinear points p, q, r, the function checks if point q lies on line segment 'pr'
        /// </summary>
        /// <param name="p"></param>
        /// <param name="q">Point to check</param>
        /// <param name="r"></param>
        private static bool OnSegment(Point p, Point q, Point r)
        {
            if (q.X <= Math.Max(p.X, r.X) && q.X >= Math.Min(p.X, r.X) &&
                    q.Y <= Math.Max(p.Y, r.Y) && q.Y >= Math.Min(p.Y, r.Y))
                return true;
            return false;
        }

        /// <summary>
        /// To find orientation of ordered triplet (p, q, r).
        /// The function returns following values
        /// 0 --> p, q and r are colinear
        /// 1 --> Clockwise
        /// 2 --> Counterclockwise
        /// </summary>
        private static int Orientation(Point p, Point q, Point r)
        {
            double val = (q.Y - p.Y) * (r.X - q.X) - (q.X - p.X) * (r.Y - q.Y);

            if ((int)val == 0) return 0;  // colinear
            return (val > 0) ? 1 : 2; // clock or counterclock wise
        }

        /// <summary>
        /// The function that returns true if line segment 'p1q1' and 'p2q2' intersect.
        /// </summary>
        private static bool IsIntersect(Point p1, Point q1, Point p2, Point q2)
        {
            // Find the four orientations needed for general and
            // special cases
            int o1 = Orientation(p1, q1, p2);
            int o2 = Orientation(p1, q1, q2);
            int o3 = Orientation(p2, q2, p1);
            int o4 = Orientation(p2, q2, q1);

            // General case

            if (o1 != o2 && o3 != o4)
            {
                if (o4 == 0)
                {
                    return false;
                }
                if (o3 == 0)
                {
                    if (o4 != 1)
                    {
                        return false;
                    }
                }
                return true;
            }

            // Special Cases
            // p1, q1 and p2 are colinear and p2 lies on segment p1q1
            if (o1 == 0 && OnSegment(p1, p2, q1)) return true;

            // p1, q1 and p2 are colinear and q2 lies on segment p1q1
            if (o2 == 0 && OnSegment(p1, q2, q1)) return true;

            // p2, q2 and p1 are colinear and p1 lies on segment p2q2
            if (o3 == 0 && OnSegment(p2, p1, q2)) return true;

            // p2, q2 and q1 are colinear and q1 lies on segment p2q2
            if (o4 == 0 && OnSegment(p2, q1, q2)) return true;

            return false; // Doesn't fall in any of the above cases
        }

        /// <summary>
        ///Returns true if the point p lies inside the polygon[] with n vertices
        /// </summary>
        public static bool IsInside(Point[] polygon, Point p)
        {
            // There must be at least 3 vertices in polygon[]
            if (polygon.Length < 3) return false;

            // Create a point for line segment from p to infinite
            Point extreme = new Point(double.MaxValue, p.Y);

            // Count intersections of the above line with sides of polygon
            int count = 0, i = 0;
            do
            {
                int next = (i + 1) % polygon.Length;

                // Check if the line segment from 'p' to 'extreme' intersects
                // with the line segment from 'polygon[i]' to 'polygon[next]'
                if (IsIntersect(polygon[i], polygon[next], p, extreme))
                {
                    // If the point 'p' is colinear with line segment 'i-next',
                    // then check if it lies on segment. If it lies, return true,
                    // otherwise false
                    if (Orientation(polygon[i], p, polygon[next]) == 0)
                        return OnSegment(polygon[i], p, polygon[next]);

                    count++;
                }
                i = next;
            } while (i != 0);

            // Return true if count is odd, false otherwise
            return count % 2 == 1;
        }
    }
}
