using System;
using System.Runtime.InteropServices;

namespace ProductionDirector.Engineering.Forecast.DataAnalysis
{
    public class Duong
    {
        [DllImport("ProductionDirector.Engineering.Forecast.DataAnalysisCppWrapper.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern IntPtr ComputeDuongForecast(int forecastMethod,
            [MarshalAs(UnmanagedType.LPArray, ArraySubType = UnmanagedType.LPWStr)] double[] x,
            [MarshalAs(UnmanagedType.LPArray, ArraySubType = UnmanagedType.LPWStr)] double[] y, int inputLength,
            [MarshalAs(UnmanagedType.LPArray, ArraySubType = UnmanagedType.LPWStr)] double[] future, int futureLength);

        [DllImport("ProductionDirector.Engineering.Forecast.DataAnalysisCppWrapper.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern double ComputeDuongEur(int forecastMethod,
            [MarshalAs(UnmanagedType.LPArray, ArraySubType = UnmanagedType.LPWStr)] double[] x,
            [MarshalAs(UnmanagedType.LPArray, ArraySubType = UnmanagedType.LPWStr)] double[] y, int inputLength,
            double tf);

        [DllImport("ProductionDirector.Engineering.Forecast.DataAnalysisCppWrapper.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern int ReleaseMemory(IntPtr intptr);

        public static double[] ComputeForecast(DuongMethodEnum method, double[] x, double[] y, double[] future)
        {
            if (x.Length != y.Length)
            {
                return null;
            }

            var result = new double[future.Length];
            IntPtr computeForecast = ComputeDuongForecast((int)method, x, y, x.Length, future, future.Length);
            Marshal.Copy(computeForecast, result, 0, future.Length);
            ReleaseMemory(computeForecast);
            return result;
        }
        public static double ComputeEur(DuongMethodEnum method, double[] x, double[] y, double tf)
        {
            if (x.Length != y.Length)
            {
                return double.MinValue;
            }

            return ComputeDuongEur((int) method, x, y, x.Length, tf);
        }
    }
}
