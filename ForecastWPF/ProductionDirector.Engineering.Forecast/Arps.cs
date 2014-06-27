using System;
using System.Runtime.InteropServices;

namespace ProductionDirector.Engineering.Forecast
{
    public class Arps
    {
        [DllImport("ProductionDirectorForecast.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern IntPtr ComputeArpsForecast(int forecastMethod,
            [MarshalAs(UnmanagedType.LPArray, ArraySubType = UnmanagedType.LPWStr)] double[] x,
            [MarshalAs(UnmanagedType.LPArray, ArraySubType = UnmanagedType.LPWStr)] double[] y, int inputLength,
            [MarshalAs(UnmanagedType.LPArray, ArraySubType = UnmanagedType.LPWStr)] double[] future, int futureLength);

        [DllImport("ProductionDirectorForecast.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern double ComputeArpsEur(int forecastMethod,
            [MarshalAs(UnmanagedType.LPArray, ArraySubType = UnmanagedType.LPWStr)] double[] x,
            [MarshalAs(UnmanagedType.LPArray, ArraySubType = UnmanagedType.LPWStr)] double[] y, int inputLength,
            double Qi, double qf);

        [DllImport("ProductionDirectorForecast.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern int ReleaseMemory(IntPtr intptr);

        public static double[] ComputeForecast(ArpsMethodEnum method, double[] x, double[] y, double[] future)
        {
            if (x.Length != y.Length)
            {
                return null;
            }

            var result = new double[future.Length];
            IntPtr computeForecast = ComputeArpsForecast((int)method, x, y, x.Length, future, future.Length);
            Marshal.Copy(computeForecast, result, 0, future.Length);
            ReleaseMemory(computeForecast);
            return result;
        }

        public static double ComputeEur(ArpsMethodEnum method, double[] x, double[] y, double Qi, double qf)
        {
            if (x.Length != y.Length)
            {
                return double.MinValue;
            }

            return ComputeArpsEur((int)method, x, y, x.Length, Qi, qf);
        }
    }
}
