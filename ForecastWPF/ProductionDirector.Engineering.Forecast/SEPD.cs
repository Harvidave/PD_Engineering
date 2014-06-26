using System;
using System.Runtime.InteropServices;

namespace ProductionDirector.Engineering.Forecast
{
    public class SEPD
    {
        [DllImport("ProductionDirectorForecast.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern IntPtr ComputeSEPD(
            [MarshalAs(UnmanagedType.LPArray, ArraySubType = UnmanagedType.LPWStr)] double[] x,
            [MarshalAs(UnmanagedType.LPArray, ArraySubType = UnmanagedType.LPWStr)] double[] y, int inputLength,
            [MarshalAs(UnmanagedType.LPArray, ArraySubType = UnmanagedType.LPWStr)] double[] future, int futureLength);

        [DllImport("ProductionDirectorForecast.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern int ReleaseMemory(IntPtr intptr);

        public static double[] ComputeForecast(double[] x, double[] y, double[] future)
        {
            if (x.Length != y.Length)
            {
                return null;
            }

            var result = new double[future.Length];
            IntPtr computeForecast = ComputeSEPD(x, y, x.Length, future, future.Length);
            Marshal.Copy(computeForecast, result, 0, future.Length);
            ReleaseMemory(computeForecast);
            return result;
        }
    }
}
