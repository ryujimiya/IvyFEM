using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; // DllImport

namespace IvyFEM.FFT
{
    public class ImportedFunctions
    {
        ///////////////////////////////////////////////////////////////////
        // OouraFFT
        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_cdft4g(int n, int isgn, double* a, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_rdft4g(int n, int isgn, double* a, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_ddct4g(int n, int isgn, double* a, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_ddst4g(int n, int isgn, double* a, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_dfct4g(int n, double* a, double* t, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_dfst4g(int n, double* a, double* t, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_cdft8g(int n, int isgn, double* a, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_rdft8g(int n, int isgn, double* a, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_ddct8g(int n, int isgn, double* a, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_ddst8g(int n, int isgn, double* a, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_dfct8g(int n, double* a, double* t, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_dfst8g(int n, double* a, double* t, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_cdftsg(int n, int isgn, double* a, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_rdftsg(int n, int isgn, double* a, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_ddctsg(int n, int isgn, double* a, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_ddstsg(int n, int isgn, double* a, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_dfctsg(int n, double* a, double* t, int* ip, double* w);

        [DllImport("OouraFFTDll.dll")]
        public static extern unsafe void OouraFFT_dfstsg(int n, double* a, double* t, int* ip, double* w);
    }
}
