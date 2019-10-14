using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.FFT
{
    public class Functions
    {
        /// <summary>
        /// 界を高速フーリエ変換する
        /// </summary>
        /// <param name="times"></param>
        /// <param name="timeDomainDatas"></param>
        /// <param name="freqs"></param>
        /// <param name="freqDomaindatas"></param>
        public static void DoFFT(
            double[] times,
            double[] timeDomainDatas,
            out double[] freqs,
            out System.Numerics.Complex[] freqDomaindatas)
        {
            freqs = null;
            freqDomaindatas = null;

            System.Diagnostics.Debug.Assert(times.Length == timeDomainDatas.Length);
            // Note: dataCnt は 2^mでなければならない
            int dataCnt = timeDomainDatas.Length;
            double dt = times[1] - times[0];

            // FFTを実行する
            {
                int n = dataCnt;
                int isgn = 1; // 1 or -1;
                double[] a = new double[n * 2];
                /*
                a[0...2 * n - 1]   :input / output data(double *)
                        input data
                            a[2 * j] = Re(x[j]), 
                            a[2 * j + 1] = Im(x[j]), 0 <= j < n
                        output data
                            a[2 * k] = Re(X[k]), 
                            a[2 * k + 1] = Im(X[k]), 0 <= k < n

                X[k] = sum_j=0^n-1 x[j]*exp(2*pi*i*j*k/n), 0<=k<n
                */
                for (int i = 0; i < n; i++)
                {
                    double data = timeDomainDatas[i];
                    a[2 * i] = data;
                    a[2 * i + 1] = 0;
                }

                int ipSize = (int)(2 + Math.Sqrt(n)) + 1;
                int[] ip = new int[ipSize];
                ip[0] = 0;
                // w[],ip[] are initialized if ip[0] == 0

                int wSize = n / 2;
                double[] w = new double[wSize];

                unsafe
                {
                    fixed (double* aP = &a[0])
                    fixed (int* ipP = &ip[0])
                    fixed (double* wP = &w[0])
                    {
                        IvyFEM.FFT.ImportedFunctions.OouraFFT_cdftsg(2 * n, isgn, aP, ipP, wP);
                    }
                }

                freqDomaindatas = new System.Numerics.Complex[dataCnt];
                for (int i = 0; i < n; i++)
                {
                    double real = a[2 * i];
                    double imag = a[2 * i + 1];
                    // フーリエ複素振幅
                    System.Numerics.Complex data = new System.Numerics.Complex(real, imag);
                    freqDomaindatas[i] = data;
                }
            }

            // 時間長さ
            double tl = dt * dataCnt;
            // 周波数刻み
            double df = 1.0 / tl;
            // 周波数アレイ
            freqs = new double[dataCnt];
            for (int i = 0; i < dataCnt; i++)
            {
                freqs[i] = i * df;
            }
        }
    }
}
