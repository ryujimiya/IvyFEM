using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class LineFE
    {
        protected double[] Calc1stN(double[] L)
        {
            double[] N = new double[2];

            // N = L
            L.CopyTo(N, 0);

            return N;
        }

        protected double[][] Calc1stNu(double[] L)
        {
            double[][] Nu = new double[2][];
            double[] a;
            double[] b;
            CalcTransMatrix(out a, out b);

            // dN/dx
            Nu[0] = b;

            return Nu;
        }

        protected double[,] Calc1stSNN()
        {
            double l = GetLineLength();
            double[,] sNN = new double[2, 2]
            {
                { 2.0, 1.0 },
                { 1.0, 2.0 }
            };
            for (int i = 0; i < sNN.GetLength(0); i++)
            {
                for (int j = 0; j < sNN.GetLength(1); j++)
                {
                    sNN[i, j] *= l / 6.0;
                }
            }
            return sNN;
        }

        protected double[,] Calc1stSNxNx()
        {
            double l = GetLineLength();
            double[,] sNxNx = new double[2, 2]
            {
                { 1.0, -1.0 },
                { -1.0, 1.0 }
            };
            for (int i = 0; i < sNxNx.GetLength(0); i++)
            {
                for (int j = 0; j < sNxNx.GetLength(1); j++)
                {
                    sNxNx[i, j] *= 1.0 / l;
                }
            }
            return sNxNx;
        }
    }
}
