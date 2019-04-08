using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class LineFE
    {
        protected double[] Calc2ndN(double[] L)
        {
            double[] N = new double[3];

            for (int i = 0; i < 2; i++)
            {
                N[i] = L[i] * (2.0 * L[i] - 1.0);
            }
            N[2] = 4.0 * L[0] * L[1];

            return N;
        }

        protected double[][] Calc2ndNu(double[] L)
        {
            double[][] Nu = new double[1][];
            double[] a;
            double[] b;
            CalcTransMatrix(out a, out b);

            // dN/dx
            Nu[0] = new double[3];
            for (int i = 0; i < 2; i++)
            {
                Nu[0][i] = b[i] * (4.0 * L[i] - 1.0);
            }
            Nu[0][2] = 4.0 * (b[0] * L[1] + b[1] * L[0]);

            return Nu;
        }

        protected double[,] Calc2ndSNN()
        {
            double l = GetLineLength();
            double[,] sNN = new double[3, 3]
            {
                { 4.0, -1.0, 2.0 },
                { -1.0, 4.0, 2.0 },
                { 2.0, 2.0, 16.0 }
            };
            for (int i = 0; i < sNN.GetLength(0); i++)
            {
                for (int j = 0; j < sNN.GetLength(1); j++)
                {
                    sNN[i, j] *= l / 30.0;
                }
            }
            return sNN;
        }

        protected double[,] Calc2ndSNxNx()
        {
            double l = GetLineLength();
            double[,] sNxNx = new double[3, 3] 
            {
                { 7.0, 1.0, -8.0 },
                { 1.0, 7.0, -8.0 },
                { -8.0, -8.0, 16.0 }
            };
            for (int i = 0; i < sNxNx.GetLength(0); i++)
            {
                for (int j = 0; j < sNxNx.GetLength(1); j++)
                {
                    sNxNx[i, j] *= 1.0 / (3.0 * l);
                }
            }
            return sNxNx;
        }
    }
}
