using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class LineFE2ndInterpolate : IInterpolate
    {
        public LineFE Owner { get; set; }

        public LineFE2ndInterpolate()
        {

        }

        public LineFE2ndInterpolate(LineFE owner)
        {
            Owner = owner;
        }

        public uint GetNodeCount()
        {
            return 3;
        }

        public double[] GetNodeL(int nodeId)
        {
            double[][] nodeL = new double[3][]
            {
                new double[] { 1.0, 0.0 },
                new double[] { 0.0, 1.0 },
                new double[] { 1.0 / 2.0, 1.0 / 2.0 }
            };
            return nodeL[nodeId];
        }

        public double[] CalcN(double[] L)
        {
            double[] N = new double[3];

            for (int i = 0; i < 2; i++)
            {
                N[i] = L[i] * (2.0 * L[i] - 1.0);
            }
            N[2] = 4.0 * L[0] * L[1];

            return N;
        }

        public double[][] CalcNu(double[] L)
        {
            double[][] Nu = new double[1][];
            double[] a;
            double[] b;
            Owner.CalcTransMatrix(out a, out b);

            // dN/dx
            Nu[0] = new double[3];
            for (int i = 0; i < 2; i++)
            {
                Nu[0][i] = b[i] * (4.0 * L[i] - 1.0);
            }
            Nu[0][2] = 4.0 * (b[0] * L[1] + b[1] * L[0]);

            return Nu;
        }

        public double[,][] CalcNuv(double[] L)
        {
            throw new NotImplementedException();
        }

        public double[] CalcSN()
        {
            double l = Owner.GetLineLength();
            double[] sN = new double[3]
            {
                (1.0 / 6.0) * l,
                (1.0 / 6.0) * l,
                (4.0 / 6.0) * l
            };
            return sN;
        }

        public double[,] CalcSNN()
        {
            double l = Owner.GetLineLength();
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

        public double[,][,] CalcSNuNv()
        {
            double[,][,] sNuNv = new double[1, 1][,];
            double[,] sNxNx = CalcSNxNx();
            sNuNv[0, 0] = sNxNx;
            return sNuNv;
        }

        protected double[,] CalcSNxNx()
        {
            double l = Owner.GetLineLength();
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

        public double[][,] CalcSNuN()
        {
            double[,] sNxN = new double[3, 3]
            {
                { -3.0, 1.0, -4.0 },
                { -1.0, 3.0, 4.0 },
                { 4.0, -4.0, 0.0 }
            };
            for (int i = 0; i < sNxN.GetLength(0); i++)
            {
                for (int j = 0; j < sNxN.GetLength(1); j++)
                {
                    sNxN[i, j] *= 1.0 / 6.0;
                }
            }
            double[][,] sNuN = { sNxN };
            return sNuN;
        }
    }
}
