using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class LineFE1stInterpolate : IInterpolate
    {
        public LineFE Owner { get; set; }

        public  LineFE1stInterpolate()
        {

        }

        public LineFE1stInterpolate(LineFE owner)
        {
            Owner = owner;
        }


        public uint GetNodeCount()
        {
            return 2;
        }

        public double[] GetNodeL(int nodeId)
        {
            double[][] nodeL = new double[2][]
            {
                new double[] { 1.0, 0.0 },
                new double[] { 0.0, 1.0 }
            };
            return nodeL[nodeId];
        }

        public double[] CalcN(double[] L)
        {
            double[] N = new double[2];

            // N = L
            L.CopyTo(N, 0);

            return N;
        }

        public double[][] CalcNu(double[] L)
        {
            double[][] Nu = new double[1][];
            double[] a;
            double[] b;
            Owner.CalcTransMatrix(out a, out b);

            // dN/dx
            Nu[0] = b;

            return Nu;
        }

        public double[,][] CalcNuv(double[] L)
        {
            throw new NotImplementedException();
        }

        public double[] CalcSN()
        {
            double l = Owner.GetLineLength();
            double[] sN = new double[2]
            {
                l / 2.0,
                l / 2.0
            };
            return sN;
        }

        public double[,] CalcSNN()
        {
            double l = Owner.GetLineLength();
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

        public double[][,] CalcSNuN()
        {
            double[,] sNxN = new double[2, 2]
            {
                { -1.0, -1.0 },
                { 1.0, 1.0 }
            };
            for (int i = 0; i < sNxN.GetLength(0); i++)
            {
                for (int j = 0; j < sNxN.GetLength(1); j++)
                {
                    sNxN[i, j] *= 1.0 / 2.0;
                }
            }
            double[][,] sNuN = { sNxN };
            return sNuN;
        }
    }
}
