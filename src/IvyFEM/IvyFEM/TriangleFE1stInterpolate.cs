using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TriangleFE1stInterpolate : IInterpolate
    {
        public TriangleFE Owner { get; set; }

        public TriangleFE1stInterpolate()
        {

        }

        public TriangleFE1stInterpolate(TriangleFE owner)
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
                new double[] { 1.0, 0.0, 0.0 },
                new double[] { 0.0, 1.0, 0.0 },
                new double[] { 0.0, 0.0, 1.0 }
            };
            return nodeL[nodeId];
        }

        public double[] CalcN(double[] L)
        {
            double[] N = new double[3];

            // N = L
            L.CopyTo(N, 0);

            return N;
        }

        public double[][] CalcNu(double[] L)
        {
            double[][] Nu = new double[2][];
            double[] a;
            double[] b;
            double[] c;
            Owner.CalcTransMatrix(out a, out b, out c);

            // dN/dx
            Nu[0] = b;

            // dN/dy
            Nu[1] = c;
            return Nu;
        }

        public double[,][] CalcNuv(double[] L)
        {
            double[,][] Nuv = new double[2, 2][];
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    double[] value = new double[3];
                    Nuv[i, j] = value;
                    value[0] = 0;
                    value[1] = 0;
                    value[2] = 0;
                }
            }
            return Nuv;
        }

        public double[] CalcSN()
        {
            double A = Owner.GetArea();
            double[] sN = new double[3]
            {
                A / 3.0,
                A / 3.0,
                A / 3.0
            };
            return sN;
        }

        public double[,] CalcSNN()
        {
            double A = Owner.GetArea();
            double[,] sNN = new double[3, 3]
            {
                { 2.0, 1.0, 1.0 },
                { 1.0, 2.0, 1.0 },
                { 1.0, 1.0, 2.0 }
            };
            for (int i = 0; i < sNN.GetLength(0); i++)
            {
                for (int j = 0; j < sNN.GetLength(1); j++)
                {
                    sNN[i, j] *= A / 12.0;
                }
            }
            return sNN;
        }

        public double[,][,] CalcSNuNv()
        {
            double A = Owner.GetArea();
            double[] a;
            double[] b;
            double[] c;
            Owner.CalcTransMatrix(out a, out b, out c);

            double[,][,] sNuNv = new double[2, 2][,];

            // sNxNx
            sNuNv[0, 0] = new double[3, 3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    sNuNv[0, 0][i, j] = A * b[i] * b[j];
                }
            }

            // sNxNy
            sNuNv[0, 1] = new double[3, 3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    sNuNv[0, 1][i, j] = A * b[i] * c[j];
                }
            }

            // sNyNx
            sNuNv[1, 0] = new double[3, 3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    sNuNv[1, 0][i, j] = A * c[i] * b[j];
                }
            }

            // sNyNy
            sNuNv[1, 1] = new double[3, 3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    sNuNv[1, 1][i, j] = A * c[i] * c[j];
                }
            }
            return sNuNv;
        }

        public double[][,] CalcSNuN()
        {
            double A = Owner.GetArea();
            double[] a;
            double[] b;
            double[] c;
            Owner.CalcTransMatrix(out a, out b, out c);

            double[][,] sNuN = new double[2][,];

            // sNxN
            sNuN[0] = new double[3, 3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    sNuN[0][i, j] = (A / 3.0) * b[i];
                }
            }

            // sNyN
            sNuN[1] = new double[3, 3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    sNuN[1][i, j] = (A / 3.0) * c[i];
                }
            }

            return sNuN;
        }
    }
}
