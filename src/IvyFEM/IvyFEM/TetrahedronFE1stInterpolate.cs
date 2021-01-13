using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TetrahedronFE1stInterpolate : IInterpolate
    {
        public TetrahedronFE Owner { get; set; }

        public TetrahedronFE1stInterpolate()
        {

        }

        public TetrahedronFE1stInterpolate(TetrahedronFE owner)
        {
            Owner = owner;
        }


        public uint GetNodeCount()
        {
            return 4;
        }

        public double[] GetNodeL(int nodeId)
        {
            double[][] nodeL = new double[4][]
            {
                new double[] { 1.0, 0.0, 0.0, 0.0},
                new double[] { 0.0, 1.0, 0.0, 0.0 },
                new double[] { 0.0, 0.0, 1.0, 0.0 },
                new double[] { 0.0, 0.0, 0.0, 1.0 }
            };
            return nodeL[nodeId];
        }

        public double[] CalcN(double[] L)
        {
            double[] N = new double[4];

            // N = L
            L.CopyTo(N, 0);

            return N;
        }

        public double[][] CalcNu(double[] L)
        {
            double[][] Nu = new double[3][];
            double[] a;
            double[] b;
            double[] c;
            double[] d;
            Owner.CalcTransMatrix(out a, out b, out c, out d);

            // dN/dx
            Nu[0] = b;

            // dN/dy
            Nu[1] = c;

            // dN/dz
            Nu[2] = d;

            return Nu;
        }

        public double[,][] CalcNuv(double[] L)
        {
            double[,][] Nuv = new double[3, 3][];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    double[] value = new double[4];
                    Nuv[i, j] = value;
                    value[0] = 0;
                    value[1] = 0;
                    value[2] = 0;
                    value[3] = 0;
                }
            }
            return Nuv;
        }

        public double[] CalcSN()
        {
            double V = Owner.GetVolume();
            double[] sN = new double[4]
            {
                V / 4.0,
                V / 4.0,
                V / 4.0,
                V / 4.0
            };
            return sN;
        }

        public double[,] CalcSNN()
        {
            double V = Owner.GetVolume();
            double[,] sNN = new double[4, 4]
            {
                { 2.0, 1.0, 1.0, 1.0 },
                { 1.0, 2.0, 1.0, 1.0 },
                { 1.0, 1.0, 2.0, 1.0 },
                { 1.0, 1.0, 1.0, 2.0 }
            };
            for (int i = 0; i < sNN.GetLength(0); i++)
            {
                for (int j = 0; j < sNN.GetLength(1); j++)
                {
                    sNN[i, j] *= V / 20.0;
                }
            }
            return sNN;
        }

        public double[,][,] CalcSNuNv()
        {
            double V = Owner.GetVolume();
            double[] a;
            double[] b;
            double[] c;
            double[] d;
            Owner.CalcTransMatrix(out a, out b, out c, out d);

            double[,][,] sNuNv = new double[3, 3][,];

            // sNxNx
            sNuNv[0, 0] = new double[4, 4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    sNuNv[0, 0][i, j] = V * b[i] * b[j];
                }
            }

            // sNxNy
            sNuNv[0, 1] = new double[4, 4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    sNuNv[0, 1][i, j] = V * b[i] * c[j];
                }
            }

            // sNxNz
            sNuNv[0, 2] = new double[4, 4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    sNuNv[0, 2][i, j] = V * b[i] * d[j];
                }
            }

            // sNyNx
            sNuNv[1, 0] = new double[4, 4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    sNuNv[1, 0][i, j] = V * c[i] * b[j];
                }
            }

            // sNyNy
            sNuNv[1, 1] = new double[4, 4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    sNuNv[1, 1][i, j] = V * c[i] * c[j];
                }
            }

            // sNyNz
            sNuNv[1, 2] = new double[4, 4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    sNuNv[1, 2][i, j] = V * c[i] * d[j];
                }
            }

            // sNzNx
            sNuNv[2, 0] = new double[4, 4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    sNuNv[2, 0][i, j] = V * d[i] * b[j];
                }
            }

            // sNzNy
            sNuNv[2, 1] = new double[4, 4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    sNuNv[2, 1][i, j] = V * d[i] * c[j];
                }
            }

            // sNzNz
            sNuNv[2, 2] = new double[4, 4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    sNuNv[2, 2][i, j] = V * d[i] * d[j];
                }
            }

            return sNuNv;
        }

        public double[][,] CalcSNuN()
        {
            double V = Owner.GetVolume();
            double[] a;
            double[] b;
            double[] c;
            double[] d;
            Owner.CalcTransMatrix(out a, out b, out c, out d);

            double[][,] sNuN = new double[3][,];

            // sNxN
            sNuN[0] = new double[4, 4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    sNuN[0][i, j] = (V / 4.0) * b[i];
                }
            }

            // sNyN
            sNuN[1] = new double[4, 4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    sNuN[1][i, j] = (V / 4.0) * c[i];
                }
            }

            // sNzN
            sNuN[2] = new double[4, 4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    sNuN[2][i, j] = (V / 4.0) * d[i];
                }
            }

            return sNuN;
        }
    }
}
