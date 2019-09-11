using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TriangleFE2ndInterpolate : IInterpolate
    {
        public TriangleFE Owner { get; set; }

        public TriangleFE2ndInterpolate()
        {

        }

        public TriangleFE2ndInterpolate(TriangleFE owner)
        {
            Owner = owner;
        }


        public uint GetNodeCount()
        {
            return 6;
        }

        public double[] GetNodeL(int nodeId)
        {
            double[][] nodeL = new double[6][]
            {
                new double[] { 1.0, 0.0, 0.0 },
                new double[] { 0.0, 1.0, 0.0 },
                new double[] { 0.0, 0.0, 1.0 },
                new double[] { 1.0 / 2.0, 1.0 / 2.0, 0.0 },
                new double[] { 0.0, 1.0 / 2.0, 1.0 / 2.0 },
                new double[] { 1.0 / 2.0, 0.0, 1.0 / 2.0 }
            };
            return nodeL[nodeId];
        }

        public double[] CalcN(double[] L)
        {
            double[] N = new double[6];

            for (int i = 0; i < 3; i++)
            {
                N[i] = L[i] * (2.0 * L[i] - 1.0);
                N[i + 3] = 4.0 * L[i] * L[(i + 1) % 3];
            }

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
            Nu[0] = new double[6];
            for (int i = 0; i < 3; i++)
            {
                Nu[0][i] = b[i] * (4.0 * L[i] - 1.0);
                Nu[0][i + 3] = 4.0 * (b[i] * L[(i + 1) % 3] + b[(i + 1) % 3] * L[i]);
            }

            // dN/dy
            Nu[1] = new double[6];
            for (int i = 0; i < 3; i++)
            {
                Nu[1][i] = c[i] * (4.0 * L[i] - 1.0);
                Nu[1][i + 3] = 4.0 * (c[i] * L[(i + 1) % 3] + c[(i + 1) % 3] * L[i]);
            }

            return Nu;
        }

        public double[,][] CalcNuv(double[] L)
        {
            double[,][] Nuv = new double[2, 2][];
            double[] a;
            double[] b;
            double[] c;
            Owner.CalcTransMatrix(out a, out b, out c);

            // d^2N/dx^2
            double[] Nxx = new double[6];
            Nuv[0, 0] = Nxx;
            {
                for (int i = 0; i < 3; i++)
                {
                    Nxx[i] = b[i] * 4.0 * b[i];
                    Nxx[i + 3] = 4.0 * (b[i] * b[(i + 1) % 3] + b[(i + 1) % 3] * b[i]);
                }
            }
            // d^2N/dxdy
            double[] Nxy = new double[6];
            Nuv[0, 1] = Nxy;
            {
                for (int i = 0; i < 3; i++)
                {
                    Nxy[i] = c[i] * 4.0 * b[i];
                    Nxy[i + 3] = 4.0 * (c[i] * b[(i + 1) % 3] + c[(i + 1) % 3] * b[i]);
                }
            }
            // d^2N/dydx
            double[] Nyx = new double[6];
            Nuv[1, 0] = Nyx;
            {
                for (int i = 0; i < 3; i++)
                {
                    Nyx[i] = b[i] * 4.0 * c[i];
                    Nyx[i + 3] = 4.0 * (b[i] * c[(i + 1) % 3] + b[(i + 1) % 3] * c[i]);
                }
            }
            // d^2N/dy^2
            double[] Nyy = new double[6];
            Nuv[1, 1] = Nyy;
            {
                for (int i = 0; i < 3; i++)
                {
                    Nyy[i] = c[i] * 4.0 * c[i];
                    Nyy[i + 3] = 4.0 * (c[i] * c[(i + 1) % 3] + c[(i + 1) % 3] * c[i]);
                }
            }
            return Nuv;
        }

        public double[] CalcSN()
        {
            double A = Owner.GetArea();
            double[] sN = new double[6]
            {
                0.0,
                0.0,
                0.0,
                A / 3.0,
                A / 3.0,
                A / 3.0
            };
            return sN;
        }

        public double[,] CalcSNN()
        {
            double A = Owner.GetArea();
            double[,] sNN = new double[6, 6]
            {
                { 6.0, -1.0, -1.0, 0.0, -4.0, 0.0 },
                { -1.0, 6.0, -1.0, 0.0, 0.0, -4.0 },
                { -1.0, -1.0, 6.0, -4.0, 0.0, 0.0 },
                { 0.0, 0.0, -4.0, 32.0, 16.0, 16.0 },
                { -4.0, 0.0, 0.0, 16.0, 32.0, 16.0 },
                { 0.0, -4.0, 0.0, 16.0, 16.0, 32.0 }
            };
            for (int i = 0; i < sNN.GetLength(0); i++)
            {
                for (int j = 0; j < sNN.GetLength(1); j++)
                {
                    sNN[i, j] *= A / 180.0;
                }
            }
            return sNN;
        }

        public double[,][,] CalcSNuNv()
        {
            double[,][,] sNuNv = new double[6, 6][,];

            double A = Owner.GetArea();
            double[] a;
            double[] b;
            double[] c;
            Owner.CalcTransMatrix(out a, out b, out c);

            // dN/duのL1, L2, L3の係数, 定数項
            double[,] cx;
            double[,] cy;
            CalcDNduCoef(a, b, c, out cx, out cy);

            // sNxNx
            sNuNv[0, 0] = new double[6, 6];
            __CalcNuNv(sNuNv[0, 0], A, cx, cx);
            // sNxNy
            sNuNv[0, 1] = new double[6, 6];
            __CalcNuNv(sNuNv[0, 1], A, cx, cy);
            // sNyNx
            sNuNv[1, 0] = new double[6, 6];
            __CalcNuNv(sNuNv[1, 0], A, cy, cx);
            // sNyNy
            sNuNv[1, 1] = new double[6, 6];
            __CalcNuNv(sNuNv[1, 1], A, cy, cy);
            return sNuNv;
        }

        protected void CalcDNduCoef(
            double[]a, double[]b, double[]c,
            out double[,] cx, out double[,] cy)
        {
            // dN/duのL1, L2, L3の係数, 定数項
            cx = new double[6, 4];
            cy = new double[6, 4];
            // dN/dx
            for (int i = 0; i < 3; i++)
            {
                //Nu[0][i] = b[i] * (4.0 * L[i] - 1.0);
                cx[i, i] = b[i] * 4.0;
                cx[i, 3] = -b[i];
                //Nu[0][i + 3] = 4.0 * (b[i] * L[(i + 1) % 3] + b[(i + 1) % 3] * L[i]);
                cx[i + 3, (i + 1) % 3] = 4.0 * b[i];
                cx[i + 3, i] = 4.0 * b[(i + 1) % 3];
            }
            // dN/dy
            for (int i = 0; i < 3; i++)
            {
                //Nu[1][i] = c[i] * (4.0 * L[i] - 1.0);
                cy[i, i] = c[i] * 4.0;
                cy[i, 3] = -c[i];
                //Nu[1][i + 3] = 4.0 * (c[i] * L[(i + 1) % 3] + c[(i + 1) % 3] * L[i]);
                cy[i + 3, (i + 1) % 3] = 4.0 * c[i];
                cy[i + 3, i] = 4.0 * c[(i + 1) % 3];
            }

        }

        protected void __CalcNuNv(double[,] sNuNv, double A, double[,]cu, double[,]cv)
        {
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    sNuNv[i, j] =
                        (A / 6.0) * (cu[i, 0] * cv[j, 0] + cu[i, 1] * cv[j, 1] + cu[i, 2] * cv[j, 2]) +
                        (A / 12.0) * (cu[i, 0] * cv[j, 1] + cu[i, 0] * cv[j, 2] +
                        cu[i, 1] * cv[j, 0] + cu[i, 1] * cv[j, 2] +
                        cu[i, 2] * cv[j, 0] + cu[i, 2] * cv[j, 1]) +
                        (A / 3.0) * (cu[i, 0] * cv[j, 3] + cu[i, 1] * cv[j, 3] + cu[i, 2] * cv[j, 3] +
                        cu[i, 3] * cv[j, 0] + cu[i, 3] * cv[j, 1] + cu[i, 3] * cv[j, 2]) +
                        A * cu[i, 3] * cv[j, 3];
                }
            }
        }

        public double[][,] CalcSNuN()
        {
            double[][,] sNuN = new double[6][,];

            double A = Owner.GetArea();
            double[] a;
            double[] b;
            double[] c;
            Owner.CalcTransMatrix(out a, out b, out c);

            // dN/duのL1, L2, L3の係数, 定数項
            double[,] cx;
            double[,] cy;
            CalcDNduCoef(a, b, c, out cx, out cy);

            // sNxN
            double[,] sNxN = new double[6, 6];
            sNuN[0] = sNxN;
            {
                for (int i = 0; i < 6; i++)
                {
                    sNxN[i, 0] = (A / 30.0) * cx[i, 0] - (A / 60.0) * (cx[i, 1] + cx[i, 2]);
                    sNxN[i, 1] = (A / 30.0) * cx[i, 1] - (A / 60.0) * (cx[i, 0] + cx[i, 2]);
                    sNxN[i, 2] = (A / 30.0) * cx[i, 2] - (A / 60.0) * (cx[i, 0] + cx[i, 1]);
                    sNxN[i, 3] = (A / 15.0) * (2.0 * cx[i, 0] + 2.0 * cx[i, 1] + cx[i, 2] + 5.0 * cx[i, 3]);
                    sNxN[i, 4] = (A / 15.0) * (cx[i, 0] + 2.0 * cx[i, 1] + 2.0 * cx[i, 2] + 5.0 * cx[i, 3]);
                    sNxN[i, 5] = (A / 15.0) * (2.0 * cx[i, 0] + cx[i, 1] + 2.0 * cx[i, 2] + 5.0 * cx[i, 3]);
                }

            }
            // sNyN
            double[,] sNyN = new double[6, 6];
            sNuN[1] = sNyN;
            {
                for (int i = 0; i < 6; i++)
                {
                    sNyN[i, 0] = (A / 30.0) * cy[i, 0] - (A / 60.0) * (cy[i, 1] + cy[i, 2]);
                    sNyN[i, 1] = (A / 30.0) * cy[i, 1] - (A / 60.0) * (cy[i, 0] + cy[i, 2]);
                    sNyN[i, 2] = (A / 30.0) * cy[i, 2] - (A / 60.0) * (cy[i, 0] + cy[i, 1]);
                    sNyN[i, 3] = (A / 15.0) * (2.0 * cy[i, 0] + 2.0 * cy[i, 1] + cy[i, 2] + 5.0 * cy[i, 3]);
                    sNyN[i, 4] = (A / 15.0) * (cy[i, 0] + 2.0 * cy[i, 1] + 2.0 * cy[i, 2] + 5.0 * cy[i, 3]);
                    sNyN[i, 5] = (A / 15.0) * (2.0 * cy[i, 0] + cy[i, 1] + 2.0 * cy[i, 2] + 5.0 * cy[i, 3]);
                }

            }
            return sNuN;
        }
    }
}
