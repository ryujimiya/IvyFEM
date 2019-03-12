using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class TriangleFE
    {
        public double[] Calc2ndN(double[] L)
        {
            double[] N = new double[6];

            for (int i = 0; i < 3; i++)
            {
                N[i] = L[i] * (2.0 * L[i] - 1.0);
                N[i + 3] = 4.0 * L[i] * L[(i + 1) % 3];
            }

            return N;
        }

        public double[][] Calc2ndNu(double[] L)
        {
            double[][] Nu = new double[2][];
            double[] a;
            double[] b;
            double[] c;
            CalcTransMatrix(out a, out b, out c);

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

        public double[] Calc2ndSN()
        {
            double A = GetArea();
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

        public double[,] Calc2ndSNN()
        {
            double A = GetArea();
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

        public double[,][,] Calc2ndSNuNv()
        {
            double[,][,] sNuNv = new double[6, 6][,];

            double A = GetArea();
            double[] a;
            double[] b;
            double[] c;
            CalcTransMatrix(out a, out b, out c);
            // dN/duのL1, L2, L3の係数, 定数項
            double[,] cx = new double[6, 4];
            double[,] cy = new double[6, 4];
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

        private void __CalcNuNv(double[,] sNuNv, double A, double[,]cu, double[,]cv)
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
    }
}
