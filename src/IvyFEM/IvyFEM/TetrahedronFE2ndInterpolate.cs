using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TetrahedronFE2ndInterpolate : IInterpolate
    {
        public TetrahedronFE Owner { get; set; }

        public TetrahedronFE2ndInterpolate()
        {

        }

        public TetrahedronFE2ndInterpolate(TetrahedronFE owner)
        {
            Owner = owner;
        }


        public uint GetNodeCount()
        {
            return 10;
        }

        public double[] GetNodeL(int nodeId)
        {
            double[][] nodeL = new double[10][]
            {
                new double[] { 1.0, 0.0, 0.0, 0.0},
                new double[] { 0.0, 1.0, 0.0, 0.0 },
                new double[] { 0.0, 0.0, 1.0, 0.0 },
                new double[] { 0.0, 0.0, 0.0, 1.0 },
                new double[] { 0.5, 0.5, 0.0, 0.0 },
                new double[] { 0.0, 0.5, 0.5, 0.0 },
                new double[] { 0.5, 0.0, 0.5, 0.0 },
                new double[] { 0.5, 0.0, 0.0, 0.5 },
                new double[] { 0.0, 0.5, 0.0, 0.5 },
                new double[] { 0.0, 0.0, 0.5, 0.5 }
            };
            return nodeL[nodeId];
        }

        public double[] CalcN(double[] L)
        {
            double[] N = new double[10];

            for (int i = 0; i < 4; i++)
            {
                N[i] = L[i] * (2.0 * L[i] - 1.0);
            }
            for (int i = 0; i < 3; i++)
            {
                N[i + 4] = 4.0 * L[i] * L[(i + 1) % 3];
                N[i + 7] = 4.0 * L[i] * L[3];
            }

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
            Nu[0] = new double[10];
            for (int i = 0; i < 4; i++)
            {
                Nu[0][i] = b[i] * (4.0 * L[i] - 1.0);
            }
            for (int i = 0; i < 3; i++)
            {
                Nu[0][i + 4] = 4.0 * (b[i] * L[(i + 1) % 3] + b[(i + 1) % 3] * L[i]);
                Nu[0][i + 7] = 4.0 * (b[i] * L[3] + b[3] * L[i]);
            }

            // dN/dy
            Nu[1] = new double[10];
            for (int i = 0; i < 4; i++)
            {
                Nu[1][i] = c[i] * (4.0 * L[i] - 1.0);
            }
            for (int i = 0; i < 3; i++)
            {
                Nu[1][i + 4] = 4.0 * (c[i] * L[(i + 1) % 3] + c[(i + 1) % 3] * L[i]);
                Nu[1][i + 7] = 4.0 * (c[i] * L[3] + c[3] * L[i]);
            }

            // dN/dz
            Nu[2] = new double[10];
            for (int i = 0; i < 4; i++)
            {
                Nu[2][i] = d[i] * (4.0 * L[i] - 1.0);
            }
            for (int i = 0; i < 3; i++)
            {
                Nu[2][i + 4] = 4.0 * (d[i] * L[(i + 1) % 3] + d[(i + 1) % 3] * L[i]);
                Nu[2][i + 7] = 4.0 * (d[i] * L[3] + d[3] * L[i]);
            }

            return Nu;
        }

        public double[,][] CalcNuv(double[] L)
        {
            double[,][] Nuv = new double[3, 3][];
            double[] a;
            double[] b;
            double[] c;
            double[] d;
            Owner.CalcTransMatrix(out a, out b, out c, out d);

            // d^2N/dx^2
            double[] Nxx = new double[10];
            Nuv[0, 0] = Nxx;
            {
                for (int i = 0; i < 4; i++)
                {
                    Nxx[i] = b[i] * 4.0 * b[i];
                }
                for (int i = 0; i < 3; i++)
                {
                    Nxx[i + 4] = 4.0 * (b[i] * b[(i + 1) % 3] + b[(i + 1) % 3] * b[i]);
                    Nxx[i + 7] = 4.0 * (b[i] * b[3] + b[3] * b[i]);
                }
            }
            // d^2N/dxdy
            double[] Nxy = new double[10];
            Nuv[0, 1] = Nxy;
            {
                for (int i = 0; i < 4; i++)
                {
                    Nxy[i] = b[i] * 4.0 * c[i];
                }
                for (int i = 0; i < 3; i++)
                {
                    Nxy[i + 4] = 4.0 * (b[i] * c[(i + 1) % 3] + b[(i + 1) % 3] * c[i]);
                    Nxy[i + 7] = 4.0 * (b[i] * c[3] + b[3] * c[i]);
                }
            }
            // d^2N/dxdz
            double[] Nxz = new double[10];
            Nuv[0, 2] = Nxz;
            {
                for (int i = 0; i < 4; i++)
                {
                    Nxz[i] = b[i] * 4.0 * d[i];
                }
                for (int i = 0; i < 3; i++)
                {
                    Nxz[i + 4] = 4.0 * (b[i] * d[(i + 1) % 3] + b[(i + 1) % 3] * d[i]);
                    Nxz[i + 7] = 4.0 * (b[i] * d[3] + b[3] * d[i]);
                }
            }
            // d^2N/dydx
            double[] Nyx = new double[10];
            Nuv[1, 0] = Nyx;
            {
                for (int i = 0; i < 4; i++)
                {
                    Nyx[i] = c[i] * 4.0 * b[i];
                }
                for (int i = 0; i < 3; i++)
                {
                    Nyx[i + 4] = 4.0 * (c[i] * b[(i + 1) % 3] + c[(i + 1) % 3] * b[i]);
                    Nyx[i + 7] = 4.0 * (c[i] * b[3] + c[3] * b[i]);
                }
            }
            // d^2N/dy^2
            double[] Nyy = new double[10];
            Nuv[1, 1] = Nyy;
            {
                for (int i = 0; i < 4; i++)
                {
                    Nyy[i] = c[i] * 4.0 * c[i];
                }
                for (int i = 0; i < 3; i++)
                {
                    Nyy[i + 4] = 4.0 * (c[i] * c[(i + 1) % 3] + c[(i + 1) % 3] * c[i]);
                    Nyy[i + 7] = 4.0 * (c[i] * c[3] + c[3] * c[i]);
                }
            }
            // d^2N/dydz
            double[] Nyz = new double[10];
            Nuv[1, 2] = Nyz;
            {
                for (int i = 0; i < 4; i++)
                {
                    Nyy[i] = c[i] * 4.0 * d[i];
                }
                for (int i = 0; i < 3; i++)
                {
                    Nyy[i + 4] = 4.0 * (c[i] * d[(i + 1) % 3] + c[(i + 1) % 3] * d[i]);
                    Nyy[i + 7] = 4.0 * (c[i] * d[3] + c[3] * d[i]);
                }
            }
            // d^2N/dzdx
            double[] Nzx = new double[10];
            Nuv[2, 0] = Nzx;
            {
                for (int i = 0; i < 4; i++)
                {
                    Nzx[i] = d[i] * 4.0 * b[i];
                }
                for (int i = 0; i < 3; i++)
                {
                    Nzx[i + 4] = 4.0 * (d[i] * b[(i + 1) % 3] + d[(i + 1) % 3] * b[i]);
                    Nzx[i + 7] = 4.0 * (d[i] * b[3] + d[3] * b[i]);
                }
            }
            // d^2N/dzdy
            double[] Nzy = new double[10];
            Nuv[2, 1] = Nzy;
            {
                for (int i = 0; i < 4; i++)
                {
                    Nzy[i] = d[i] * 4.0 * c[i];
                }
                for (int i = 0; i < 3; i++)
                {
                    Nzy[i + 4] = 4.0 * (d[i] * c[(i + 1) % 3] + d[(i + 1) % 3] * c[i]);
                    Nzy[i + 7] = 4.0 * (d[i] * c[3] + d[3] * c[i]);
                }
            }
            // d^2N/dz^2
            double[] Nzz = new double[10];
            Nuv[2, 2] = Nzz;
            {
                for (int i = 0; i < 4; i++)
                {
                    Nzz[i] = d[i] * 4.0 * d[i];
                }
                for (int i = 0; i < 3; i++)
                {
                    Nzz[i + 4] = 4.0 * (d[i] * d[(i + 1) % 3] + d[(i + 1) % 3] * d[i]);
                    Nzz[i + 7] = 4.0 * (d[i] * d[3] + d[3] * d[i]);
                }
            }

            return Nuv;
        }

        public double[] CalcSN()
        {
            double V = Owner.GetVolume();
            double[] sN = new double[10]
            {
                -(1.0 / 20.0) * V,
                -(1.0 / 20.0) * V,
                -(1.0 / 20.0) * V,
                -(1.0 / 20.0) * V,
                (1.0 /5.0) * V,
                (1.0 /5.0) * V,
                (1.0 /5.0) * V,
                (1.0 /5.0) * V,
                (1.0 /5.0) * V,
                (1.0 /5.0) * V
            };
            return sN;
        }

        public double[,] CalcSNN()
        {
            double V = Owner.GetVolume();
            double[,] sNN = new double[10, 10]
            {
                { 6.0, 1.0, 1.0, 1.0, -4.0, -6.0, -4.0, -4.0, -6.0, -6.0 },
                { 1.0, 6.0, 1.0, 1.0, -4.0, -4.0, -6.0, -6.0, -4.0, -6.0 },
                { 1.0, 1.0, 6.0, 1.0, -6.0, -4.0, -4.0, -6.0, -6.0, -4.0 },
                { 1.0, 1.0, 1.0, 6.0, -6.0, -6.0, -6.0, -4.0, -4.0, -4.0 },
                { -4.0, -4.0, -6.0, -6.0, 32.0, 16.0, 16.0, 16.0, 16.0, 8.0 },
                { -6.0, -4.0, -4.0, -6.0, 16.0, 32.0, 16.0, 8.0, 16.0, 16.0 },
                { -4.0, -6.0, -4.0, -6.0, 16.0, 16.0, 32.0, 16.0, 8.0, 16.0 },
                { -4.0, -6.0, -6.0, -4.0, 16.0, 8.0, 16.0, 32.0, 16.0, 16.0 },
                { -6.0, -4.0, -6.0, -4.0, 16.0, 16.0, 8.0, 16.0, 32.0, 16.0 },
                { -6.0, -6.0, -4.0, -4.0, 8.0, 16.0, 16.0, 16.0, 16.0, 32.0 }
            };
            for (int i = 0; i < sNN.GetLength(0); i++)
            {
                for (int j = 0; j < sNN.GetLength(1); j++)
                {
                    sNN[i, j] *= V / 420.0;
                }
            }
            return sNN;
        }

        public double[,][,] CalcSNuNv()
        {
            double[,][,] sNuNv = new double[3, 3][,];

            double V = Owner.GetVolume();
            double[] a;
            double[] b;
            double[] c;
            double[] d;
            Owner.CalcTransMatrix(out a, out b, out c, out d);

            // dN/duのL1, L2, L3, L4の係数, 定数項
            double[,] cx;
            double[,] cy;
            double[,] cz;
            CalcDNduCoef(a, b, c, d, out cx, out cy, out cz);

            // sNxNx
            sNuNv[0, 0] = new double[10, 10];
            __CalcNuNv(sNuNv[0, 0], V, cx, cx);
            // sNxNy
            sNuNv[0, 1] = new double[10, 10];
            __CalcNuNv(sNuNv[0, 1], V, cx, cy);
            // sNxNz
            sNuNv[0, 2] = new double[10, 10];
            __CalcNuNv(sNuNv[0, 2], V, cx, cz);
            // sNyNx
            sNuNv[1, 0] = new double[10, 10];
            __CalcNuNv(sNuNv[1, 0], V, cy, cx);
            // sNyNy
            sNuNv[1, 1] = new double[10, 10];
            __CalcNuNv(sNuNv[1, 1], V, cy, cy);
            // sNyNz
            sNuNv[1, 2] = new double[10, 10];
            __CalcNuNv(sNuNv[1, 2], V, cy, cy);
            // sNzNx
            sNuNv[2, 0] = new double[10, 10];
            __CalcNuNv(sNuNv[2, 0], V, cy, cx);
            // sNzNy
            sNuNv[2, 1] = new double[10, 10];
            __CalcNuNv(sNuNv[2, 1], V, cy, cy);
            // sNzNz
            sNuNv[2, 2] = new double[10, 10];
            __CalcNuNv(sNuNv[2, 2], V, cy, cy);
            return sNuNv;
        }

        protected void CalcDNduCoef(
            double[] a, double[] b, double[] c, double[] d,
            out double[,] cx, out double[,] cy, out double[,] cz)
        {
            // dN/duのL1, L2, L3, L4の係数, 定数項
            cx = new double[10, 5];
            cy = new double[10, 5];
            cz = new double[10, 5];

            // dN/dx
            for (int i = 0; i < 4; i++)
            {
                //Nu[0][i] = b[i] * (4.0 * L[i] - 1.0);
                cx[i, i] = b[i] * 4.0;
                cx[i, 4] = -b[i];
            }
            for (int i = 0; i < 3; i++)
            {
                //Nu[0][i + 4] = 4.0 * (b[i] * L[(i + 1) % 3] + b[(i + 1) % 3] * L[i]);
                cx[i + 4, (i + 1) % 3] = 4.0 * b[i];
                cx[i + 4, i] = 4.0 * b[(i + 1) % 3];
                //Nu[0][i + 7] = 4.0 * (b[i] * L[3] + b[3] * L[i]);
                cx[i + 7, 3] = 4.0 * b[i];
                cx[i + 7, i] = 4.0 * b[3];
            }

            // dN/dy
            for (int i = 0; i < 4; i++)
            {
                //Nu[1][i] = c[i] * (4.0 * L[i] - 1.0);
                cy[i, i] = c[i] * 4.0;
                cy[i, 4] = -c[i];
            }
            for (int i = 0; i < 3; i++)
            {
                //Nu[1][i + 4] = 4.0 * (c[i] * L[(i + 1) % 3] + c[(i + 1) % 3] * L[i]);
                cy[i + 4, (i + 1) % 3] = 4.0 * c[i];
                cy[i + 4, i] = 4.0 * c[(i + 1) % 3];
                //Nu[1][i + 7] = 4.0 * (c[i] * L[3] + c[3] * L[i]);
                cy[i + 7, 3] = 4.0 * c[i];
                cy[i + 7, i] = 4.0 * c[3];
            }

            // dN/dz
            for (int i = 0; i < 4; i++)
            {
                //Nu[2][i] = d[i] * (4.0 * L[i] - 1.0);
                cz[i, i] = d[i] * 4.0;
                cz[i, 4] = -d[i];
            }
            for (int i = 0; i < 3; i++)
            {
                //Nu[2][i + 4] = 4.0 * (d[i] * L[(i + 1) % 3] + d[(i + 1) % 3] * L[i]);
                cz[i + 4, (i + 1) % 3] = 4.0 * d[i];
                cz[i + 4, i] = 4.0 * d[(i + 1) % 3];
                //Nu[2][i + 7] = 4.0 * (d[i] * L[3] + d[3] * L[i]);
                cz[i + 7, 3] = 4.0 * d[i];
                cz[i + 7, i] = 4.0 * d[3];
            }
        }

        protected void __CalcNuNv(double[,] sNuNv, double V, double[,] cu, double[,] cv)
        {
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    sNuNv[i, j] =
                        (V / 10.0) * (
                        cu[i, 0] * cv[j, 0] + cu[i, 1] * cv[j, 1] + cu[i, 2] * cv[j, 2] + cu[i, 3] * cv[j, 3]) +
                        (V / 20.0) * (
                        cu[i, 0] * cv[j, 1] + cu[i, 0] * cv[j, 2] + cu[i, 0] * cv[j, 3] +
                        cu[i, 1] * cv[j, 0] + cu[i, 1] * cv[j, 2] + cu[i, 1] * cv[j, 3] +
                        cu[i, 2] * cv[j, 0] + cu[i, 2] * cv[j, 1] + cu[i, 2] * cv[j, 3] +
                        cu[i, 3] * cv[j, 0] + cu[i, 3] * cv[j, 1] + cu[i, 3] * cv[j, 2]) +
                        (V / 4.0) * (
                        cu[i, 0] * cv[j, 4] + cu[i, 1] * cv[j, 4] + cu[i, 2] * cv[j, 4] + cu[i, 3] * cv[j, 4] +
                        cu[i, 4] * cv[j, 0] + cu[i, 4] * cv[j, 1] + cu[i, 4] * cv[j, 2] + cu[i, 4] * cv[j, 3]) +
                        V * cu[i, 4] * cv[j, 4];
                }
            }
        }

        public double[][,] CalcSNuN()
        {
            double V = Owner.GetVolume();
            double[] a;
            double[] b;
            double[] c;
            double[] d;
            Owner.CalcTransMatrix(out a, out b, out c, out d);

            // dN/duのL1, L2, L3, L4の係数, 定数項
            double[,] cx;
            double[,] cy;
            double[,] cz;
            CalcDNduCoef(a, b, c, d, out cx, out cy, out cz);

            double[][,] sNuN = new double[3][,];

            // sNxN
            sNuN[0] = new double[10, 10];
            __CalcNuN(sNuN[0], V, cx);

            // sNyN
            sNuN[1] = new double[10, 10];
            __CalcNuN(sNuN[1], V, cy);

            // sNzN
            sNuN[2] = new double[10, 10];
            __CalcNuN(sNuN[2], V, cz);

            return sNuN;
        }

        protected void __CalcNuN(double[,] sNuN, double V, double[,] cu)
        {
            for (int i = 0; i < 10; i++)
            {
                sNuN[i, 0] =
                    2.0 * cu[i, 0] * (V / 20.0) + 2.0 * cu[i, 1] * (V / 60.0) +
                    2.0 * cu[i, 2] * (V / 60.0) + 2.0 * cu[i, 3] * (V / 60.0) +
                    2.0 * cu[i, 4] * (V / 10.0) -
                    cu[i, 0] * (V / 10.0) - cu[i, 1] * (V / 20.0) -
                    cu[i, 2] * (V / 20.0) - cu[i, 3] * (V / 20.0) -
                    cu[i, 4] * (V / 4.0);
                sNuN[i, 1] =
                    2.0 * cu[i, 0] * (V / 60.0) + 2.0 * cu[i, 1] * (V / 20.0) +
                    2.0 * cu[i, 2] * (V / 60.0) + 2.0 * cu[i, 3] * (V / 60.0) +
                    2.0 * cu[i, 4] * (V / 10.0) -
                    cu[i, 0] * (V / 20.0) - cu[i, 1] * (V / 10.0) -
                    cu[i, 2] * (V / 20.0) - cu[i, 3] * (V / 20.0) -
                    cu[i, 4] * (V / 4.0);
                sNuN[i, 2] =
                    2.0 * cu[i, 0] * (V / 60.0) + 2.0 * cu[i, 1] * (V / 60.0) +
                    2.0 * cu[i, 2] * (V / 20.0) + 2.0 * cu[i, 3] * (V / 60.0) +
                    2.0 * cu[i, 4] * (V / 10.0) -
                    cu[i, 0] * (V / 20.0) - cu[i, 1] * (V / 20.0) -
                    cu[i, 2] * (V / 10.0) - cu[i, 3] * (V / 20.0) -
                    cu[i, 4] * (V / 4.0);
                sNuN[i, 3] =
                    2.0 * cu[i, 0] * (V / 60.0) + 2.0 * cu[i, 1] * (V / 60.0) +
                    2.0 * cu[i, 2] * (V / 60.0) + 2.0 * cu[i, 3] * (V / 20.0) +
                    2.0 * cu[i, 4] * (V / 10.0) -
                    cu[i, 0] * (V / 20.0) - cu[i, 1] * (V / 20.0) -
                    cu[i, 2] * (V / 20.0) - cu[i, 3] * (V / 10.0) -
                    cu[i, 4] * (V / 4.0);
                sNuN[i, 4] =
                    4.0 * cu[i, 0] * (V / 60.0) + 4.0 * cu[i, 1] * (V / 60.0) +
                    4.0 * cu[i, 2] * (V / 120.0) + 4.0 * cu[i, 3] * (V / 120.0) +
                    4.0 * cu[i, 4] * (V / 20.0);
                sNuN[i, 5] =
                    4.0 * cu[i, 0] * (V / 120.0) + 4.0 * cu[i, 1] * (V / 60.0) +
                    4.0 * cu[i, 2] * (V / 60.0) + 4.0 * cu[i, 3] * (V / 120.0) +
                    4.0 * cu[i, 4] * (V / 20.0);
                sNuN[i, 6] =
                    4.0 * cu[i, 0] * (V / 60.0) + 4.0 * cu[i, 1] * (V / 120.0) +
                    4.0 * cu[i, 2] * (V / 60.0) + 4.0 * cu[i, 3] * (V / 120.0) +
                    4.0 * cu[i, 4] * (V / 20.0);
                sNuN[i, 7] =
                    4.0 * cu[i, 0] * (V / 60.0) + 4.0 * cu[i, 1] * (V / 120.0) +
                    4.0 * cu[i, 2] * (V / 120.0) + 4.0 * cu[i, 3] * (V / 60.0) +
                    4.0 * cu[i, 4] * (V / 20.0);
                sNuN[i, 8] =
                    4.0 * cu[i, 0] * (V / 120.0) + 4.0 * cu[i, 1] * (V / 60.0) +
                    4.0 * cu[i, 2] * (V / 120.0) + 4.0 * cu[i, 3] * (V / 60.0) +
                    4.0 * cu[i, 4] * (V / 20.0);
                sNuN[i, 9] =
                    4.0 * cu[i, 0] * (V / 120.0) + 4.0 * cu[i, 1] * (V / 120.0) +
                    4.0 * cu[i, 2] * (V / 60.0) + 4.0 * cu[i, 3] * (V / 60.0) +
                    4.0 * cu[i, 4] * (V / 20.0);
            }
        }
    }
}
