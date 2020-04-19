using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DFEMUtils
    {
        private static void CalcFieldConsistentTLFrameUAtLocalX(double x, double l0, double[] ul,
            out double[] u, out double[] up, out double[] upp)
        {
            System.Diagnostics.Debug.Assert(ul.Length == 7);
            double ux1 = ul[0];
            double uy1 = ul[1];
            double tz1 = ul[2];
            double ux2 = ul[3];
            double uy2 = ul[4];
            double tz2 = ul[5];
            double ux3 = ul[6];
            double tantz1 = Math.Tan(tz1);
            double tantz2 = Math.Tan(tz2);

            // u
            double ux = ux1 - (x / l0) * (3.0 * ux1 + ux2 - 4.0 * ux3) +
                ((x * x) / (l0 * l0)) * 2.0 * (ux1 + ux2 - 2.0 * ux3);
            double uy = uy1 + (x / l0) * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * tantz1 -
                ((x * x) / (l0 * l0)) * (
                3.0 * (uy1 - uy2) +
                2.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * tantz1 +
                (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * tantz2
                ) +
                ((x * x * x) / (l0 * l0 * l0)) * (
                2.0 * (uy1 - uy2) +
                (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * tantz1 +
                (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * tantz2
                );
            u = new double[] { ux, uy };

            // u'
            double upx = -(1.0 / l0) * (3.0 * ux1 + ux2 - 4.0 * ux3) +
                (x / (l0 * l0)) * 4.0 * (ux1 + ux2 - 2.0 * ux3);
            double upy = (1.0 / l0) * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * tantz1 -
                (x / (l0 * l0)) * (
                6.0 * (uy1 - uy2) +
                4.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * tantz1 +
                2.0 * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * tantz2
                ) +
                ((x * x) / (l0 * l0 * l0)) * (
                6.0 * (uy1 - uy2) +
                3.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * tantz1 +
                3.0 * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * tantz2
                );
            up = new double[] { upx, upy };

            // u''
            double uppx = (1.0 / (l0 * l0)) * 4.0 * (ux1 + ux2 - 2.0 * ux3);
            double uppy = -(1.0 / (l0 * l0)) * (
                6.0 * (uy1 - uy2) +
                4.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * tantz1 +
                2.0 * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * tantz2
                ) +
                (x / (l0 * l0 * l0)) * (
                12.0 * (uy1 - uy2) +
                6.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * tantz1 +
                6.0 * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * tantz2
                );
            upp = new double[] { uppx, uppy };
        }

        private static void CalcFieldConsistentTLFrameNsAndNpsAndNppsAtLocalX(double x, double l0, double[] ul,
            out double[][] Ns, out double[][] Nps, out double[][] Npps)
        {
            System.Diagnostics.Debug.Assert(ul.Length == 7);
            double ux1 = ul[0];
            double uy1 = ul[1];
            double tz1 = ul[2];
            double ux2 = ul[3];
            double uy2 = ul[4];
            double tz2 = ul[5];
            double ux3 = ul[6];
            double tantz1 = Math.Tan(tz1);
            double tantz2 = Math.Tan(tz2);
            double costz1 = Math.Cos(tz1);
            double costz2 = Math.Cos(tz2);

            //------------------------------------------------
            double Nx1 = 1.0 - 3.0 * (x / l0) + 2.0 * ((x * x) / (l0 * l0));
            double Nx2 = 0.0;
            double Nx3 = 0.0;
            double Nx4 = -(x / l0) + 2.0 * ((x * x) / (l0 * l0));
            double Nx5 = 0.0;
            double Nx6 = 0.0;
            double Nx7 = 4.0 * (x / l0) - 4.0 * ((x * x) / (l0 * l0));
            double[] Nx = { Nx1, Nx2, Nx3, Nx4, Nx5, Nx6, Nx7 };
            //
            double Ny1 = -(x / l0) * 3.0 * tantz1 +
                ((x * x) / (l0 * l0)) * (6.0 * tantz1 - tantz2) -
                ((x * x * x) / (l0 * l0 * l0)) * (3.0 * tantz1 - tantz2);
            double Ny2 = 1.0 - ((x * x) / (l0 * l0)) * 3.0 + ((x * x * x) / (l0 * l0 * l0)) * 2.0;
            double Ny3 = (x / l0) * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (1.0 / (costz1 * costz1)) -
                ((x * x) / (l0 * l0)) * 2.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (1.0 / (costz1 * costz1)) +
                ((x * x * x) / (l0 * l0 * l0)) * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (1.0 / (costz1 * costz1));
            double Ny4 = -(x / l0) * tantz1 +
                ((x * x) / (l0 * l0)) * (2.0 * tantz1 - 3.0 * tantz2) -
                ((x * x * x) / (l0 * l0 * l0)) * (tantz1 - 3.0 * tantz2);
            double Ny5 = ((x * x) / (l0 * l0)) * 3.0 - ((x * x * x) / (l0 * l0 * l0)) * 2.0;
            double Ny6 = -((x * x) / (l0 * l0)) * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (1.0 / (costz2 * costz2)) +
                ((x * x * x) / (l0 * l0 * l0)) * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (1.0 / (costz2 * costz2));
            double Ny7 = (x / l0) * 4.0 * tantz1 -
                ((x * x) / (l0 * l0)) * (8.0 * tantz1 - 4.0 * tantz2) +
                ((x * x * x) / (l0 * l0 * l0)) * (4.0 * tantz1 - 4.0 * tantz2);
            double[] Ny = { Ny1, Ny2, Ny3, Ny4, Ny5, Ny6, Ny7 };
            //
            Ns = new double[2][] { Nx, Ny };

            //------------------------------------------------
            double Npx1 = -3.0 * (1.0 / l0) + 4.0 * (x / (l0 * l0));
            double Npx2 = 0.0;
            double Npx3 = 0.0;
            double Npx4 = -(1.0 / l0) + 4.0 * (x / (l0 * l0));
            double Npx5 = 0.0;
            double Npx6 = 0.0;
            double Npx7 = 4.0 * (1.0 / l0) - 8.0 * (x / (l0 * l0));
            double[] Npx = { Npx1, Npx2, Npx3, Npx4, Npx5, Npx6, Npx7 };
            //
            double Npy1 = -(1.0 / l0) * 3.0 * tantz1 +
                (x / (l0 * l0)) * (12.0 * tantz1 - 2.0 * tantz2) -
                ((x * x) / (l0 * l0 * l0)) * (9.0 * tantz1 - 3.0 * tantz2);
            double Npy2 = -(x / (l0 * l0)) * 6.0 + ((x * x) / (l0 * l0 * l0)) * 6.0;
            double Npy3 = (1.0 / l0) * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (1.0 / (costz1 * costz1)) -
                (x / (l0 * l0)) * 4.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (1.0 / (costz1 * costz1)) +
                ((x * x) / (l0 * l0 * l0)) * 3.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (1.0 / (costz1 * costz1));
            double Npy4 = -(1.0 / l0) * tantz1 +
                (x / (l0 * l0)) * (4.0 * tantz1 - 6.0 * tantz2) -
                ((x * x) / (l0 * l0 * l0)) * (3.0 * tantz1 - 9.0 * tantz2);
            double Npy5 = (x / (l0 * l0)) * 6.0 - ((x * x) / (l0 * l0 * l0)) * 6.0;
            double Npy6 = -(x / (l0 * l0)) * 2.0 * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (1.0 / (costz2 * costz2)) +
                ((x * x) / (l0 * l0 * l0)) * 3.0 * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (1.0 / (costz2 * costz2));
            double Npy7 = (1.0 / l0) * 4.0 * tantz1 -
                (x / (l0 * l0)) * (16.0 * tantz1 - 8.0 * tantz2) +
                ((x * x) / (l0 * l0 * l0)) * (12.0 * tantz1 - 12.0 * tantz2);
            double[] Npy = { Npy1, Npy2, Npy3, Npy4, Npy5, Npy6, Npy7 };
            //
            Nps = new double[2][] { Npx, Npy };

            //------------------------------------------------
            double Nppx1 = 4.0 / (l0 * l0);
            double Nppx2 = 0.0;
            double Nppx3 = 0.0;
            double Nppx4 = 4.0 / (l0 * l0);
            double Nppx5 = 0.0;
            double Nppx6 = 0.0;
            double Nppx7 = -8.0 / (l0 * l0);
            double[] Nppx = { Nppx1, Nppx2, Nppx3, Nppx4, Nppx5, Nppx6, Nppx7 };
            //
            double Nppy1 = (1.0 / (l0 * l0)) * (12.0 * tantz1 - 2.0 * tantz2) -
                (x / (l0 * l0 * l0)) * (18.0 * tantz1 - 6.0 * tantz2);
            double Nppy2 = -6.0 / (l0 * l0) + (x / (l0 * l0 * l0)) * 12.0;
            double Nppy3 = -(1.0 / (l0 * l0)) * 4.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (1.0 / (costz1 * costz1)) +
                (x / (l0 * l0 * l0)) * 6.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (1.0 / (costz1 * costz1));
            double Nppy4 = (1.0 / (l0 * l0)) * (4.0 * tantz1 - 6.0 * tantz2) -
                (x / (l0 * l0 * l0)) * (6.0 * tantz1 - 18.0 * tantz2);
            double Nppy5 = 6.0 / (l0 * l0) - (x / (l0 * l0 * l0)) * 12.0;
            double Nppy6 = -(1.0 / (l0 * l0)) * 2.0 * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (1.0 / (costz2 * costz2)) +
                (x / (l0 * l0 * l0)) * 6.0 * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (1.0 / (costz2 * costz2));
            double Nppy7 = -(1.0 / (l0 * l0)) * (16.0 * tantz1 - 8.0 * tantz2) +
                (x / (l0 * l0 * l0)) * (24.0 * tantz1 - 24.0 * tantz2);
            double[] Nppy = { Nppy1, Nppy2, Nppy3, Nppy4, Nppy5, Nppy6, Nppy7 };
            //
            Npps = new double[2][] { Nppx, Nppy };
        }

        private static void CalcFieldConsistentTLFrameNtAndNptAtLocalX(
            double[] u, double[] up, double[] upp,
            double[][] Ns, double[][] Nps, double[][] Npps,
            out double[] Nt, out double[] Npt)
        {
            double upx = up[0];
            double upy = up[1];
            double uppx = upp[0];
            double uppy = upp[1];
            double[] Npx = Nps[0];
            double[] Npy = Nps[1];
            double[] Nppx = Npps[0];
            double[] Nppy = Npps[1];

            double c1 = (1.0 + upx) * (1.0 + upx) + upy * upy;
            Nt = new double[7];
            {
                double d1 = upy / c1;
                double d2 = (1.0 + upx) / c1;
                for (int i = 0; i < 7; i++)
                {
                    Nt[i] = -d1 * Npx[i] + d2 * Npy[i];
                }
            }
            // 
            Npt = new double[7];
            {
                double d1 = (
                    -2.0 * (1.0 + upx) * uppx * upy +
                    (1.0 + upx) * (1.0 + upx) * uppy -
                    upy * upy * uppy
                    ) / (c1 * c1);
                double d2 = upy / c1;
                double d3 = (-1.0 * (1.0 + upx) * (1.0 + upx) * uppx -
                    2.0 * upy * uppy * (1.0 + upx) +
                    upy * upy * uppx) / (c1 * c1);
                double d4 = (1.0 + upx) / c1;
                for (int i = 0; i < 7; i++)
                {
                    Npt[i] = -d1 * Npx[i] - d2 * Nppx[i] +
                        d3 * Npy[i] + d4 * Nppy[i];
                }
            }
        }

        private static void CalcFieldConsistentTLFrameNpussAndNppussAtLocalX(double x, double l0, double[] ul,
            out double[][][] Npuss, out double[][][] Nppuss)
        {
            System.Diagnostics.Debug.Assert(ul.Length == 7);
            double ux1 = ul[0];
            double uy1 = ul[1];
            double tz1 = ul[2];
            double ux2 = ul[3];
            double uy2 = ul[4];
            double tz2 = ul[5];
            double ux3 = ul[6];
            double tantz1 = Math.Tan(tz1);
            double tantz2 = Math.Tan(tz2);
            double costz1 = Math.Cos(tz1);
            double costz2 = Math.Cos(tz2);
            double sintz1 = Math.Sin(tz1);
            double sintz2 = Math.Sin(tz2);

            //---------------------------------------------------------------
            double[][] Npxus = new double[7][];
            for (int i = 0; i < 7; i++)
            {
                Npxus[i] = new double[7];
            }

            //
            double Npy1ux1 = 0.0;
            double Npy1uy1 = 0.0;
            double Npy1tz1 = -(1.0 / l0) * 3.0 * (1.0 / (costz1 * costz1)) +
                (x / (l0 * l0)) * 12.0 * (1.0 / (costz1 * costz1)) -
                ((x * x) / (l0 * l0 * l0)) * 9.0 * (1.0 / (costz1 * costz1));
            double Npy1ux2 = 0.0;
            double Npy1uy2 = 0.0;
            double Npy1tz2 = -(x / (l0 * l0)) * (2.0 / (costz2 * costz2)) +
                ((x * x) / (l0 * l0 * l0)) * (3.0 / (costz2 * costz2));
            double Npy1ux3 = 0.0;
            double[] Npy1u = { Npy1ux1, Npy1uy1, Npy1tz1, Npy1ux2, Npy1uy2, Npy1tz2, Npy1ux3 };
            //
            double[] Npy2u = new double[7];
            //
            double Npy3ux1 = -(1.0 / l0) * (3.0 / (costz1 * costz1)) +
                (x / (l0 * l0)) * (12.0 / (costz1 * costz1)) -
                ((x * x) / (l0 * l0 * l0)) * (9.0 / (costz1 * costz1));
            double Npy3uy1 = 0.0;
            double Npy3tz1 =
                (1.0 / l0) * 2.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (sintz1 / (costz1 * costz1 * costz1)) -
                (x / (l0 * l0)) * 8.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (sintz1 / (costz1 * costz1 * costz1)) +
                ((x * x) / (l0 * l0 * l0)) * 6.0 *
                (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (sintz1 / (costz1 * costz1 * costz1));
            double Npy3ux2 = -(1.0 / l0) * (1.0 / (costz1 * costz1)) +
                (x / (l0 * l0)) * (4.0 / (costz1 * costz1)) -
                ((x * x) / (l0 * l0 * l0)) * (3.0 / (costz1 * costz1));
            double Npy3uy2 = 0.0;
            double Npy3tz2 = 0.0;
            double Npy3ux3 = (1.0 / l0) * (4.0 / (costz1 * costz1)) -
                (x / (l0 * l0)) * (16.0 / (costz1 * costz1)) +
                ((x * x) / (l0 * l0 * l0)) * (12.0 / (costz1 * costz1));
            double[] Npy3u = { Npy3ux1, Npy3uy1, Npy3tz1, Npy3ux2, Npy3uy2, Npy3tz2, Npy3ux3 };
            //
            double Npy4ux1 = 0.0;
            double Npy4uy1 = 0.0;
            double Npy4tz1 = -(1.0 / l0) * (1.0 / (costz1 * costz1)) +
                (x / (l0 * l0)) * 4.0 * (1.0 / (costz1 * costz1)) -
                ((x * x) / (l0 * l0 * l0)) * 3.0 * (1.0 / (costz1 * costz1));
            double Npy4ux2 = 0.0;
            double Npy4uy2 = 0.0;
            double Npy4tz2 = -(x / (l0 * l0)) * (6.0 / (costz2 * costz2)) +
                ((x * x) / (l0 * l0 * l0)) * (9.0 / (costz2 * costz2));
            double Npy4ux3 = 0.0;
            double[] Npy4u = { Npy4ux1, Npy4uy1, Npy4tz1, Npy4ux2, Npy4uy2, Npy4tz2, Npy4ux3 };
            //
            double[] Npy5u = new double[7];
            //
            double Npy6ux1 = -(x / (l0 * l0)) * (2.0 / (costz2 * costz2)) +
                ((x * x) / (l0 * l0 * l0)) * (3.0 / (costz2 * costz2));
            double Npy6uy1 = 0.0;
            double Npy6tz1 = 0.0;
            double Npy6ux2 = -(x / (l0 * l0)) * (6.0 / (costz2 * costz2)) +
                ((x * x) / (l0 * l0 * l0)) * (9.0 / (costz2 * costz2));
            double Npy6uy2 = 0.0;
            double Npy6tz2 =
                -(x / (l0 * l0)) * 4.0 * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (sintz2 / (costz2 * costz2 * costz2)) +
                ((x * x) / (l0 * l0 * l0)) * 6.0 *
                (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (sintz2 / (costz2 * costz2 * costz2));
            double Npy6ux3 = (x / (l0 * l0)) * (8.0 / (costz2 * costz2)) -
                ((x * x) / (l0 * l0 * l0)) * (12.0 / (costz2 * costz2));
            double[] Npy6u = { Npy6ux1, Npy6uy1, Npy6tz1, Npy6ux2, Npy6uy2, Npy6tz2, Npy6ux3 };
            //
            double Npy7ux1 = 0.0;
            double Npy7uy1 = 0.0;
            double Npy7tz1 = (1.0 / l0) * (4.0 / (costz1 * costz1)) -
                (x / (l0 * l0)) * (16.0 / (costz1 * costz1)) +
                ((x * x) / (l0 * l0 * l0)) * (12.0 / (costz1 * costz1));
            double Npy7ux2 = 0.0;
            double Npy7uy2 = 0.0;
            double Npy7tz2 = (x / (l0 * l0)) * (8.0 / (costz2 * costz2)) -
                ((x * x) / (l0 * l0 * l0)) * (12.0 / (costz2 * costz2));
            double Npy7ux3 = 0.0;
            double[] Npy7u = { Npy7ux1, Npy7uy1, Npy7tz1, Npy7ux2, Npy7uy2, Npy7tz2, Npy7ux3 };

            double[][] Npyus = { Npy1u, Npy2u, Npy3u, Npy4u, Npy5u, Npy6u, Npy7u };

            Npuss = new double[2][][];
            Npuss[0] = Npxus;
            Npuss[1] = Npyus;

            //---------------------------------------------------------------
            double[][] Nppxus = new double[7][];
            for (int i = 0; i < 7; i++)
            {
                Nppxus[i] = new double[7];
            }

            //
            double Nppy1ux1 = 0.0;
            double Nppy1uy1 = 0.0;
            double Nppy1tz1 = (1.0 / (l0 * l0)) * (12.0 / (costz1 * costz1)) -
                (x / (l0 * l0 * l0)) * (18.0 / (costz1 * costz1));
            double Nppy1ux2 = 0.0;
            double Nppy1uy2 = 0.0;
            double Nppy1tz2 = -(1.0 / (l0 * l0)) * (2.0 / (costz2 * costz2)) +
                (x / (l0 * l0 * l0)) * (6.0 / (costz2 * costz2));
            double Nppy1ux3 = 0.0;
            double[] Nppy1u = { Nppy1ux1, Nppy1uy1, Nppy1tz1, Nppy1ux2, Nppy1uy2, Nppy1tz2, Nppy1ux3 };
            //
            double[] Nppy2u = new double[7];
            //
            double Nppy3ux1 = (1.0 / (l0 * l0)) * (12.0 / (costz1 * costz1)) -
                (x / (l0 * l0 * l0)) * (18.0 / (costz1 * costz1));
            double Nppy3uy1 = 0.0;
            double Nppy3tz1 =
                -(1.0 / (l0 * l0)) * 8.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (sintz1 / (costz1 * costz1 * costz1)) +
                (x / (l0 * l0 * l0)) * 12.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (sintz1 / (costz1 * costz1 * costz1));
            double Nppy3ux2 = (1.0 / (l0 * l0)) * (8.0 / (costz1 * costz1)) -
                (x / (l0 * l0 * l0)) * (6.0 / (costz1 * costz1));
            double Nppy3uy2 = 0.0;
            double Nppy3tz2 = 0.0;
            double Nppy3ux3 = -(1.0 / (l0 * l0)) * (16.0 / (costz1 * costz1)) +
                (x / (l0 * l0 * l0)) * (24.0 / (costz1 * costz1));
            double[] Nppy3u = { Nppy3ux1, Nppy3uy1, Nppy3tz1, Nppy3ux2, Nppy3uy2, Nppy3tz2, Nppy3ux3 };
            //
            double Nppy4ux1 = 0.0;
            double Nppy4uy1 = 0.0;
            double Nppy4tz1 = (1.0 / (l0 * l0)) * (4.0 / (costz1 * costz1)) -
                (x / (l0 * l0 * l0)) * (6.0 / (costz1 * costz1));
            double Nppy4ux2 = 0.0;
            double Nppy4uy2 = 0.0;
            double Nppy4tz2 = -(1.0 / (l0 * l0)) * (6.0 / (costz2 * costz2)) +
                (x / (l0 * l0 * l0)) * (18.0 / (costz2 * costz2));
            double Nppy4ux3 = 0.0;
            double[] Nppy4u = { Nppy4ux1, Nppy4uy1, Nppy4tz1, Nppy4ux2, Nppy4uy2, Nppy4tz2, Nppy4ux3 };
            //
            double[] Nppy5u = new double[7];
            //
            double Nppy6ux1 = -(1.0 / (l0 * l0)) * (2.0 / (costz2 * costz2)) +
                (x / (l0 * l0 * l0)) * (6.0 / (costz2 * costz2));
            double Nppy6uy1 = 0.0;
            double Nppy6tz1 = 0.0;
            double Nppy6ux2 = -(1.0 / (l0 * l0)) * (6.0 / (costz2 * costz2)) +
                (x / (l0 * l0 * l0)) * (18.0 / (costz2 * costz2));
            double Nppy6uy2 = 0.0;
            double Nppy6tz2 =
                -(1.0 / (l0 * l0)) * 4.0 * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (sintz2 / (costz2 * costz2 * costz2)) +
                (x / (l0 * l0 * l0)) * 12.0 * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (sintz2 / (costz2 * costz2 * costz2));
            double Nppy6ux3 = -(1.0 / (l0 * l0)) * (8.0 / (costz2 * costz2)) -
                (x / (l0 * l0 * l0)) * (24.0 / (costz2 * costz2));
            double[] Nppy6u = { Nppy6ux1, Nppy6uy1, Nppy6tz1, Nppy6ux2, Nppy6uy2, Nppy6tz2, Nppy6ux3 };
            //
            double Nppy7ux1 = 0.0;
            double Nppy7uy1 = 0.0;
            double Nppy7tz1 = -(1.0 / (l0 * l0)) * (16.0 / (costz1 * costz1)) +
                (x / (l0 * l0 * l0)) * (24.0 / (costz1 * costz1));
            double Nppy7ux2 = 0.0;
            double Nppy7uy2 = 0.0;
            double Nppy7tz2 = (1.0 / (l0 * l0)) * (8.0 / (costz2 * costz2)) -
                (x / (l0 * l0 * l0)) * (24.0 / (costz2 * costz2));
            double Nppy7ux3 = 0.0;
            double[] Nppy7u = { Nppy7ux1, Nppy7uy1, Nppy7tz1, Nppy7ux2, Nppy7uy2, Nppy7tz2, Nppy7ux3 };

            double[][] Nppyus = { Nppy1u, Nppy2u, Nppy3u, Nppy4u, Nppy5u, Nppy6u, Nppy7u };

            Nppuss = new double[2][][];
            Nppuss[0] = Nppxus;
            Nppuss[1] = Nppyus;
        }

        private static void CalcFieldConsistentTLFrameUpusAndUppusAtLocalX(double x, double l0, double[] ul,
            out double[][] upus, out double[][] uppus)
        {
            System.Diagnostics.Debug.Assert(ul.Length == 7);
            double ux1 = ul[0];
            double uy1 = ul[1];
            double tz1 = ul[2];
            double ux2 = ul[3];
            double uy2 = ul[4];
            double tz2 = ul[5];
            double ux3 = ul[6];
            double tantz1 = Math.Tan(tz1);
            double tantz2 = Math.Tan(tz2);
            double costz1 = Math.Cos(tz1);
            double costz2 = Math.Cos(tz2);

            //------------------------------------------------------------------
            double upxux1 = -(1.0 / l0) * 3.0 + (x / (l0 * l0)) * 4.0;
            double upxuy1 = 0.0;
            double upxtz1 = 0.0;
            double upxux2 = -(1.0 / l0) + (x / (l0 * l0)) * 4.0;
            double upxuy2 = 0.0;
            double upxtz2 = 0.0;
            double upxux3 = (1.0 / l0) * 4.0 - (x / (l0 * l0)) * 8.0;
            double[] upxu = { upxux1, upxuy1, upxtz1, upxux2, upxuy2, upxtz2, upxux3 };
            //
            double upyux1 = -(1.0 / l0) * 3.0 * tantz1 +
                (x / (l0 * l0)) * (12.0 * tantz1 - 2.0 * tantz2) -
                ((x * x) / (l0 * l0 * l0)) * (9.0 * tantz1 - 3.0 * tantz2);
            double upyuy1 = -(x / (l0 * l0)) * 6.0 + ((x * x) / (l0 * l0 * l0)) * 6.0;
            double upytz1 = (1.0 / l0) * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (1.0 / (costz1 * costz1)) -
                (x / (l0 * l0)) * 4.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (1.0 / (costz1 * costz1)) +
                ((x * x) / (l0 * l0 * l0)) * 3.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (1.0 / (costz1 * costz1));
            double upyux2 = -(1.0 / l0) * tantz1 +
                (x / (l0 * l0)) * (4.0 * tantz1 - 6.0 * tantz2) -
                ((x * x) / (l0 * l0 * l0)) * (3.0 * tantz1 - 9.0 * tantz2);
            double upyuy2 = (x / (l0 * l0)) * 6.0 - ((x * x) / (l0 * l0 * l0)) * 6.0;
            double upytz2 = -(x / (l0 * l0)) * 2.0 * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (1.0 / (costz2 * costz2)) +
                ((x * x) / (l0 * l0 * l0)) * 3.0 * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (1.0 / (costz2 * costz2));
            double upyux3 = (1.0 / l0) * 4.0 * tantz1 -
                (x / (l0 * l0)) * (16.0 * tantz1 - 8.0 * tantz2) +
                ((x * x) / (l0 * l0 * l0)) * (12.0 * tantz1 - 12.0 * tantz2);
            double[] upyu = { upyux1, upyuy1, upytz1, upyux2, upyuy2, upytz2, upyux3 };

            upus = new double[2][];
            upus[0] = upxu;
            upus[1] = upyu;

            //------------------------------------------------------------------
            double uppxux1 = 4.0 / (l0 * l0);
            double uppxuy1 = 0.0;
            double uppxtz1 = 0.0;
            double uppxux2 = 4.0 / (l0 * l0);
            double uppxuy2 = 0.0;
            double uppxtz2 = 0.0;
            double uppxux3 = -8.0 / (l0 * l0);
            double[] uppxu = { uppxux1, uppxuy1, uppxtz1, uppxux2, uppxuy2, uppxtz2, uppxux3 };
            //
            double uppyux1 = (1.0 / (l0 * l0)) * (12.0 * tantz1 - 2.0 * tantz2) -
                (x / (l0 * l0 * l0)) * (18.0 * tantz1 - 6.0 * tantz2);
            double uppyuy1 = -(1.0 / (l0 * l0)) * 6.0 + (x / (l0 * l0 * l0)) * 12.0;
            double uppytz1 = -(1.0 / (l0 * l0)) * 4.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (1.0 / (costz1 * costz1)) +
                (x / (l0 * l0 * l0)) * 6.0 * (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (1.0 / (costz1 * costz1));
            double uppyux2 = (1.0 / (l0 * l0)) * (4.0 * tantz1 - 6.0 * tantz2) -
                (x / (l0 * l0 * l0)) * (6.0 * tantz1 - 18.0 * tantz2);
            double uppyuy2 = (1.0 / (l0 * l0)) * 6.0 - (x / (l0 * l0 * l0)) * 12.0;
            double uppytz2 = -(1.0 / (l0 * l0)) * 2.0 * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (1.0 / (costz2 * costz2)) +
                (x / (l0 * l0 * l0)) * 6.0 * (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (1.0 / (costz2 * costz2));
            double uppyux3 = -(1.0 / (l0 * l0)) * (16.0 * tantz1 - 8.0 * tantz2) +
                (x / (l0 * l0 * l0)) * (24.0 * tantz1 - 24.0 * tantz2);
            double[] uppyu = { uppyux1, uppyuy1, uppytz1, uppyux2, uppyuy2, uppytz2, uppyux3 };

            uppus = new double[2][];
            uppus[0] = uppxu;
            uppus[1] = uppyu;
        }

        private static void CalcFieldConsistentTLFrameNptusAtLocalX(
            double[] u, double[] up, double[] upp,
            double[][] Nps, double[][] Npps,
            double[][][] Npuss, double[][][] Nppuss,
            double[][] upus, double[][] uppus,
            out double[][] Nptus)
        {
            double upx = up[0];
            double upy = up[1];
            double uppx = upp[0];
            double uppy = upp[1];
            double[] Npx = Nps[0];
            double[] Npy = Nps[1];
            double[] Nppx = Npps[0];
            double[] Nppy = Npps[1];
            double[][] Npxus = Npuss[0];
            double[][] Npyus = Npuss[1];
            double[][] Nppxus = Nppuss[0];
            double[][] Nppyus = Nppuss[1];
            double[] upxu = upus[0];
            double[] upyu = upus[1];
            double[] uppxu = uppus[0];
            double[] uppyu = uppus[1];

            Nptus = new double[7][];

            double c1 = (1.0 + upx) * (1.0 + upx) + upy * upy;
            double d2 = (1.0 / (c1 * c1)) * (
                -2.0 * (1.0 + upx) * uppx * upy +
                (1.0 + upx) * (1.0 + upx) * uppy -
                upy * upy * uppy
                );
            double d4 = (1.0 / c1) * upy;
            double d6 = (1.0 / (c1 * c1)) * (
                -(1.0 + upx) * (1.0 + upx) * uppx -
                2.0 * upy * uppy * (1.0 + upx) +
                upy * upy * uppx
                );
            double d8 = (1.0 / c1) * (1.0 + upx);
            for (int row = 0; row < 7; row++)
            {
                double[] Nptu = new double[7];
                Nptus[row] = Nptu;
                for (int col = 0; col < 7; col++)
                {
                    double d1 = -2.0 * (1.0 / (c1 * c1 * c1)) * (
                        2.0 * (1.0 + upx) * upxu[col] + 2.0 * upy * upyu[col]
                        ) * (
                        -2.0 * (1.0 + upx) * uppx * upy +
                        (1.0 + upx) * (1.0 + upx) * uppy -
                        upy * upy * uppy
                        ) +
                        (1.0 / (c1 * c1)) * (
                        -2.0 * upxu[col] * uppx * upy -
                        2.0 * (1.0 + upx) * uppxu[col] * upy -
                        2.0 * (1.0 + upx) * uppx * upyu[col] +
                        2.0 * (1.0 + upx) * upxu[col] * uppy +
                        (1.0 + upx) * (1.0 + upx) * uppyu[col] -
                        2.0 * upy * upyu[col] * uppy -
                        upy * upy * uppyu[col]
                        );
                    double d3 = -(1.0 / (c1 * c1)) * (
                        2.0 * (1.0 + upx) * upxu[col] + 2.0 * upy * upyu[col]
                        ) * upy +
                        (1.0 / c1) * upyu[col];
                    double d5 = -2.0 * (1.0 / (c1 * c1 * c1)) * (
                        2.0 * (1.0 + upx) * upxu[col] + 2.0 + upy * upyu[col]
                        ) * (
                        -(1.0 + upx) * (1.0 + upx) * uppx -
                        2.0 * upy * uppy * (1.0 + upx) +
                        upy * upy * uppx
                        ) +
                        (1.0 / (c1 * c1)) * (
                        -2.0 * (1.0 + upx) * upxu[col] * uppx -
                        (1.0 + upx) * (1.0 + upx) * uppxu[col] -
                        2.0 * upyu[col] * uppy * (1.0 + upx) -
                        2.0 * upy * uppyu[col] * (1.0 + upx) -
                        2.0 * upy * uppy * upxu[col] +
                        2.0 * upy * upyu[col] * uppx +
                        upy * upy * uppxu[col]
                        );
                    double d7 = -(1.0 / (c1 * c1)) * (
                        2.0 * (1.0 + upx) * upxu[col] + 2.0 * upy * upyu[col]
                        ) * (1.0 + upx) +
                        (1.0 / c1) * upxu[col];

                    Nptu[col] = -d1 * Npx[row] - d2 * Npxus[row][col] -
                        d3 * Nppx[row] - d4 * Nppxus[row][col] +
                        d5 * Npy[row] + d6 * Npyus[row][col] +
                        d7 * Nppy[row] + d8 * Nppyus[row][col];                        
                }
            }
        }

        public static void CalcFieldConsistentTLFrameKl(
            double E, double Ae, double I,
            double l0, double[] ul,
            out double[] fl, out IvyFEM.Lapack.DoubleMatrix kl)
        {
            fl = new double[7];
            kl = new IvyFEM.Lapack.DoubleMatrix(7, 7);

            IntegrationPoints ip = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);
            System.Diagnostics.Debug.Assert(ip.Ls.Length == 5);
            for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
            {
                double[] L = ip.Ls[ipPt];
                double lineLen = l0;
                double weight = ip.Weights[ipPt];
                double detJWeight = (lineLen / 2.0) * weight;

                double x = l0 * L[1];

                double[] u;
                double[] up;
                double[] upp;
                CalcFieldConsistentTLFrameUAtLocalX(x, l0, ul, out u, out up, out upp);
                double ux = u[0];
                double uy = u[1];
                double upx = up[0];
                double upy = up[1];
                double uppx = upp[0];
                double uppy = upp[1];

                double[][] Ns;
                double[][] Nps;
                double[][] Npps;
                CalcFieldConsistentTLFrameNsAndNpsAndNppsAtLocalX(x, l0, ul, out Ns, out Nps, out Npps);
                double[] Nx = Ns[0];
                double[] Ny = Ns[1];
                double[] Npx = Nps[0];
                double[] Npy = Nps[1];
                double[] Nppx = Npps[0];
                double[] Nppy = Npps[1];

                double[] Nt;
                double[] Npt;
                CalcFieldConsistentTLFrameNtAndNptAtLocalX(u, up, upp, Ns, Nps, Npps, out Nt, out Npt);

                double[][][] Npuss;
                double[][][] Nppuss;
                CalcFieldConsistentTLFrameNpussAndNppussAtLocalX(x, l0, ul, out Npuss, out Nppuss);
                double[][] Npxus = Npuss[0];
                double[][] Npyus = Npuss[1];
                double[][] Nppxus = Nppuss[0];
                double[][] Nppyus = Nppuss[1];

                double[][] upus;
                double[][] uppus;
                CalcFieldConsistentTLFrameUpusAndUppusAtLocalX(x, l0, ul, out upus, out uppus);
                double[] upxu = upus[0];
                double[] upyu = upus[1];
                double[] uppxu = uppus[0];
                double[] uppyu = uppus[1];

                double[][] Nptus;
                CalcFieldConsistentTLFrameNptusAtLocalX(u, up, upp, Nps, Npps, Npuss, Nppuss, upus, uppus, out Nptus);

                // fl
                {
                    double c1 = (1.0 + upx) * (1.0 + upx) + upy * upy;
                    double d1 = 1.0 + upx;
                    double d2 = upy;
                    double d3 = (1.0 + upx) / Math.Sqrt(c1);
                    double d4 = upy / Math.Sqrt(c1);
                    double d5 = ((1.0 + upx) * uppy - upy * uppx) / c1;
                    for (int row = 0; row < 7; row++)
                    {
                        fl[row] +=
                            detJWeight * E * Ae * d1 * Npx[row] +
                            detJWeight * E * Ae * d2 * Npy[row] -
                            detJWeight * E * Ae * d3 * Npx[row] -
                            detJWeight * E * Ae * d4 * Npy[row] +
                            detJWeight * E * I * d5 * Npt[row];
                    }
                }
                // kl
                {
                    double c1 = (1.0 + upx) * (1.0 + upx) + upy * upy;
                    double d2 = 1.0 + upx;
                    double d4 = upy;
                    double d6 = (1.0 + upx) / Math.Sqrt(c1);
                    double d8 = upy / Math.Sqrt(c1);
                    double d10 = ((1.0 + upx) * uppy - upy * uppx) / c1;
                    for (int row = 0; row < 7; row++)
                    {
                        for (int col = 0; col < 7; col++)
                        {
                            double d1 = upxu[col];
                            double d3 = upyu[col];
                            double d5 = (-(1.0 + upx) * upy * upyu[col] + upy * upy * upxu[col]) / Math.Pow(c1, 3.0 / 2.0);
                            double d7 = ((1.0 + upx) * (1.0 + upx) * upyu[col] -
                                (1.0 + upx) * upy * upxu[col]) / Math.Pow(c1, 3.0 / 2.0);
                            double d9 = ((1.0 + upx) * (1.0 + upx) * (1.0 + upx) * uppyu[col] +
                                (1.0 + upx) * (1.0 + upx) * (
                                -uppy * upxu[col] - uppx * upyu[col] - upy * uppxu[col]
                                ) +
                                (1.0 + upx) * (
                                -2.0 * uppy * upy * upyu[col] + 2.0 * upy * uppx * upxu[col] + upy * upy * uppyu[col]
                                ) +
                                upy * upy * (uppx * upyu[col] + uppy * upxu[col]) -
                                upy * upy * upy * uppxu[col]
                                ) / (c1 * c1);

                            kl[row, col] +=
                                detJWeight * E * Ae * (d1 * Npx[row] + d2 * Npxus[row][col]) +
                                detJWeight * E * Ae * (d3 * Npy[row] + d4 * Npyus[row][col]) -
                                detJWeight * E * Ae * (d5 * Npx[row] + d6 * Npxus[row][col]) -
                                detJWeight * E * Ae * (d7 * Npy[row] + d8 * Npyus[row][col]) +
                                detJWeight * E * I * (d9 * Npt[row] + d10 * Nptus[row][col]);
                        }
                    }
                }
            }
        }
    }
}
