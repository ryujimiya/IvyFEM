using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DFEMUtils
    {
        private static double[] CalcFieldConsistentTLFrameUAtLocalX(double x, double l0, double[] ul)
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

            double[] u = new double[] { ux, uy };
            return u;
        }

        private static double[] CalcFieldConsistentTLFrameUpAtLocalX(double x, double l0, double[] ul)
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
            double[] up = new double[] { upx, upy };
            return up;
        }

        private static double[] CalcFieldConsistentTLFrameUppAtLocalX(double x, double l0, double[] ul)
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

            double[] upp = new double[] { uppx, uppy };
            return upp;
        }

        private static double[][] CalcFieldConsistentTLFrameUpusAtLocalX(double x, double l0, double[] ul)
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

            double[][] upus = new double[2][] { upxu, upyu };
            return upus;
        }

        private static double[][] CalcFieldConsistentTLFrameUppusAtLocalX(double x, double l0, double[] ul)
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

            double[][] uppus = new double[2][] { uppxu, uppyu };
            return uppus;
        }

        private static double[][,] CalcFieldConsistentTLFrameUpuusAtLocalX(double x, double l0, double[] ul)
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

            //
            double[,] upxuu = new double[7, 7];

            //
            double[,] upyuu = new double[7, 7];
            double upyux1tz1 = -(1.0 / l0) * (3.0 / (costz1 * costz1)) +
                (x / (l0 * l0)) * (12.0 / (costz1 * costz1)) -
                (x * x / (l0 * l0 * l0)) * (9.0 / (costz1 * costz1));
            double upyux1tz2 = -(x / (l0 * l0)) * (2.0 / (costz2 * costz2)) +
                (x * x / (l0 * l0 * l0)) * (3.0 / (costz2 * costz2));
            upyuu[0, 2] = upyux1tz1;
            upyuu[0, 5] = upyux1tz2;
            double upytz1ux1 = -(1.0 / l0) * (3.0 / (costz1 * costz1)) +
                (x / (l0 * l0)) * (12.0 / (costz1 * costz1)) -
                (x * x / (l0 * l0 * l0)) * (9.0 / (costz1 * costz1));
            double upytz1tz1 = (1.0 / l0) * 2.0 *
                (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (sintz1 / (costz1 * costz1 * costz1)) -
                (x / (l0 * l0)) * 8.0 *
                (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (sintz1 / (costz1 * costz1 * costz1)) +
                (x * x / (l0 * l0 * l0)) * 6.0 *
                (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (sintz1 / (costz1 * costz1 * costz1));
            double upytz1ux2 = -(1.0 / l0) * (1.0 / (costz1 * costz1)) +
                (x / (l0 * l0)) * (4.0 / (costz1 * costz1)) -
                (x * x / (l0 * l0 * l0)) * (3.0 / (costz1 * costz1));
            double upytz1ux3 = (1.0 / l0) * (4.0 / (costz1 * costz1)) -
                (x / (l0 * l0)) * (16.0 / (costz1 * costz1)) +
                (x * x / (l0 * l0 * l0)) * (12.0 / (costz1 * costz1));
            upyuu[2, 0] = upytz1ux1;
            upyuu[2, 2] = upytz1tz1;
            upyuu[2, 3] = upytz1ux2;
            upyuu[2, 6] = upytz1ux3;
            double upyux2tz1 = -(1.0 / l0) * (1.0 / (costz1 * costz1)) +
                (x / (l0 * l0)) * (4.0 / (costz1 * costz1)) -
                (x * x / (l0 * l0 * l0)) * (3.0 / (costz1 * costz1));
            double upyux2tz2 = -(x / (l0 * l0)) * (6.0 / (costz2 * costz2)) +
                (x * x / (l0 * l0 * l0)) * (9.0 / (costz2 * costz2));
            upyuu[3, 2] = upyux2tz1;
            upyuu[3, 5] = upyux2tz2;
            double upytz2ux1 = -(x / (l0 * l0)) * (2.0 / (costz2 * costz2)) +
                (x * x / (l0 * l0 * l0)) * (3.0 / (costz2 * costz2));
            double upytz2ux2 = -(x / (l0 * l0)) * (6.0 / (costz2 * costz2)) +
                (x * x / (l0 * l0 * l0)) * (9.0 / (costz2 * costz2));
            double upytz2tz2 = -(x / (l0 * l0)) * 4.0 *
                (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (sintz2 / (costz2 * costz2 * costz2)) +
                (x * x / (l0 * l0 * l0)) * 6.0 *
                (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (sintz2 / (costz2 * costz2 * costz2));
            double upytz2ux3 = (x / (l0 * l0)) * (8.0 / (costz2 * costz2)) -
                (x * x / (l0 * l0 * l0)) * (12.0 / (costz2 * costz2));
            upyuu[5, 0] = upytz2ux1;
            upyuu[5, 3] = upytz2ux2;
            upyuu[5, 5] = upytz2tz2;
            upyuu[5, 6] = upytz2ux3;
            double upyux3tz1 = (1.0 / l0) * (4.0 / (costz1 * costz1)) -
                (x / (l0 * l0)) * (16.0 / (costz1 * costz1)) +
                (x * x / (l0 * l0 * l0)) * (12.0 / (costz1 * costz1));
            double upyux3tz2 = (x / (l0 * l0)) * (8.0 / (costz2 * costz2)) -
                (x * x / (l0 * l0 * l0)) * (12.0 / (costz2 * costz2));
            upyuu[6, 2] = upyux3tz1;
            upyuu[6, 5] = upyux3tz2;

            double[][,] upuus = new double[2][,] { upxuu, upyuu };
            return upuus;
        }

        private static double[][] CalcFieldConsistentTLFrameNsAtLocalX(double x, double l0, double[] ul)
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

            double[][] Ns = new double[2][] { Nx, Ny };
            return Ns;
        }

        private static double[][] CalcFieldConsistentTLFrameNpsAtLocalX(double x, double l0, double[] ul)
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

            double[][] Nps = new double[2][] { Npx, Npy };
            return Nps;
        }

        private static double[][] CalcFieldConsistentTLFrameNppsAtLocalX(double x, double l0, double[] ul)
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

            double[][] Npps = new double[2][] { Nppx, Nppy };
            return Npps;
        }

        private static double[][][] CalcFieldConsistentTLFrameNussAtLocalX(double x, double l0, double[] ul)
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
            double[][] Nxus = new double[7][];
            for (int i = 0; i < 7; i++)
            {
                Nxus[i] = new double[7];
            }

            //
            double Ny1ux1 = 0.0;
            double Ny1uy1 = 0.0;
            double Ny1tz1 = -(x / l0) * (3.0 / (costz1 * costz1)) +
                (x * x / (l0 * l0)) * (6.0 / (costz1 * costz1)) -
                (x * x * x / (l0 * l0 * l0)) * (3.0 / (costz1 * costz1));
            double Ny1ux2 = 0.0;
            double Ny1uy2 = 0.0;
            double Ny1tz2 = -(x * x / (l0 * l0)) * (1.0 / (costz2 * costz2)) +
                (x * x * x / (l0 * l0 * l0)) * (1.0 / (costz2 * costz2));
            double Ny1ux3 = 0.0;
            double[] Ny1u = { Ny1ux1, Ny1uy1, Ny1tz1, Ny1ux2, Ny1uy2, Ny1tz2, Ny1ux3 };
            //
            double[] Ny2u = new double[7];
            //
            double Ny3ux1 = -(x / l0) * (3.0 / (costz1 * costz1)) +
                (x * x / (l0 * l0)) * (6.0 / (costz1 * costz1)) -
                (x * x * x / (l0 * l0 * l0)) * (3.0 / (costz1 * costz1));
            double Ny3uy1 = 0.0;
            double Ny3tz1 = (x / l0) * 2.0 *
                (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (sintz1 / (costz1 * costz1 * costz1)) -
                (x * x / (l0 * l0)) * 4.0 *
                (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (sintz1 / (costz1 * costz1 * costz1)) +
                (x * x * x / (l0 * l0 * l0)) * 2.0 *
                (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) * (sintz1 / (costz1 * costz1 * costz1));
            double Ny3ux2 = -(x / l0) * (1.0 / (costz1 * costz1)) +
                (x * x / (l0 * l0)) * (2.0 / (costz1 * costz1)) -
                (x * x * x / (l0 * l0 * l0)) * (1.0 / (costz1 * costz1));
            double Ny3uy2 = 0.0;
            double Ny3tz2 = 0.0;
            double Ny3ux3 = (x / l0) * (4.0 / (costz1 * costz1)) -
                (x * x / (l0 * l0)) * (8.0 / (costz1 * costz1)) +
                (x * x * x / (l0 * l0 * l0)) * (4.0 / (costz1 * costz1));
            double[] Ny3u = { Ny3ux1, Ny3uy1, Ny3tz1, Ny3ux2, Ny3uy2, Ny3tz2, Ny3ux3 };
            //
            double Ny4ux1 = 0.0;
            double Ny4uy1 = 0.0;
            double Ny4tz1 = -(x / l0) * (1.0 / (costz1 * costz1)) +
                (x * x / (l0 * l0)) * (2.0 / (costz1 * costz1)) -
                (x * x * x / (l0 * l0 * l0)) * (1.0 / (costz1 * costz1));
            double Ny4ux2 = 0.0;
            double Ny4uy2 = 0.0;
            double Ny4tz2 = -(x * x / (l0 * l0)) * (3.0 / (costz2 * costz2)) +
                (x * x * x / (l0 * l0 * l0)) * (3.0 / (costz2 * costz2));
            double Ny4ux3 = 0.0;
            double[] Ny4u = { Ny4ux1, Ny4uy1, Ny4tz1, Ny4ux2, Ny4uy2, Ny4tz2, Ny4ux3 };
            //
            double[] Ny5u = new double[7];
            //
            double Ny6ux1 = -(x * x / (l0 * l0)) * (1.0 / (costz2 * costz2)) +
                (x * x * x / (l0 * l0 * l0)) * (1.0 / (costz2 * costz2));
            double Ny6uy1 = 0.0;
            double Ny6tz1 = 0.0;
            double Ny6ux2 = -(x * x / (l0 * l0)) * (3.0 / (costz2 * costz2)) +
                (x * x * x / (l0 * l0 * l0)) * (3.0 / (costz2 * costz2));
            double Ny6uy2 = 0.0;
            double Ny6tz2 = -(x * x / (l0 * l0)) * 2.0 *
                (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (sintz2 / (costz2 * costz2 * costz2)) +
                (x * x * x / (l0 * l0 * l0)) * 2.0 *
                (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) * (sintz2 / (costz2 * costz2 * costz2));
            double Ny6ux3 = (x * x / (l0 * l0)) * (4.0 / (costz2 * costz2)) -
                (x * x * x / (l0 * l0 * l0)) * (4.0 / (costz2 * costz2));
            double[] Ny6u = { Ny6ux1, Ny6uy1, Ny6tz1, Ny6ux2, Ny6uy2, Ny6tz2, Ny6ux3 };
            //
            double Ny7ux1 = 0.0;
            double Ny7uy1 = 0.0;
            double Ny7tz1 = (x / l0) * (4.0 / (costz1 * costz1)) -
                (x * x / (l0 * l0)) * (8.0 / (costz1 * costz1)) +
                (x * x * x / (l0 * l0 * l0)) * (4.0 / (costz1 * costz1));
            double Ny7ux2 = 0.0;
            double Ny7uy2 = 0.0;
            double Ny7tz2 = (x * x / (l0 * l0)) * (4.0 / (costz2 * costz2)) -
                (x * x * x / (l0 * l0 * l0)) * (4.0 / (costz2 * costz2));
            double Ny7ux3 = 0.0;
            double[] Ny7u = { Ny7ux1, Ny7uy1, Ny7tz1, Ny7ux2, Ny7uy2, Ny7tz2, Ny7ux3 };

            double[][] Nyus = { Ny1u, Ny2u, Ny3u, Ny4u, Ny5u, Ny6u, Ny7u };

            double[][][] Nuss = new double[2][][] { Nxus, Nyus };
            return Nuss;
        }

        private static double[][][] CalcFieldConsistentTLFrameNpussAtLocalX(double x, double l0, double[] ul)
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

            double[][][] Npuss = new double[2][][] { Npxus, Npyus };
            return Npuss;
        }

        private static double[][][] CalcFieldConsistentTLFrameNppussAtLocalX(double x, double l0, double[] ul)
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

            double[][][] Nppuss = new double[2][][] { Nppxus, Nppyus };
            return Nppuss;
        }

        private static double[][][,] CalcFieldConsistentTLFrameNuussAtLocalX(double x, double l0, double[] ul)
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
            double[][,] Nxuus = new double[7][,];
            for (int i = 0; i < 7; i++)
            {
                Nxuus[i] = new double[7, 7];
            }

            //
            double[,] Ny1uu = new double[7, 7];
            double Ny1tz1tz1 = -(x / l0) * 6.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x * x / (l0 * l0)) * 12.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x * x * x / (l0 * l0 * l0)) * 6.0 * (sintz1 / (costz1 * costz1 * costz1));
            Ny1uu[2, 2] = Ny1tz1tz1;

            double Ny1tz2tz2 = -(x * x / (l0 * l0)) * 2.0 * (sintz2 / (costz2 * costz2 * costz2)) +
                (x * x * x / (l0 * l0 * l0)) * 2.0 * (sintz2 / (costz2 * costz2 * costz2));
            Ny1uu[5, 5] = Ny1tz2tz2;
            //
            double[,] Ny2uu = new double[7, 7];
            //
            double[,] Ny3uu = new double[7, 7];
            double Ny3ux1tz1 = -(x / l0) * 6.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x * x / (l0 * l0)) * 12.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x * x * x / (l0 * l0 * l0)) * 6.0 * (sintz1 / (costz1 * costz1 * costz1));
            Ny3uu[0, 2] = Ny3ux1tz1;
            double Ny3tz1ux1 = -(x / l0) * 6.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x * x / (l0 * l0)) * 12.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x * x * x / (l0 * l0 * l0)) * 6.0 * (sintz1 / (costz1 * costz1 * costz1));
            double Ny3tz1tz1 = (x / l0) * 2.0 *
                (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) *
                ((2.0 * sintz1 * sintz1 + 1.0) / (costz1 * costz1 * costz1 * costz1)) -
                (x * x / (l0 * l0)) * 4.0 *
                (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) *
                ((2.0 * sintz1 * sintz1 + 1.0) / (costz1 * costz1 * costz1 * costz1)) +
                (x * x * x / (l0 * l0 * l0)) * 2.0 *
                (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) *
                ((2.0 * sintz1 * sintz1 + 1.0) / (costz1 * costz1 * costz1 * costz1));
            double Ny3tz1ux2 = -(x / l0) * 2.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x * x / (l0 * l0)) * 4.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x * x * x / (l0 * l0 * l0)) * 2.0 * (sintz1 / (costz1 * costz1 * costz1));
            double Ny3tz1ux3 = (x / l0) * 8.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x * x / (l0 * l0)) * 16.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x * x * x / (l0 * l0 * l0)) * 8.0 * (sintz1 / (costz1 * costz1 * costz1));
            Ny3uu[2, 0] = Ny3tz1ux1;
            Ny3uu[2, 2] = Ny3tz1tz1;
            Ny3uu[2, 3] = Ny3tz1ux2;
            Ny3uu[2, 6] = Ny3tz1ux3;
            double Ny3ux2tz1 = -(x / l0) * 2.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x * x / (l0 * l0)) * 4.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x * x * x / (l0 * l0 * l0)) * 2.0 * (sintz1 / (costz1 * costz1 * costz1));
            Ny3uu[3, 2] = Ny3ux2tz1;
            double Ny3ux3tz1 = (x / l0) * 8.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x * x / (l0 * l0)) * 16.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x * x * x / (l0 * l0 * l0)) * 8.0 * (sintz1 / (costz1 * costz1 * costz1));
            Ny3uu[6, 2] = Ny3ux3tz1;
            //
            double[,] Ny4uu = new double[7, 7];
            double Ny4tz1tz1 = -(x / l0) * 2.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x * x / (l0 * l0)) * 4.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x * x * x / (l0 * l0 * l0)) * 2.0 * (sintz1 / (costz1 * costz1 * costz1));
            Ny4uu[2, 2] = Ny4tz1tz1;
            double Ny4tz2tz2 = -(x * x / (l0 * l0)) * 6.0 * (sintz2 / (costz2 * costz2 * costz2)) +
                (x * x * x / (l0 * l0 * l0)) * 6.0 * (sintz2 / (costz2 * costz2 * costz2));
            Ny4uu[5, 5] = Ny4tz2tz2;
            //
            double[,] Ny5uu = new double[7, 7];
            //
            double[,] Ny6uu = new double[7, 7];
            double Ny6ux1tz2 = -(x * x / (l0 * l0)) * 2.0 * (sintz2 / (costz2 * costz2 * costz2)) +
                (x * x * x / (l0 * l0 * l0)) * 2.0 * (sintz2 / (costz2 * costz2 * costz2));
            Ny6uu[0, 5] = Ny6ux1tz2;
            double Ny6ux2tz2 = -(x * x / (l0 * l0)) * 6.0 * (sintz2 / (costz2 * costz2 * costz2)) +
                (x * x * x / (l0 * l0 * l0)) * 6.0 * (sintz2 / (costz2 * costz2 * costz2));
            Ny6uu[3, 5] = Ny6ux2tz2;
            double Ny6tz2ux1 = -(x * x / (l0 * l0)) * 2.0 * (sintz2 / (costz2 * costz2 * costz2)) +
                (x * x * x / (l0 * l0 * l0)) * 2.0 * (sintz2 / (costz2 * costz2 * costz2));
            double Ny6tz2ux2 = -(x * x / (l0 * l0)) * 6.0 * (sintz2 / (costz2 * costz2 * costz2)) +
                (x * x * x / (l0 * l0 * l0)) * 6.0 * (sintz2 / (costz2 * costz2 * costz2));
            double Ny6tz2tz2 = -(x * x / (l0 * l0)) * 2.0 *
                (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) *
                ((2.0 * sintz2 * sintz2 + 1.0) / (costz2 * costz2 * costz2 * costz2)) +
                (x * x * x / (l0 * l0 * l0)) * 2.0 *
                (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) *
                ((2.0 * sintz2 * sintz2 + 1.0) / (costz2 * costz2 * costz2 * costz2));
            double Ny6tz2ux3 = (x * x / (l0 * l0)) * 8.0 * (sintz2 / (costz2 * costz2 * costz2)) -
                (x * x * x / (l0 * l0 * l0)) * 8.0 * (sintz2 / (costz2 * costz2 * costz2));
            Ny6uu[5, 0] = Ny6tz2ux1;
            Ny6uu[5, 3] = Ny6tz2ux2;
            Ny6uu[5, 5] = Ny6tz2tz2;
            Ny6uu[5, 6] = Ny6tz2ux3;
            double Ny6ux3tz2 = (x * x / (l0 * l0)) * 8.0 * (sintz2 / (costz2 * costz2 * costz2)) -
                (x * x * x / (l0 * l0 * l0)) * 8.0 * (sintz2 / (costz2 * costz2 * costz2));
            Ny6uu[6, 5] = Ny6ux3tz2;
            //
            double[,] Ny7uu = new double[7, 7];
            double Ny7tz1tz1 = (x / l0) * 8.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x * x / (l0 * l0)) * 16.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x * x * x / (l0 * l0 * l0)) * 8.0 * (sintz1 / (costz1 * costz1 * costz1));
            Ny7uu[2, 2] = Ny7tz1tz1;
            double Ny7tz2tz2 = (x * x / (l0 * l0)) * 8.0 * (sintz2 / (costz2 * costz2 * costz2)) -
                (x * x * x / (l0 * l0 * l0)) * 8.0 * (sintz2 / (costz2 * costz2 * costz2));
            Ny7uu[5, 5] = Ny7tz2tz2;

            double[][,] Nyuus = { Ny1uu, Ny2uu, Ny3uu, Ny4uu, Ny5uu, Ny6uu, Ny7uu };

            double[][][,] Nuuss = new double[2][][,] { Nxuus, Nyuus };
            return Nuuss;
        }

        private static double[][][,] CalcFieldConsistentTLFrameNpuussAtLocalX(double x, double l0, double[] ul)
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

            //
            double[][,] Npxuus = new double[7][,];
            for (int i = 0; i < 7; i++)
            {
                Npxuus[i] = new double[7, 7];
            }

            //
            double[,] Npy1uu = new double[7, 7];
            double Npy1tz1tz1 = -(1.0 / l0) * 6.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x / (l0 * l0)) * 24.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x * x / (l0 * l0 * l0)) * 18.0 * (sintz1 / (costz1 * costz1 * costz1));
            Npy1uu[2, 2] = Npy1tz1tz1;
            double Npy1tz2tz2 = -(x / (l0 * l0)) * 4.0 * (sintz2 / (costz2 * costz2 * costz2)) +
                (x * x / (l0 * l0 * l0)) * 6.0 * (sintz2 / (costz2 * costz2 * costz2));
            Npy1uu[5, 5] = Npy1tz2tz2;

            double[,] Npy2uu = new double[7, 7];

            double[,] Npy3uu = new double[7, 7];
            double Npy3ux1tz1 = -(1.0 / l0) * 6.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x / (l0 * l0)) * 24.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x * x / (l0 * l0 * l0)) * 18.0 * (sintz1 / (costz1 * costz1 * costz1));
            Npy3uu[0, 2] = Npy3ux1tz1;
            double Npy3tz1ux1 = -(1.0 / l0) * 6.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x / (l0 * l0)) * 24.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x * x / (l0 * l0 * l0)) * 18.0 * (sintz1 / (costz1 * costz1 * costz1));
            double Npy3tz1tz1 = (1.0 / l0) * 2.0 *
                (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) *
                ((2.0 * sintz1 * sintz1 + 1.0) / (costz1 * costz1 * costz1 * costz1)) -
                (x / (l0 * l0)) * 8.0 *
                (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) *
                ((2.0 * sintz1 * sintz1 + 1.0) / (costz1 * costz1 * costz1 * costz1)) +
                (x * x / (l0 * l0 * l0)) * 6.0 *
                (l0 - 3.0 * ux1 - ux2 + 4.0 * ux3) *
                ((2.0 * sintz1 * sintz1 + 1.0) / (costz1 * costz1 * costz1 * costz1));
            double Npy3tz1ux2 = -(1.0 / l0) * 2.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x / (l0 * l0)) * 8.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x * x / (l0 * l0 * l0)) * 6.0 * (sintz1 / (costz1 * costz1 * costz1));
            double Npy3tz1ux3 = (1.0 / l0) * 8.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x / (l0 * l0)) * 32.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x * x / (l0 * l0 * l0)) * 24.0 * (sintz1 / (costz1 * costz1 * costz1));
            Npy3uu[2, 0] = Npy3tz1ux1;
            Npy3uu[2, 2] = Npy3tz1tz1;
            Npy3uu[2, 3] = Npy3tz1ux2;
            Npy3uu[2, 6] = Npy3tz1ux3;
            double Npy3ux2tz1 = -(1.0 / l0) * 2.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x / (l0 * l0)) * 8.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x * x / (l0 * l0 * l0)) * 6.0 * (sintz1 / (costz1 * costz1 * costz1));
            Npy3uu[3, 2] = Npy3ux2tz1;
            double Npy3ux3tz1 = (1.0 / l0) * 8.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x / (l0 * l0)) * 32.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x * x / (l0 * l0 * l0)) * 24.0 * (sintz1 / (costz1 * costz1 * costz1));
            Npy3uu[6, 2] = Npy3ux3tz1;

            double[,] Npy4uu = new double[7, 7];
            double Npy4tz1tz1 = -(1.0 / l0) * 2.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x / (l0 * l0)) * 8.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x * x / (l0 * l0 * l0)) * 6.0 * (sintz1 / (costz1 * costz1 * costz1));
            Npy4uu[2, 2] = Npy4tz1tz1;
            double Npy4tz2tz2 = -(x / (l0 * l0)) * 12.0 * (sintz2 / (costz2 * costz2 * costz2)) -
                (x * x / (l0 * l0 * l0)) * 18.0 * (sintz2 / (costz2 * costz2 * costz2));
            Npy4uu[5, 5] = Npy4tz2tz2;

            double[,] Npy5uu = new double[7, 7];

            double[,] Npy6uu = new double[7, 7];
            double Npy6ux1tz2 = -(x / (l0 * l0)) * 4.0 * (sintz2 / (costz2 * costz2 * costz2)) +
                (x * x / (l0 * l0 * l0)) * 6.0 * (sintz2 / (costz2 * costz2 * costz2));
            Npy6uu[0, 5] = Npy6ux1tz2;
            double Npy6ux2tz2 = -(x / (l0 * l0)) * 12.0 * (sintz2 / (costz2 * costz2 * costz2)) +
                (x * x / (l0 * l0 * l0)) * 18.0 * (sintz2 / (costz2 * costz2 * costz2));
            Npy6uu[3, 5] = Npy6ux2tz2;
            double Npy6tz2ux1 = -(x / (l0 * l0)) * 4.0 * (sintz2 / (costz2 * costz2 * costz2)) +
                (x * x / (l0 * l0 * l0)) * 6.0 * (sintz2 / (costz2 * costz2 * costz2));
            double Npy6tz2ux2 = -(x / (l0 * l0)) * 12.0 * (sintz2 / (costz2 * costz2 * costz2)) +
                (x * x / (l0 * l0 * l0)) * 18.0 * (sintz2 / (costz2 * costz2 * costz2));
            double Npy6tz2tz2 = -(x / (l0 * l0)) * 4.0 *
                (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) *
                ((2.0 * sintz2 * sintz2 + 1.0) / (costz2 * costz2 * costz2 * costz2)) +
                (x * x / (l0 * l0 * l0)) * 6.0 *
                (l0 + ux1 + 3.0 * ux2 - 4.0 * ux3) *
                ((2.0 * sintz2 * sintz2 + 1.0) / (costz2 * costz2 * costz2 * costz2));
            double Npy6tz2ux3 = (x / (l0 * l0)) * 16.0 * (sintz2 / (costz2 * costz2 * costz2)) -
                (x * x / (l0 * l0 * l0)) * 24.0 * (sintz2 / (costz2 * costz2 * costz2));
            Npy6uu[5, 0] = Npy6tz2ux1;
            Npy6uu[5, 3] = Npy6tz2ux2;
            Npy6uu[5, 5] = Npy6tz2tz2;
            Npy6uu[5, 6] = Npy6tz2ux3;
            double Npy6ux3tz2 = (x / (l0 * l0)) * 16.0 * (sintz2 / (costz2 * costz2 * costz2)) -
                (x * x / (l0 * l0 * l0)) * 24.0 * (sintz2 / (costz2 * costz2 * costz2));
            Npy6uu[6, 5] = Npy6ux3tz2;

            double[,] Npy7uu = new double[7, 7];
            double Npy7tz1tz1 = (1.0 / l0) * 8.0 * (sintz1 / (costz1 * costz1 * costz1)) -
                (x / (l0 * l0)) * 32.0 * (sintz1 / (costz1 * costz1 * costz1)) +
                (x * x / (l0 * l0 * l0)) * 24.0 * (sintz1 / (costz1 * costz1 * costz1));
            Npy7uu[2, 2] = Npy7tz1tz1;
            double Npy7tz2tz2 = (x / (l0 * l0)) * 16.0 * (sintz2 / (costz2 * costz2 * costz2)) -
                (x * x / (l0 * l0 * l0)) * 24.0 * (sintz2 / (costz2 * costz2 * costz2));
            Npy7uu[5, 5] = Npy7tz2tz2;

            double[][,] Npyuus = { Npy1uu, Npy2uu, Npy3uu, Npy4uu, Npy5uu, Npy6uu, Npy7uu };

            double[][][,] Npuuss = new double[2][][,] { Npxuus, Npyuus };
            return Npuuss;
        }


        private static double[] CalcFieldConsistentTLFrameNtAtLocalX(double[] up, double[][] Nps)
        {
            double upx = up[0];
            double upy = up[1];
            double[] Npx = Nps[0];
            double[] Npy = Nps[1];

            //
            double c1 = (1.0 + upx) * (1.0 + upx) + upy * upy;
            double[] Nt = new double[7];
            {
                double d1 = upy / c1;
                double d2 = (1.0 + upx) / c1;
                for (int i = 0; i < 7; i++)
                {
                    Nt[i] = -d1 * Npx[i] + d2 * Npy[i];
                }
            }
            return Nt;
        }

        private static double[] CalcFieldConsistentTLFrameNptAtLocalX(
            double[] up, double[] upp, double[][] Nps, double[][] Npps)
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
            // 
            double[] Npt = new double[7];
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
            return Npt;
        }

        private static double[][] CalcFieldConsistentTLFrameNtusAtLocalX(
            double[] up, double[][] upus, double[][] Nps, double[][][] Npuss)
        {
            double upx = up[0];
            double upy = up[1];
            double[] upxu = upus[0];
            double[] upyu = upus[1];
            double[] Npx = Nps[0];
            double[] Npy = Nps[1];
            double[][] Npxus = Npuss[0];
            double[][] Npyus = Npuss[1];

            double[][] Ntus = new double[7][];

            double c1 = (1.0 + upx) * (1.0 + upx) + upy * upy;
            double d2 = (1 / c1) * upy;
            double d4 = (1.0 / c1) * (1.0 + upx);
            for (int iNode = 0; iNode < 7; iNode++)
            {
                double[] Ntu = new double[7];
                Ntus[iNode] = Ntu;
                for (int row = 0; row < 7; row++)
                {
                    double d1 = -(1.0 / (c1 * c1)) * (
                        2.0 * (1.0 + upx) * upxu[row] + 2.0 * upy * upyu[row]
                        ) * upy +
                        (1.0 / c1) * upyu[row];
                    double d3 = -(1.0 / (c1 * c1)) * (
                        2.0 * (1.0 + upx) * upxu[row] + 2.0 * upy * upyu[row]
                        ) * (1.0 + upx) +
                        (1.0 / c1) * upxu[row];

                    Ntu[row] = -d1 * Npx[iNode] - d2 * Npxus[iNode][row] +
                        d3 * Npy[iNode] + d4 * Npyus[iNode][row];
                }
            }
            return Ntus;
        }

        private static double[][] CalcFieldConsistentTLFrameNptusAtLocalX(
            double[] up, double[] upp, double[][] upus, double[][] uppus,
            double[][] Nps, double[][] Npps,
            double[][][] Npuss, double[][][] Nppuss)
        {
            double upx = up[0];
            double upy = up[1];
            double uppx = upp[0];
            double uppy = upp[1];
            double[] upxu = upus[0];
            double[] upyu = upus[1];
            double[] uppxu = uppus[0];
            double[] uppyu = uppus[1];
            double[] Npx = Nps[0];
            double[] Npy = Nps[1];
            double[] Nppx = Npps[0];
            double[] Nppy = Npps[1];
            double[][] Npxus = Npuss[0];
            double[][] Npyus = Npuss[1];
            double[][] Nppxus = Nppuss[0];
            double[][] Nppyus = Nppuss[1];

            double[][] Nptus = new double[7][];

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
            for (int iNode = 0; iNode < 7; iNode++)
            {
                double[] Nptu = new double[7];
                Nptus[iNode] = Nptu;
                for (int row = 0; row < 7; row++)
                {
                    double d1 = -2.0 * (1.0 / (c1 * c1 * c1)) * (
                        2.0 * (1.0 + upx) * upxu[row] + 2.0 * upy * upyu[row]
                        ) * (
                        -2.0 * (1.0 + upx) * uppx * upy +
                        (1.0 + upx) * (1.0 + upx) * uppy -
                        upy * upy * uppy
                        ) +
                        (1.0 / (c1 * c1)) * (
                        -2.0 * upxu[row] * uppx * upy -
                        2.0 * (1.0 + upx) * uppxu[row] * upy -
                        2.0 * (1.0 + upx) * uppx * upyu[row] +
                        2.0 * (1.0 + upx) * upxu[row] * uppy +
                        (1.0 + upx) * (1.0 + upx) * uppyu[row] -
                        2.0 * upy * upyu[row] * uppy -
                        upy * upy * uppyu[row]
                        );
                    double d3 = -(1.0 / (c1 * c1)) * (
                        2.0 * (1.0 + upx) * upxu[row] + 2.0 * upy * upyu[row]
                        ) * upy +
                        (1.0 / c1) * upyu[row];
                    double d5 = -2.0 * (1.0 / (c1 * c1 * c1)) * (
                        2.0 * (1.0 + upx) * upxu[row] + 2.0 + upy * upyu[row]
                        ) * (
                        -(1.0 + upx) * (1.0 + upx) * uppx -
                        2.0 * upy * uppy * (1.0 + upx) +
                        upy * upy * uppx
                        ) +
                        (1.0 / (c1 * c1)) * (
                        -2.0 * (1.0 + upx) * upxu[row] * uppx -
                        (1.0 + upx) * (1.0 + upx) * uppxu[row] -
                        2.0 * upyu[row] * uppy * (1.0 + upx) -
                        2.0 * upy * uppyu[row] * (1.0 + upx) -
                        2.0 * upy * uppy * upxu[row] +
                        2.0 * upy * upyu[row] * uppx +
                        upy * upy * uppxu[row]
                        );
                    double d7 = -(1.0 / (c1 * c1)) * (
                        2.0 * (1.0 + upx) * upxu[row] + 2.0 * upy * upyu[row]
                        ) * (1.0 + upx) +
                        (1.0 / c1) * upxu[row];

                    Nptu[row] = -d1 * Npx[iNode] - d2 * Npxus[iNode][row] -
                        d3 * Nppx[iNode] - d4 * Nppxus[iNode][row] +
                        d5 * Npy[iNode] + d6 * Npyus[iNode][row] +
                        d7 * Nppy[iNode] + d8 * Nppyus[iNode][row];
                }
            }
            return Nptus;
        }

        private static double[][,] CalcFieldConsistentTLFrameNtuusAtLocalX(
            double[] up, double[][] upus, double[][,] upuus,
            double[][] Nps, double[][][] Npuss, double[][][,] Npuuss)
        {
            double upx = up[0];
            double upy = up[1];
            double[] upxu = upus[0];
            double[] upyu = upus[1];
            double[,] upxuu = upuus[0];
            double[,] upyuu = upuus[1];
            double[] Npx = Nps[0];
            double[] Npy = Nps[1];
            double[][] Npxus = Npuss[0];
            double[][] Npyus = Npuss[1];
            double[][,] Npxuus = Npuuss[0];
            double[][,] Npyuus = Npuuss[1];

            double[][,] Ntuus = new double[7][,];

            double c1 = (1.0 + upx) * (1.0 + upx) + upy * upy;
            double d4 = (1.0 / c1) * upy;
            double d8 = (1.0 / c1) * (1.0 + upx);
            for (int iNode = 0; iNode < 7; iNode++)
            {
                double[,] Ntuu = new double[7, 7];
                Ntuus[iNode] = Ntuu;
                for (int l = 0; l < 7; l++)
                {
                    double d2 = -(1.0 / (c1 * c1)) * (
                        2.0 * (1.0 + upx) * upxu[l] + 2.0 * upy * upyu[l]
                        ) * upy +
                        (1.0 / c1) * upyu[l];
                    double d6 = -(1.0 / (c1 * c1)) * (
                        2.0 * (1.0 + upx) * upxu[l] + 2.0 * upy * upyu[l]
                        ) * (1.0 + upx) +
                        (1.0 / c1) * upxu[l];

                    for (int m = 0; m < 7; m++)
                    {
                        double d1 = 2.0 * (1.0 / (c1 * c1 * c1)) *
                            2.0 * c1 * (2.0 * (1.0 + upx) * upxu[l] + 2.0 * upy * upyu[l]) *
                            (2.0 * (1.0 + upx) * upxu[m] + 2.0 * upy * upyu[m]) * upy -
                            (1.0 / (c1 * c1)) *
                            (2.0 * upxu[l] * upxu[m] + 2.0 * (1.0 + upx) * upxuu[l, m] +
                            2.0 * upyu[l] * upyu[m] + 2.0 * upy * upyuu[l, m]) * upy -
                            (1.0 / (c1 * c1)) * (2.0 * (1.0 + upx) * upxu[m] + 2.0 * upy * upyu[m]) * upyu[l] -
                            (1.0 / (c1 * c1)) * (2.0 * (1.0 + upx) * upxu[l] + 2.0 * upy * upyu[l]) * upyu[m] +
                            (1.0 / c1) * upyuu[l, m];
                        double d3 = -(1.0 / (c1 * c1)) * (
                            2.0 * (1.0 + upx) * upxu[m] + 2.0 * upy * upyu[m]
                            ) * upy +
                            (1.0 / c1) * upyu[m];
                        double d5 = 2.0 * (1.0 / (c1 * c1 * c1)) *
                            2.0 * c1 * (2.0 * (1.0 + upx) * upxu[l] + 2.0 * upy * upyu[l]) *
                            (2.0 * (1.0 + upx) * upxu[m] + 2.0 * upy * upyu[m]) * (1.0 + upx) -
                            (1.0 / (c1 * c1)) *
                            (2.0 * upxu[l] * upxu[m] + 2.0 * (1.0 + upx) * upxuu[l, m] +
                            2.0 * upyu[l] * upyu[m] + 2.0 * upy * upyuu[l, m]) * (1.0 + upx) -
                            (1.0 / (c1 * c1)) *
                            (2.0 * (1.0 + upx) * upxu[m] + 2.0 * upy * upyu[m]) * upxu[l] -
                            (1.0 / (c1 * c1)) *
                            (2.0 * (1.0 + upx) * upxu[l] + 2.0 * upy * upyu[l]) * upxu[m] +
                            (1.0 / c1) * upxuu[l, m];
                        double d7 = -(1.0 / (c1 * c1)) * (
                            2.0 * (1.0 + upx) * upxu[m] + 2.0 * upy * upyu[m]
                            ) * (1.0 + upx) +
                            (1.0 / c1) * upxu[m];

                        Ntuu[l, m] = -d1 * Npx[iNode] - d2 * Npxus[iNode][m] -
                            d3 * Npxus[iNode][l] - d4 * Npxuus[iNode][l, m] +
                            d5 * Npy[iNode] + d6 * Npyus[iNode][m] +
                            d7 * Npyus[iNode][l] + d8 * Npyuus[iNode][l, m];
                    }
                }
            }
            return Ntuus;
        }

        public static void CalcFieldConsistentTLFrameKl(
            double l0, double E, double Ae, double Iz, double[] ul,
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

                double[] u = CalcFieldConsistentTLFrameUAtLocalX(x, l0, ul);
                double[] up = CalcFieldConsistentTLFrameUpAtLocalX(x, l0, ul);
                double[] upp = CalcFieldConsistentTLFrameUppAtLocalX(x, l0, ul);
                double ux = u[0];
                double uy = u[1];
                double upx = up[0];
                double upy = up[1];
                double uppx = upp[0];
                double uppy = upp[1];
                double[][] upus = CalcFieldConsistentTLFrameUpusAtLocalX(x, l0, ul);
                double[][] uppus = CalcFieldConsistentTLFrameUppusAtLocalX(x, l0, ul);
                double[] upxu = upus[0];
                double[] upyu = upus[1];
                double[] uppxu = uppus[0];
                double[] uppyu = uppus[1];

                double[][] Nps = CalcFieldConsistentTLFrameNpsAtLocalX(x, l0, ul);
                double[][] Npps = CalcFieldConsistentTLFrameNppsAtLocalX(x, l0, ul);
                double[] Npx = Nps[0];
                double[] Npy = Nps[1];
                double[] Nppx = Npps[0];
                double[] Nppy = Npps[1];
                double[][][] Npuss = CalcFieldConsistentTLFrameNpussAtLocalX(x, l0, ul);
                double[][][] Nppuss = CalcFieldConsistentTLFrameNppussAtLocalX(x, l0, ul);
                double[][] Npxus = Npuss[0];
                double[][] Npyus = Npuss[1];
                double[][] Nppxus = Nppuss[0];
                double[][] Nppyus = Nppuss[1];

                double[] Npt = CalcFieldConsistentTLFrameNptAtLocalX(up, upp, Nps, Npps);
                double[][] Nptus = CalcFieldConsistentTLFrameNptusAtLocalX(
                    up, upp, upus, uppus, Nps, Npps, Npuss, Nppuss);

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
                            detJWeight * E * Iz * d5 * Npt[row];
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
                                detJWeight * E * Iz * (d9 * Npt[row] + d10 * Nptus[row][col]);
                        }
                    }
                }
            }
        }

        private static void CalcFieldConsistentTLFrameMAndMuAndMuu(
            double l0, double rho, double Ae, double Iz, double[] ul,
            out double[,] ml, out double[,,] mlu, out double[,,,] mluu)
        {
            ml = new double[7, 7];
            mlu = new double[7, 7, 7];
            mluu = new double[7, 7, 7, 7];

            IntegrationPoints ip = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);
            System.Diagnostics.Debug.Assert(ip.Ls.Length == 5);
            for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
            {
                double[] L = ip.Ls[ipPt];
                double lineLen = l0;
                double weight = ip.Weights[ipPt];
                double detJWeight = (lineLen / 2.0) * weight;

                double x = l0 * L[1];

                double[] u = CalcFieldConsistentTLFrameUAtLocalX(x, l0, ul);
                double[] up = CalcFieldConsistentTLFrameUpAtLocalX(x, l0, ul);
                double[] upp = CalcFieldConsistentTLFrameUppAtLocalX(x, l0, ul);
                double ux = u[0];
                double uy = u[1];
                double upx = up[0];
                double upy = up[1];
                double uppx = upp[0];
                double uppy = upp[1];
                double[][] upus = CalcFieldConsistentTLFrameUpusAtLocalX(x, l0, ul);
                double[] upxu = upus[0];
                double[] upyu = upus[1];
                double[][,] upuus = CalcFieldConsistentTLFrameUpuusAtLocalX(x, l0, ul);
                double[,] upxuu = upuus[0];
                double[,] upyuu = upuus[1];

                double[][] Ns = CalcFieldConsistentTLFrameNsAtLocalX(x, l0, ul);
                double[][] Nps = CalcFieldConsistentTLFrameNpsAtLocalX(x, l0, ul);
                double[] Nx = Ns[0];
                double[] Ny = Ns[1];
                double[] Npx = Nps[0];
                double[] Npy = Nps[1];
                double[][][] Nuss = CalcFieldConsistentTLFrameNussAtLocalX(x, l0, ul);
                double[][][] Npuss = CalcFieldConsistentTLFrameNpussAtLocalX(x, l0, ul);
                double[][] Nxus = Nuss[0];
                double[][] Nyus = Nuss[1];
                double[][] Npxus = Npuss[0];
                double[][] Npyus = Npuss[1];
                double[][][,] Nuuss = CalcFieldConsistentTLFrameNuussAtLocalX(x, l0, ul);
                double[][][,] Npuuss = CalcFieldConsistentTLFrameNpuussAtLocalX(x, l0, ul);
                double[][,] Nxuus = Nuuss[0];
                double[][,] Nyuus = Nuuss[1];
                double[][,] Npxuus = Npuuss[0];
                double[][,] Npyuus = Npuuss[1];

                double[] Nt = CalcFieldConsistentTLFrameNtAtLocalX(up, Nps);
                double[][] Ntus = CalcFieldConsistentTLFrameNtusAtLocalX(up, upus, Nps, Npuss);
                double[][,] Ntuus = CalcFieldConsistentTLFrameNtuusAtLocalX(up, upus, upuus, Nps, Npuss, Npuuss);

                // M
                for (int i = 0; i < 7; i++)
                {
                    for (int j = 0; j < 7; j++)
                    {
                        double mxValue = detJWeight * rho * Ae * Nx[i] * Nx[j];
                        double myValue = detJWeight * rho * Ae * Ny[i] * Ny[j];
                        double mtValue = detJWeight * rho * Iz * Nt[i] * Nt[j];
                        ml[i, j] += mxValue + myValue + mtValue;
                    }
                }
                // Mu
                for (int i = 0; i < 7; i++)
                {
                    for (int j = 0; j < 7; j++)
                    {
                        for (int l = 0; l < 7; l++)
                        {
                            double mxuValue = detJWeight * rho * Ae * (Nxus[i][l] * Nx[j] + Nx[i] * Nxus[j][l]);
                            double myuValue = detJWeight * rho * Ae * (Nyus[i][l] * Ny[j] + Ny[i] * Nyus[j][l]);
                            double mtuVaule = detJWeight * rho * Iz * (Ntus[i][l] * Nt[j] + Nt[i] * Ntus[j][l]);
                            mlu[i, j, l] += mxuValue + myuValue + mtuVaule;
                        }
                    }
                }
                // Muu
                for (int i = 0; i < 7; i++)
                {
                    for (int j = 0; j < 7; j++)
                    {
                        for (int l = 0; l < 7; l++)
                        {
                            for (int m = 0; m < 7; m++)
                            {
                                double mxuuValue = detJWeight * rho * Ae * (
                                    Nxuus[i][l, m] * Nx[j] + Nxus[i][l] * Nxus[j][m] +
                                    Nxus[i][m] * Nxus[j][l] + Nx[i] * Nxuus[j][l, m]);
                                double myuuValue = detJWeight * rho * Ae * (
                                    Nyuus[i][l, m] * Ny[j] + Nyus[i][l] * Nyus[j][m] +
                                    Nyus[i][m] * Nyus[j][l] + Ny[i] * Nyuus[j][l, m]);
                                double mtuuValue = detJWeight * rho * Iz * (
                                    Ntuus[i][l, m] * Nt[j] + Ntus[i][l] * Ntus[j][m] +
                                    Ntus[i][m] * Ntus[j][l] + Nt[i] * Ntuus[j][l, m]);
                                mluu[i, j, l, m] += mxuuValue + myuuValue + mtuuValue;
                            }
                        }
                    }
                }
            }
        }

        public static void CalcFieldConsistentTLFrameMl(
            double l0, double rho, double Ae, double Iz,
            double[] ul, double[] velUl, double[] accUl,
            out double[]fkl,
            out IvyFEM.Lapack.DoubleMatrix ml,
            out IvyFEM.Lapack.DoubleMatrix ckl,
            out IvyFEM.Lapack.DoubleMatrix kkl)
        {
            double[,] _ml;
            double[,,] _mlu;
            double[,,,] _mluu;
            CalcFieldConsistentTLFrameMAndMuAndMuu(l0, rho, Ae, Iz, ul, out _ml, out _mlu, out _mluu);

            // fk
            fkl = new double[7];
            for (int i = 0; i < 7; i++)
            {
                double fkValue = 0.0;
                for (int m = 0; m < 7; m++)
                {
                    fkValue += _ml[i, m] * accUl[m];
                }
                for (int l = 0; l < 7; l++)
                {
                    for (int m = 0; m < 7; m++)
                    {
                        fkValue += _mlu[i, m, l] * velUl[l] * velUl[m];
                        fkValue += -(1.0 / 2.0) * velUl[l] * _mlu[l, m, i] * velUl[m];
                    }
                }
                fkl[i] = fkValue;
            }

            // ml
            ml = new IvyFEM.Lapack.DoubleMatrix(7, 7);
            for (int i = 0; i < 7; i++)
            {
                for (int j = 0; j < 7; j++)
                {
                    ml[i, j] = _ml[i, j];
                }
            }

            // ckl
            ckl = new IvyFEM.Lapack.DoubleMatrix(7, 7);
            for (int i = 0; i < 7; i++)
            {
                for (int j = 0; j < 7; j++)
                {
                    double velmValue = 0.0;
                    double c1Value = 0.0;
                    double transc1Value = 0.0;
                    for (int l = 0; l < 7; l++)
                    {
                        velmValue += _mlu[i, j, l] * velUl[l];
                        c1Value += _mlu[i, l, j] * velUl[l];
                        transc1Value += _mlu[j, l, i] * velUl[l];
                    }
                    ckl[i, j] = velmValue + c1Value + transc1Value;
                }
            }

            // kkl
            kkl = new IvyFEM.Lapack.DoubleMatrix(7, 7);
            for (int i = 0; i < 7; i++)
            {
                for (int j = 0; j < 7; j++)
                {
                    double kkl1Value = 0.0;
                    double kkl2Value = 0.0;
                    double kkl3Value = 0.0;
                    for (int l = 0; l < 7; l++)
                    {
                        kkl1Value += _mlu[i, l, j] * accUl[l];
                    }
                    for (int l = 0; l < 7; l++)
                    {
                        for (int m = 0; m < 7; m++)
                        {
                            kkl2Value += _mluu[i, m, j, l] * velUl[l] * velUl[m];
                            kkl3Value += (1.0 / 2.0) * velUl[l] * _mluu[l, m, i, j] * velUl[m];
                        }
                    }
                    kkl[i, j] = kkl1Value + kkl2Value - kkl3Value;
                }
            }
        }
    }
}
