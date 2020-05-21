using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic3DFEMUtils
    {
        public static IvyFEM.Lapack.DoubleMatrix CalcMindlinPlateKl(
            TriangleFE d1TriFE, double h, IvyFEM.Lapack.DoubleMatrix Cb, IvyFEM.Lapack.DoubleMatrix Cs, double kappa)
        {
            double Ae = d1TriFE.GetArea();

            //------------------------------------
            // CST membrane element
            IvyFEM.Lapack.DoubleMatrix Kml;
            {
                // Nuを求める 定数なのでLはどこでもいい
                double[] L = { 0.5, 0.5, 0.5 }; // 重心
                double[][] Nu = d1TriFE.CalcNu(L);
                double[] Nx = Nu[0];
                double[] Ny = Nu[1];
                IvyFEM.Lapack.DoubleMatrix Bm = new IvyFEM.Lapack.DoubleMatrix(3, 6);
                for (int i = 0; i < 3; i++)
                {
                    Bm[0, 2 * i] = Nx[i];
                    Bm[1, 2 * i + 1] = Ny[i];
                    Bm[2, 2 * i] = Ny[i];
                    Bm[2, 2 * i + 1] = Nx[i];
                }
                var transBm = IvyFEM.Lapack.DoubleMatrix.Transpose(Bm);
                Kml = IvyFEM.Lapack.DoubleMatrix.Scal((transBm * Cb * Bm), (h * Ae));
            }

            //------------------------------------
            // Mindlin bending element
            IvyFEM.Lapack.DoubleMatrix Kbl = new IvyFEM.Lapack.DoubleMatrix(9, 9);

            // 選択的低減積分
            IntegrationPoints ipK = TriangleFE.GetIntegrationPoints(TriangleIntegrationPointCount.Point1);
            System.Diagnostics.Debug.Assert(ipK.Ls.Length == 1);
            for (int ipPt = 0; ipPt < ipK.PointCount; ipPt++)
            {
                double[] L = ipK.Ls[ipPt];
                double[] N = d1TriFE.CalcN(L);
                double[][] Nu = d1TriFE.CalcNu(L);
                double[] Nx = Nu[0];
                double[] Ny = Nu[1];
                double detJ = d1TriFE.GetDetJacobian(L);
                double weight = ipK.Weights[ipPt];
                double detJWeight = (1.0 / 2.0) * weight * detJ;

                IvyFEM.Lapack.DoubleMatrix Bs = new IvyFEM.Lapack.DoubleMatrix(2, 9);
                double[] Bs1 = { Nx[0], 0.0, N[0], Nx[1], 0.0, N[1], Nx[2], 0.0, N[2] };
                double[] Bs2 = { Ny[0], -N[0], 0.0, Ny[1], -N[1], 0.0, Ny[2], -N[2], 0.0 };
                for (int j = 0; j < 9; j++)
                {
                    Bs[0, j] = Bs1[j];
                    Bs[1, j] = Bs2[j];
                }
                var transBs = IvyFEM.Lapack.DoubleMatrix.Transpose(Bs);
                var Kslip = transBs * Cs * Bs;
                Kslip = IvyFEM.Lapack.DoubleMatrix.Scal(Kslip, kappa * h);
                for (int i = 0; i < 9; i++)
                {
                    for (int j = 0; j < 9; j++)
                    {
                        Kbl[i, j] += detJWeight * Kslip[i, j];
                    }
                }
            }

            IntegrationPoints ipM = TriangleFE.GetIntegrationPoints(TriangleIntegrationPointCount.Point3);
            System.Diagnostics.Debug.Assert(ipM.Ls.Length == 3);
            for (int ipPt = 0; ipPt < ipM.PointCount; ipPt++)
            {
                double[] L = ipM.Ls[ipPt];
                double[] N = d1TriFE.CalcN(L);
                double[][] Nu = d1TriFE.CalcNu(L);
                double[] Nx = Nu[0];
                double[] Ny = Nu[1];
                double detJ = d1TriFE.GetDetJacobian(L);
                double weight = ipM.Weights[ipPt];
                double detJWeight = (1.0 / 2.0) * weight * detJ;

                IvyFEM.Lapack.DoubleMatrix Bb = new IvyFEM.Lapack.DoubleMatrix(3, 9);
                double[] Bb1 = { 0.0, 0.0, -Nx[0], 0.0, 0.0, -Nx[1], 0.0, 0.0, -Nx[2] };
                double[] Bb2 = { 0.0, Ny[0], 0.0, 0.0, Ny[1], 0.0, 0.0, Ny[2], 0.0 };
                double[] Bb3 = { 0.0, Nx[0], -Ny[0], 0.0, Nx[1], -Ny[1], 0.0, Nx[2], -Ny[2] };
                for (int j = 0; j < 9; j++)
                {
                    Bb[0, j] = Bb1[j];
                    Bb[1, j] = Bb2[j];
                    Bb[2, j] = Bb3[j];
                }
                var transBb = IvyFEM.Lapack.DoubleMatrix.Transpose(Bb);
                var Kblip = transBb * Cb * Bb;
                Kblip = IvyFEM.Lapack.DoubleMatrix.Scal(Kblip, h * h * h / 12.0);
                for (int i = 0; i < 9; i++)
                {
                    for (int j = 0; j < 9; j++)
                    {
                        Kbl[i, j] += detJWeight * Kblip[i, j];
                    }
                }
            }

            //------------------------------------
            // CST + Mindlin
            int dof = 6;
            int mdof = 2;
            int bdof = 3;
            int tzdof = 1;
            int offsetb = mdof;
            int offsettz = mdof + bdof;
            System.Diagnostics.Debug.Assert(dof == mdof + bdof + tzdof);
            IvyFEM.Lapack.DoubleMatrix Kl = new IvyFEM.Lapack.DoubleMatrix(18, 18);
            for (int i = 0; i < mdof; i++)
            {
                for (int j = 0; j < mdof; j++)
                {
                    Kl[i, j] = Kml[i, j];
                    Kl[i, j + dof] = Kml[i, j + mdof];
                    Kl[i, j + dof * 2] = Kml[i, j + mdof * 2];
                    Kl[i + dof, j] = Kml[i + mdof, j];
                    Kl[i + dof, j + dof] = Kml[i + mdof, j + mdof];
                    Kl[i + dof, j + dof * 2] = Kml[i + mdof, j + mdof * 2];
                    Kl[i + dof * 2, j] = Kml[i + mdof * 2, j];
                    Kl[i + dof * 2, j + dof] = Kml[i + mdof * 2, j + mdof];
                    Kl[i + dof * 2, j + dof * 2] = Kml[i + mdof * 2, j + mdof * 2];
                }
            }
            for (int i = 0; i < bdof; i++)
            {
                for (int j = 0; j < bdof; j++)
                {
                    Kl[i + offsetb, j + offsetb] = Kbl[i, j];
                    Kl[i + offsetb, j + offsetb + dof] = Kbl[i, j + bdof];
                    Kl[i + offsetb, j + offsetb + dof * 2] = Kbl[i, j + bdof * 2];
                    Kl[i + offsetb + dof, j + offsetb] = Kbl[i + bdof, j];
                    Kl[i + offsetb + dof, j + offsetb + dof] = Kbl[i + bdof, j + bdof];
                    Kl[i + offsetb + dof, j + offsetb + dof * 2] = Kbl[i + bdof, j + bdof * 2];
                    Kl[i + offsetb + dof * 2, j + offsetb] = Kbl[i + bdof * 2, j];
                    Kl[i + offsetb + dof * 2, j + offsetb + dof] = Kbl[i + bdof * 2, j + bdof];
                    Kl[i + offsetb + dof * 2, j + offsetb + dof * 2] = Kbl[i + bdof * 2, j + bdof * 2];
                }
            }
            for (int i = 0; i < 3; i++)
            {
                // fictitious stiffness value
                double f = 0.0;
                for (int k = 0; k < (mdof + bdof); k++)
                {
                    double diagVal = Kl[k + i * dof, k + i * dof];
                    System.Diagnostics.Debug.Assert(diagVal >= IvyFEM.Constants.PrecisionLowerLimit);
                    if (diagVal > f)
                    {
                        f = diagVal;
                    }
                }
                f *= 1.0e-3;
                Kl[i * dof + offsettz, i * dof + offsettz] = f;
            }
            return Kl;
        }

        public static IvyFEM.Lapack.DoubleMatrix CalcMindlinPlateMl(
            TriangleFE d1TriFE, double h, double rho)
        {
            double Ae = d1TriFE.GetArea();

            IvyFEM.Lapack.DoubleMatrix Ml = new IvyFEM.Lapack.DoubleMatrix(18, 18);

            // CST+ Mindlin
            double rhom = rho * h;
            double rhob = rho * h * h * h / 12.0;
            IntegrationPoints ip = TriangleFE.GetIntegrationPoints(TriangleIntegrationPointCount.Point3);
            System.Diagnostics.Debug.Assert(ip.Ls.Length == 3);
            for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
            {
                double[] L = ip.Ls[ipPt];
                double[] N = d1TriFE.CalcN(L);
                double[][] Nu = d1TriFE.CalcNu(L);
                double[] Nx = Nu[0];
                double[] Ny = Nu[1];
                double detJ = d1TriFE.GetDetJacobian(L);
                double weight = ip.Weights[ipPt];
                double detJWeight = (1.0 / 2.0) * weight * detJ;

                double[] NxVec = new double[18]
                {
                    N[0], 0.0, 0.0, 0.0, 0.0, 0.0,
                    N[1], 0.0, 0.0, 0.0, 0.0, 0.0,
                    N[2], 0.0, 0.0, 0.0, 0.0, 0.0
                };
                double[] NyVec = new double[18]
                {
                    0.0, N[0], 0.0, 0.0, 0.0, 0.0,
                    0.0, N[1], 0.0, 0.0, 0.0, 0.0,
                    0.0, N[2], 0.0, 0.0, 0.0, 0.0
                };
                double[] NzVec = new double[18]
                {
                    0.0, 0.0, N[0], 0.0, 0.0, 0.0,
                    0.0, 0.0, N[1], 0.0, 0.0, 0.0,
                    0.0, 0.0, N[2], 0.0, 0.0, 0.0
                };
                double[] NtxVec = new double[18]
                {
                    0.0, 0.0, 0.0, N[0], 0.0, 0.0,
                    0.0, 0.0, 0.0, N[1], 0.0, 0.0,
                    0.0, 0.0, 0.0, N[2], 0.0, 0.0
                };
                double[] NtyVec = new double[18]
                {
                    0.0, 0.0, 0.0, 0.0, N[0], 0.0,
                    0.0, 0.0, 0.0, 0.0, N[1], 0.0,
                    0.0, 0.0, 0.0, 0.0, N[2], 0.0
                };
                var NxVecMat = new IvyFEM.Lapack.DoubleMatrix(NxVec, 18, 1, false);
                var NyVecMat = new IvyFEM.Lapack.DoubleMatrix(NyVec, 18, 1, false);
                var NzVecMat = new IvyFEM.Lapack.DoubleMatrix(NzVec, 18, 1, false);
                var NtxVecMat = new IvyFEM.Lapack.DoubleMatrix(NtxVec, 18, 1, false);
                var NtyVecMat = new IvyFEM.Lapack.DoubleMatrix(NtyVec, 18, 1, false);
                var transNxVecMat = IvyFEM.Lapack.DoubleMatrix.Transpose(NxVecMat);
                var transNyVecMat = IvyFEM.Lapack.DoubleMatrix.Transpose(NyVecMat);
                var transNzVecMat = IvyFEM.Lapack.DoubleMatrix.Transpose(NzVecMat);
                var transNtxVecMat = IvyFEM.Lapack.DoubleMatrix.Transpose(NtxVecMat);
                var transNtyVecMat = IvyFEM.Lapack.DoubleMatrix.Transpose(NtyVecMat);

                var NxNxMat = NxVecMat * transNxVecMat;
                var NyNyMat = NyVecMat * transNyVecMat;
                var NzNzMat = NzVecMat * transNzVecMat;
                var NtxNtxMat = NtxVecMat * transNtxVecMat;
                var NtyNtyMat = NtyVecMat * transNtyVecMat;

                IvyFEM.Lapack.DoubleMatrix Mlip = new IvyFEM.Lapack.DoubleMatrix(18, 18);
                for (int i = 0; i < 18; i++)
                {
                    for (int j = 0; j < 18; j++)
                    {
                        Mlip[i, j] = rhom * (NxNxMat[i, j] + NyNyMat[i, j] + NzNzMat[i, j]) +
                            rhob * (NtxNtxMat[i, j] + NtyNtyMat[i, j]);
                    }
                }

                for (int i = 0; i < 18; i++)
                {
                    for (int j = 0; j < 18; j++)
                    {
                        Ml[i, j] += detJWeight * Mlip[i, j];
                    }
                }
            }
            return Ml;
        }
    }
}
