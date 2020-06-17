using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic3DFEMUtils
    {
        protected static void CalcDKTHFunctionCoeffs(
            OpenTK.Vector2d[] pt2Ds,
            out double[] ha, out double[] hb, out double[] hc, out double[] hd, out double[] he)
        {
            // 辺0,1,2
            double[] xab = new double[3];
            double[] yab = new double[3];
            double[] lc = new double[3];
            for (int ic = 0; ic < 3; ic++)
            {
                int ia = (ic + 1) % 3;
                int ib = (ic + 2) % 3;
                OpenTK.Vector2d pta = pt2Ds[ia];
                OpenTK.Vector2d ptb = pt2Ds[ib];
                double _xab = pta.X - ptb.X;
                double _yab = pta.Y - ptb.Y;
                double _lc = Math.Sqrt(_xab * _xab + _yab * _yab);
                xab[ic] = _xab;
                yab[ic] = _yab;
                lc[ic] = _lc;
            }
            // 節点3,4,5 --> 辺2,0,1
            ha = new double[3];
            hb = new double[3];
            hc = new double[3];
            hd = new double[3];
            he = new double[3];
            for (int i = 0; i < 3; i++)
            {
                int ic = (i + 2) % 3;
                double _xab = xab[ic];
                double _yab = yab[ic];
                double _lc = lc[ic];
                ha[i] = -(1.0 / (_lc * _lc)) * _xab;
                hb[i] = (3.0 / 4.0) * (1.0 / (_lc * _lc)) * _xab * _yab;
                hc[i] = (1.0 / (_lc * _lc)) * ((1.0 / 4.0) * _xab * _xab - (1.0 / 2.0) * _yab * _yab);
                hd[i] = -(1.0 / (_lc * _lc)) * _yab;
                he[i] = (1.0 / (_lc * _lc)) * (-(1.0 / 2.0) * _xab * _xab + (1.0 / 4.0) * _yab * _yab);
            }
        }

        protected static void CalcDKTHFunctions(
            double[] N, double[] Nx, double[] Ny,
            double[] ha, double[] hb, double[] hc, double[] hd, double[] he,
            out double[] Hx, out double[] Hy,
            out double[] Hxx, out double[] Hxy, out double[] Hyx, out double[] Hyy)
        {
            Hx = new double[9];
            Hy = new double[9];
            Hxx = new double[9];
            Hxy = new double[9];
            Hyx = new double[9];
            Hyy = new double[9];
            for (int i = 0; i < 3; i++)
            {
                int i4 = i + 3;
                int i6 = (i + 2) % 3 + 3;

                // 形状関数
                Hx[i * 3 + 0] = (3.0 / 2.0) * (ha[i4 - 3] * N[i4] - ha[i6 - 3] * N[i6]);
                Hx[i * 3 + 1] = hb[i4 - 3] * N[i4] + hb[i6 - 3] * N[i6];
                Hx[i * 3 + 2] = N[i] - hc[i4 - 3] * N[i4] - hc[i6 - 3] * N[i6];
                //
                Hy[i * 3 + 0] = (3.0 / 2.0) * (hd[i4 - 3] * N[i4] - hd[i6 - 3] * N[i6]);
                Hy[i * 3 + 1] = -N[i] + he[i4 - 3] * N[i4] + he[i6 - 3] * N[i6];
                Hy[i * 3 + 2] = -hb[i4 - 3] * N[i4] - hb[i6 - 3] * N[i6];

                // 微分
                Hxx[i * 3 + 0] = (3.0 / 2.0) * (ha[i4 - 3] * Nx[i4] - ha[i6 - 3] * Nx[i6]);
                Hxx[i * 3 + 1] = hb[i4 - 3] * Nx[i4] + hb[i6 - 3] * Nx[i6];
                Hxx[i * 3 + 2] = Nx[i] - hc[i4 - 3] * Nx[i4] - hc[i6 - 3] * Nx[i6];
                Hxy[i * 3 + 0] = (3.0 / 2.0) * (ha[i4 - 3] * Ny[i4] - ha[i6 - 3] * Ny[i6]);
                Hxy[i * 3 + 1] = hb[i4 - 3] * Ny[i4] + hb[i6 - 3] * Ny[i6];
                Hxy[i * 3 + 2] = Ny[i] - hc[i4 - 3] * Ny[i4] - hc[i6 - 3] * Ny[i6];
                //
                Hyx[i * 3 + 0] = (3.0 / 2.0) * (hd[i4 - 3] * Nx[i4] - hd[i6 - 3] * Nx[i6]);
                Hyx[i * 3 + 1] = -Nx[i] + he[i4 - 3] * Nx[i4] + he[i6 - 3] * Nx[i6];
                Hyx[i * 3 + 2] = -hb[i4 - 3] * Nx[i4] - hb[i6 - 3] * Nx[i6];
                Hyy[i * 3 + 0] = (3.0 / 2.0) * (hd[i4 - 3] * Ny[i4] - hd[i6 - 3] * Ny[i6]);
                Hyy[i * 3 + 1] = -Ny[i] + he[i4 - 3] * Ny[i4] + he[i6 - 3] * Ny[i6];
                Hyy[i * 3 + 2] = -hb[i4 - 3] * Ny[i4] - hb[i6 - 3] * Ny[i6];
            }
        }

        public static IvyFEM.Lapack.DoubleMatrix CalcDKTPlateKl(
            TriangleFE d1TriFE, double h, IvyFEM.Lapack.DoubleMatrix C)
        {
            double Ae = d1TriFE.GetArea();
            // DKT要素の補間用に2次三角形要素を作る
            int[] vertexCoIds = d1TriFE.VertexCoordIds;
            int[] nodeCoIds = new int[6];
            for (int i = 0; i < 3; i++)
            {
                nodeCoIds[i] = vertexCoIds[i];
                nodeCoIds[i + 3] = -1; // 中点の節点は存在しない
            }
            TriangleFE triFE2nd = new TriangleFE(2, FiniteElementType.ScalarLagrange);
            triFE2nd.QuantityId = d1TriFE.QuantityId;
            triFE2nd.World = d1TriFE.World;
            triFE2nd.SetVertexCoordIds(vertexCoIds);
            triFE2nd.SetNodeCoordIds(nodeCoIds);

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
                Kml = IvyFEM.Lapack.DoubleMatrix.Scal((transBm * C * Bm), (h * Ae));
            }

            //------------------------------------
            // DKT bending element
            IvyFEM.Lapack.DoubleMatrix Kbl = new IvyFEM.Lapack.DoubleMatrix(9, 9);
            var Db = IvyFEM.Lapack.DoubleMatrix.Scal(C, (h * h * h / 12.0));
            OpenTK.Vector2d[] pt2Ds = triFE2nd.ProjectVertexsFrom3D();
            double[] ha;
            double[] hb;
            double[] hc;
            double[] hd;
            double[] he;
            CalcDKTHFunctionCoeffs(pt2Ds, out ha, out hb, out hc, out hd, out he);

            IntegrationPoints ip = TriangleFE.GetIntegrationPoints(TriangleIntegrationPointCount.Point3);
            System.Diagnostics.Debug.Assert(ip.Ls.Length == 3);
            for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
            {
                double[] L = ip.Ls[ipPt];
                double[] N = triFE2nd.CalcN(L);
                double[][] Nu = triFE2nd.CalcNu(L);
                double[] Nx = Nu[0];
                double[] Ny = Nu[1];
                double detJ = triFE2nd.GetDetJacobian(L);
                double weight = ip.Weights[ipPt];
                double detJWeight = (1.0 / 2.0) * weight * detJ;

                double[] Hx;
                double[] Hy;
                double[] Hxx;
                double[] Hxy;
                double[] Hyx;
                double[] Hyy;
                CalcDKTHFunctions(N, Nx, Ny, ha, hb, hc, hd, he, out Hx, out Hy, out Hxx, out Hxy, out Hyx, out Hyy);

                IvyFEM.Lapack.DoubleMatrix Bb = new IvyFEM.Lapack.DoubleMatrix(3, 9);
                for (int j = 0; j < 9; j++)
                {
                    Bb[0, j] = Hxx[j];
                    Bb[1, j] = Hyy[j];
                    Bb[2, j] = Hxy[j] + Hyx[j];
                }
                var transBb = IvyFEM.Lapack.DoubleMatrix.Transpose(Bb);
                var Kblip = transBb * Db * Bb;
                for (int i = 0; i < 9; i++)
                {
                    for (int j = 0; j < 9; j++)
                    {
                        Kbl[i, j] += detJWeight * Kblip[i, j];
                    }
                }
            }

            //------------------------------------
            // CST + DKT
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
            for (int iNode = 0; iNode < 3; iNode++)
            {
                // fictitious stiffness value
                double f = 0.0;
                for (int k = 0; k < (mdof + bdof); k++)
                {
                    double diagVal = Kl[iNode * dof + k, iNode * dof + k];
                    System.Diagnostics.Debug.Assert(diagVal >= IvyFEM.Constants.PrecisionLowerLimit);
                    if (diagVal > f)
                    {
                        f = diagVal;
                    }
                }
                f *= 1.0e-3;
                Kl[iNode * dof + offsettz, iNode * dof + offsettz] = f;
            }
            return Kl;
        }

        public static IvyFEM.Lapack.DoubleMatrix CalcDKTPlateMl(
            TriangleFE d1TriFE, double h, double rho)
        {
            double Ae = d1TriFE.GetArea();
            // DKT要素の補間用に2次三角形要素を作る
            int[] vertexCoIds = d1TriFE.VertexCoordIds;
            int[] nodeCoIds = new int[6];
            for (int i = 0; i < 3; i++)
            {
                nodeCoIds[i] = vertexCoIds[i];
                nodeCoIds[i + 3] = -1; // 中点の節点は存在しない
            }
            TriangleFE triFE2nd = new TriangleFE(2, FiniteElementType.ScalarLagrange);
            triFE2nd.QuantityId = d1TriFE.QuantityId;
            triFE2nd.World = d1TriFE.World;
            triFE2nd.SetVertexCoordIds(vertexCoIds);
            triFE2nd.SetNodeCoordIds(nodeCoIds);

            IvyFEM.Lapack.DoubleMatrix Ml = new IvyFEM.Lapack.DoubleMatrix(18, 18);
            OpenTK.Vector2d[] pt2Ds = triFE2nd.ProjectVertexsFrom3D();
            double[] ha;
            double[] hb;
            double[] hc;
            double[] hd;
            double[] he;
            CalcDKTHFunctionCoeffs(pt2Ds, out ha, out hb, out hc, out hd, out he);

            // CST+DKT
            double rhom = rho * h;
            double rhomb = 0.0;
            double rhob = rho * h * h * h / 12.0;
            // Point3だと不十分
            IntegrationPoints ip = TriangleFE.GetIntegrationPoints(TriangleIntegrationPointCount.Point7);
            System.Diagnostics.Debug.Assert(ip.Ls.Length == 7);
            for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
            {
                double[] L = ip.Ls[ipPt];
                double[] N = triFE2nd.CalcN(L);
                double[][] Nu = triFE2nd.CalcNu(L);
                double[] Nx = Nu[0];
                double[] Ny = Nu[1];
                double detJ = triFE2nd.GetDetJacobian(L);
                double weight = ip.Weights[ipPt];
                double detJWeight = (1.0 / 2.0) * weight * detJ;

                double[] Hx;
                double[] Hy;
                double[] Hxx;
                double[] Hxy;
                double[] Hyx;
                double[] Hyy;
                CalcDKTHFunctions(N, Nx, Ny, ha, hb, hc, hd, he, out Hx, out Hy, out Hxx, out Hxy, out Hyx, out Hyy);

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
                double[] HxVec = new double[18]
                {
                    0.0, 0.0, Hx[0], Hx[1], Hx[2], 0.0,
                    0.0, 0.0, Hx[3], Hx[4], Hx[5], 0.0,
                    0.0, 0.0, Hx[6], Hx[7], Hx[8], 0.0
                };
                double[] HyVec = new double[18]
                {
                    0.0, 0.0, Hy[0], Hy[1], Hy[2], 0.0,
                    0.0, 0.0, Hy[3], Hy[4], Hy[5], 0.0,
                    0.0, 0.0, Hy[6], Hy[7], Hy[8], 0.0
                };
                var NxVecMat = new IvyFEM.Lapack.DoubleMatrix(NxVec, 18, 1, false);
                var NyVecMat = new IvyFEM.Lapack.DoubleMatrix(NyVec, 18, 1, false);
                var NzVecMat = new IvyFEM.Lapack.DoubleMatrix(NzVec, 18, 1, false);
                var HxVecMat = new IvyFEM.Lapack.DoubleMatrix(HxVec, 18, 1, false);
                var HyVecMat = new IvyFEM.Lapack.DoubleMatrix(HyVec, 18, 1, false);
                var transNxVecMat = IvyFEM.Lapack.DoubleMatrix.Transpose(NxVecMat);
                var transNyVecMat = IvyFEM.Lapack.DoubleMatrix.Transpose(NyVecMat);
                var transNzVecMat = IvyFEM.Lapack.DoubleMatrix.Transpose(NzVecMat);
                var transHxVecMat = IvyFEM.Lapack.DoubleMatrix.Transpose(HxVecMat);
                var transHyVecMat = IvyFEM.Lapack.DoubleMatrix.Transpose(HyVecMat);

                var NxNxMat = NxVecMat * transNxVecMat;
                var NyNyMat = NyVecMat * transNyVecMat;
                var NzNzMat = NzVecMat * transNzVecMat;
                var HxNxMat = HxVecMat * transNxVecMat;
                var HyNyMat = HyVecMat * transNyVecMat;
                var NxHxMat = NxVecMat * transHxVecMat;
                var NyHyMat = NyVecMat * transHyVecMat;
                var HxHxMat = HxVecMat * transHxVecMat;
                var HyHyMat = HyVecMat * transHyVecMat;
                                
                IvyFEM.Lapack.DoubleMatrix Mlip = new IvyFEM.Lapack.DoubleMatrix(18, 18);
                for (int i = 0; i < 18; i++)
                {
                    for (int j = 0; j < 18; j++)
                    {
                        Mlip[i, j] = rhom * (NxNxMat[i, j] + NyNyMat[i, j] + NzNzMat[i, j]) +
                            rhomb * (HxNxMat[i, j] + NyNyMat[i, j] + NxHxMat[i, j] + NyHyMat[i, j]) +
                            rhob * (HxHxMat[i, j] + HyHyMat[i, j]);
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
