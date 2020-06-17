
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Animation;

namespace IvyFEM
{
    public partial class Elastic3DFEMUtils
    {
        protected static void CalcMITCLinearDirectorVectors(
            OpenTK.Vector3d[] xPts, out OpenTK.Vector3d[] Vns, out OpenTK.Vector3d[] V1s, out OpenTK.Vector3d[] V2s)
        {
            int nodeCnt = xPts.Length;
            System.Diagnostics.Debug.Assert(nodeCnt == 3);
            OpenTK.Vector3d normal = CadUtils3D.TriNormal(xPts[0], xPts[1], xPts[2]);
            OpenTK.Vector3d e1 = new OpenTK.Vector3d(1.0, 0.0, 0.0);
            OpenTK.Vector3d e2 = new OpenTK.Vector3d(0.0, 1.0, 0.0);
            OpenTK.Vector3d e3 = new OpenTK.Vector3d(0.0, 0.0, 1.0);
            // director vectors etc
            Vns = new OpenTK.Vector3d[nodeCnt];
            V1s = new OpenTK.Vector3d[nodeCnt];
            V2s = new OpenTK.Vector3d[nodeCnt];
            for (int i = 0; i < nodeCnt; i++)
            {
                OpenTK.Vector3d vn = normal;
                OpenTK.Vector3d v1;
                OpenTK.Vector3d v2;
                if (OpenTK.Vector3d.Cross(e2, vn).Length >= Constants.PrecisionLowerLimit)
                {
                    v1 = OpenTK.Vector3d.Normalize(OpenTK.Vector3d.Cross(e2, vn));
                    v2 = OpenTK.Vector3d.Cross(vn, v1);
                }
                else
                {
                    v1 = e3;
                    v2 = e1;
                }
                Vns[i] = vn;
                V1s[i] = v1;
                V2s[i] = v2;
            }
        }

        protected static OpenTK.Vector3d[] CalcMITCLinearCovariantBasisVectors(
            OpenTK.Vector3d[] xPts, double h, OpenTK.Vector3d[] Vns, double[] r)
        {
            int nodeCnt = xPts.Length;
            System.Diagnostics.Debug.Assert(nodeCnt == 3);
            double r1 = r[0];
            double r2 = r[1];
            double r3 = r[2];
            double[] Ns = new double[3];
            Ns[0] = 1.0 - r1 - r2;
            Ns[1] = r1;
            Ns[2] = r2;
            // partial r1,r2,r3
            double[] N1r = new double[3];
            double[] N2r = new double[3];
            double[] N3r = new double[3];
            N1r[0] = -1.0;
            N1r[1] = -1.0;
            N1r[2] = 0.0;
            N2r[0] = 1.0;
            N2r[1] = 0.0;
            N2r[2] = 0.0;
            N3r[0] = 0.0;
            N3r[1] = 1.0;
            N3r[2] = 0.0;
            double[][] Nrs = new double[3][] { N1r, N2r, N3r };

            OpenTK.Vector3d[] gs = new OpenTK.Vector3d[3];
            for (int i = 0; i < 3; i++)
            {
                double deltai3 = i == 2 ? 1.0 : 0.0;
                gs[i] = new OpenTK.Vector3d();
                for (int k = 0; k < nodeCnt; k++)
                {
                    gs[i] += Nrs[k][i] * xPts[k] + (r3 / 2.0) * Nrs[k][i] * h * Vns[k] +
                        (1.0 / 2.0) * deltai3 * Ns[k] * h * Vns[k];
                }
            }
            return gs;
        }

        protected static OpenTK.Vector3d[] CalcMITCLinearContravariantBasisVectors(OpenTK.Vector3d[] gs)
        {
            var g1 = gs[0];
            var g2 = gs[1];
            var g3 = gs[2];
            double v = OpenTK.Vector3d.Dot(g1, OpenTK.Vector3d.Cross(g2, g3));
            OpenTK.Vector3d h1 = OpenTK.Vector3d.Cross(g2, g3) / v;
            OpenTK.Vector3d h2 = OpenTK.Vector3d.Cross(g3, g1) / v;
            OpenTK.Vector3d h3 = OpenTK.Vector3d.Cross(g1, g2) / v;
            OpenTK.Vector3d[] hs = { h1, h2, h3 };
            return hs;
        }

        protected static OpenTK.Vector3d[] CalcMITCLinearLocalBasicVectors(OpenTK.Vector3d[] gs)
        {
            var g1 = gs[0];
            var g2 = gs[1];
            var g3 = gs[2];
            var l3 = OpenTK.Vector3d.Normalize(g3);
            var l1 = OpenTK.Vector3d.Normalize(OpenTK.Vector3d.Cross(g2, l3));
            var l2 = OpenTK.Vector3d.Cross(l3, l1);
            OpenTK.Vector3d[] lVecs = { l1, l2, l3 };
            return lVecs;
        }

        protected static OpenTK.Vector3d[][] CalcMITCLinearUiCoeffVectors(
            int dof, int nodeCnt,
            double h, OpenTK.Vector3d[] Vns, OpenTK.Vector3d[] V1s, OpenTK.Vector3d[] V2s, double[] r)
        {
            System.Diagnostics.Debug.Assert(dof == 6);
            System.Diagnostics.Debug.Assert(nodeCnt == 3);
            double r1 = r[0];
            double r2 = r[1];
            double r3 = r[2];
            double[] Ns = new double[3];
            Ns[0] = 1.0 - r1 - r2;
            Ns[1] = r1;
            Ns[2] = r2;
            // partial r1,r2,r3
            double[] N1r = new double[3];
            double[] N2r = new double[3];
            double[] N3r = new double[3];
            N1r[0] = -1.0;
            N1r[1] = -1.0;
            N1r[2] = 0.0;
            N2r[0] = 1.0;
            N2r[1] = 0.0;
            N2r[2] = 0.0;
            N3r[0] = 0.0;
            N3r[1] = 1.0;
            N3r[2] = 0.0;
            double[][] Nrs = new double[3][] { N1r, N2r, N3r };

            OpenTK.Vector3d e1 = new OpenTK.Vector3d(1.0, 0.0, 0.0);
            OpenTK.Vector3d e2 = new OpenTK.Vector3d(0.0, 1.0, 0.0);
            OpenTK.Vector3d e3 = new OpenTK.Vector3d(0.0, 0.0, 1.0);
            OpenTK.Vector3d[] es = { e1, e2, e3 };

            OpenTK.Vector3d[][] coeffAis = new OpenTK.Vector3d[3][];
            for (int i = 0; i < 3; i++)
            {
                coeffAis[i] = new OpenTK.Vector3d[nodeCnt * dof];
                double deltai3 = i == 2 ? 1.0 : 0.0;
                for (int k = 0; k < nodeCnt; k++)
                {
                    for (int d = 0; d < 3; d++)
                    {
                        coeffAis[i][k * dof + d] = Nrs[k][i] * es[d];
                    }
                    coeffAis[i][k * dof + 3] = -(r3 / 2.0) * Nrs[k][i] * h * V2s[k] -
                        (1.0 / 2.0) * deltai3 * Ns[k] * h * V2s[k];
                    coeffAis[i][k * dof + 4] = (r3 / 2.0) * Nrs[k][i] * h * V1s[k] +
                        (1.0 / 2.0) * deltai3 * Ns[k] * h * V1s[k];
                    coeffAis[i][k * dof + 5] = new OpenTK.Vector3d();
                }
            }
            return coeffAis;
        }

        protected static double[,,] CalcMITCLinearStrainCoeff(
            int dof, int nodeCnt, OpenTK.Vector3d[] gs, OpenTK.Vector3d[][] coeffAis)
        {
            System.Diagnostics.Debug.Assert(dof == 6);
            System.Diagnostics.Debug.Assert(nodeCnt == 3);
            double[,,] smallBmn = new double[3, 3, nodeCnt * dof];
            for (int m = 0; m < 3; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    for (int k = 0; k < nodeCnt * dof; k++)
                    {
                        double bmnkValue = (1.0 / 2.0) * (
                            OpenTK.Vector3d.Dot(gs[m], coeffAis[n][k]) +
                            OpenTK.Vector3d.Dot(gs[n], coeffAis[m][k]));
                        smallBmn[m, n, k] = bmnkValue;
                    }
                }
            }
            return smallBmn;
        }

        protected static OpenTK.Vector3d[] CalcMITCLinearUCoeffVectors(
            int dof, int nodeCnt, double h, OpenTK.Vector3d[] Vns, OpenTK.Vector3d[] V1s, OpenTK.Vector3d[] V2s, double[] r)
        {
            System.Diagnostics.Debug.Assert(dof == 6);
            System.Diagnostics.Debug.Assert(nodeCnt == 3);
            double r1 = r[0];
            double r2 = r[1];
            double r3 = r[2];
            double[] Ns = new double[3];
            Ns[0] = 1.0 - r1 - r2;
            Ns[1] = r1;
            Ns[2] = r2;

            OpenTK.Vector3d e1 = new OpenTK.Vector3d(1.0, 0.0, 0.0);
            OpenTK.Vector3d e2 = new OpenTK.Vector3d(0.0, 1.0, 0.0);
            OpenTK.Vector3d e3 = new OpenTK.Vector3d(0.0, 0.0, 1.0);
            OpenTK.Vector3d[] es = { e1, e2, e3 };

            OpenTK.Vector3d[] coeffBs = new OpenTK.Vector3d[nodeCnt * dof];
            for (int k = 0; k < nodeCnt; k++)
            {
                for (int d = 0; d < 3; d++)
                {
                    coeffBs[k * dof + d] = Ns[k] * es[d];
                }
                coeffBs[k * dof + 3] = -(r3 / 2.0) * Ns[k] * h * V2s[k];
                coeffBs[k * dof + 4] = (r3 / 2.0) * Ns[k] * h * V1s[k];
                coeffBs[k * dof + 5] = new OpenTK.Vector3d();
            }
            return coeffBs;
        }

        public static IvyFEM.Lapack.DoubleMatrix CalcMITCLinearTransferMatrix(
            int dof, int nodeCnt, OpenTK.Vector3d[] Vns, OpenTK.Vector3d[] V1s, OpenTK.Vector3d[] V2s)
        {
            IvyFEM.Lapack.DoubleMatrix T = new Lapack.DoubleMatrix(nodeCnt * dof, nodeCnt * dof);
            for (int k = 0; k < nodeCnt; k++)
            {
                double[] doubleV1 = { V1s[k].X, V1s[k].Y, V1s[k].Z };
                double[] doubleV2 = { V2s[k].X, V2s[k].Y, V2s[k].Z };
                double[] doubleVn = { Vns[k].X, Vns[k].Y, Vns[k].Z };
                IvyFEM.Lapack.DoubleMatrix T3 = new IvyFEM.Lapack.DoubleMatrix(3, 3);
                for (int j = 0; j < 3; j++)
                {
                    T3[0, j] = doubleV1[j];
                    T3[1, j] = doubleV2[j];
                    T3[2, j] = doubleVn[j];
                }
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        T[i + 6 * k, j + 6 * k] = i == j ? 1.0 : 0.0;
                    }
                }
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        T[i + 6 * k + 3, j + 6 * k + 3] = T3[i, j];
                    }
                }
            }
            return T;
        }

        public static IvyFEM.Lapack.DoubleMatrix CalcMITCLinearPlateKe(
            TriangleFE d1TriFE, double h, OpenTK.Vector3d[] xPts,
            IvyFEM.Lapack.DoubleMatrix Cb, IvyFEM.Lapack.DoubleMatrix Cs, double kappa)
        {
            int nodeCnt = xPts.Length;
            System.Diagnostics.Debug.Assert(nodeCnt == 3);

            OpenTK.Vector3d[] Vns;
            OpenTK.Vector3d[] V1s;
            OpenTK.Vector3d[] V2s;
            CalcMITCLinearDirectorVectors(xPts, out Vns, out V1s, out V2s);

            int uDof = 3;
            int aDof = 1;
            int bDof = 1;
            int zDof = 1;
            int dof = uDof + aDof + bDof + zDof;
            int aOffset = uDof;
            int bOffset = aOffset + aDof;
            int zOffset = bOffset + bDof;
            System.Diagnostics.Debug.Assert(dof == 6);

            IvyFEM.Lapack.DoubleMatrix localKe = new IvyFEM.Lapack.DoubleMatrix(nodeCnt * dof, nodeCnt * dof);

            //------------------------------------
            // sampling for MITC
            double[][] sampleRs = new double[3][];
            sampleRs[0] = new double[3] { 1.0 / 2.0, 0.0, 0.0 };
            sampleRs[1] = new double[3] { 0.0, 1.0 / 2.0, 0.0 };
            sampleRs[2] = new double[3] { 1.0 / 2.0, 1.0 / 2.0, 0.0 };
            double[][,,] sampleSmallBmns = new double[3][,,];
            for (int iSample = 0; iSample < 3; iSample++)
            {
                double[] sampleR = sampleRs[iSample];
                OpenTK.Vector3d[] gs = CalcMITCLinearCovariantBasisVectors(xPts, h, Vns, sampleR);
                OpenTK.Vector3d[][] coeffAis = CalcMITCLinearUiCoeffVectors(
                    dof, nodeCnt, h, Vns, V1s, V2s, sampleR);
                double[,,] smallBmn = CalcMITCLinearStrainCoeff(dof, nodeCnt, gs, coeffAis);
                sampleSmallBmns[iSample] = smallBmn;
            }
            double[] sampleSmallB13Pt1 = new double[nodeCnt * dof];
            double[] sampleSmallB23Pt2 = new double[nodeCnt * dof];
            double[] sampleSmallB13Pt3 = new double[nodeCnt * dof];
            double[] sampleSmallB23Pt3 = new double[nodeCnt * dof];
            for (int k = 0; k < nodeCnt * dof; k++)
            {
                sampleSmallB13Pt1[k] = sampleSmallBmns[0][0, 2, k];
                sampleSmallB23Pt2[k] = sampleSmallBmns[1][1, 2, k];
                sampleSmallB13Pt3[k] = sampleSmallBmns[2][0, 2, k];
                sampleSmallB23Pt3[k] = sampleSmallBmns[2][1, 2, k];
            }
            //------------------------------------

            IntegrationPoints ipZ = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point2);
            IntegrationPoints ipXY = TriangleFE.GetIntegrationPoints(TriangleIntegrationPointCount.Point4);
            System.Diagnostics.Debug.Assert(ipZ.Ls.Length == 2);
            System.Diagnostics.Debug.Assert(ipXY.Ls.Length == 4);
            for (int ipZPt = 0; ipZPt < ipZ.PointCount; ipZPt++)
            {
                double[] LZPt = ipZ.Ls[ipZPt];
                double lineLenZPt = h;
                double weightZPt = ipZ.Weights[ipZPt];
                double detJWeightZPt = (lineLenZPt / 2.0) * weightZPt;

                double r3 = 2.0 * LZPt[1] - 1.0;

                for (int ipXYPt = 0; ipXYPt < ipXY.PointCount; ipXYPt++)
                {
                    double[] L = ipXY.Ls[ipXYPt];
                    double detJ = d1TriFE.GetDetJacobian(L);
                    double weight = ipXY.Weights[ipXYPt];
                    double detJWeight = (1.0 / 2.0) * weight * detJ;

                    double r1 = L[1];
                    double r2 = L[2];
                    double[] r = new double[3] { r1, r2, r3 };

                    // from displacements
                    OpenTK.Vector3d[] gs = CalcMITCLinearCovariantBasisVectors(xPts, h, Vns, r);
                    OpenTK.Vector3d[] hs = CalcMITCLinearContravariantBasisVectors(gs);
                    OpenTK.Vector3d[] lVecs = CalcMITCLinearLocalBasicVectors(gs);
                    OpenTK.Vector3d[][] coeffAis = CalcMITCLinearUiCoeffVectors(dof, nodeCnt, h, Vns, V1s, V2s, r);
                    double[,,] smallBmn = CalcMITCLinearStrainCoeff(dof, nodeCnt, gs, coeffAis);

                    //-------------------------------------
                    // MITC
                    for (int k = 0; k < nodeCnt * dof; k++)
                    {
                        double coeffC = sampleSmallB13Pt3[k] - sampleSmallB13Pt1[k] -
                            (sampleSmallB23Pt3[k] - sampleSmallB23Pt2[k]);
                        // rt
                        double hSmallB13 = sampleSmallB13Pt1[k] + coeffC * r2;
                        double hSmallB23 = sampleSmallB23Pt2[k] - coeffC * r1;

                        // rt
                        smallBmn[0, 2, k] = hSmallB13;
                        smallBmn[2, 0, k] = hSmallB13;
                        // st
                        smallBmn[1, 2, k] = hSmallB23;
                        smallBmn[2, 1, k] = hSmallB23;
                    }
                    //-------------------------------------

                    double[,,] Bij = new double[3, 3, nodeCnt * dof];
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            for (int k = 0; k < nodeCnt * dof; k++)
                            {
                                double bijkValue = 0.0;
                                for (int m = 0; m < 3; m++)
                                {
                                    for (int n = 0; n < 3; n++)
                                    {
                                        bijkValue += OpenTK.Vector3d.Dot(lVecs[i], hs[m]) *
                                            OpenTK.Vector3d.Dot(lVecs[j], hs[n]) * smallBmn[m, n, k];
                                    }
                                }
                                Bij[i, j, k] = bijkValue;
                            }
                        }
                    }

                    IvyFEM.Lapack.DoubleMatrix Bb = new IvyFEM.Lapack.DoubleMatrix(3, nodeCnt * dof);
                    IvyFEM.Lapack.DoubleMatrix Bs = new IvyFEM.Lapack.DoubleMatrix(2, nodeCnt * dof);
                    for (int k = 0; k < nodeCnt * dof; k++)
                    {
                        // ε11
                        Bb[0, k] = Bij[0, 0, k];
                        // ε22
                        Bb[1, k] = Bij[1, 1, k];
                        // 2ε12
                        Bb[2, k] = 2.0 * Bij[0, 1, k];
                       
                        // 2ε13
                        Bs[0, k] = 2.0 * Bij[0, 2, k];
                        // 2ε23
                        Bs[1, k] = 2.0 * Bij[1, 2, k];
                    }
                    var transBb = IvyFEM.Lapack.DoubleMatrix.Transpose(Bb);
                    var transBs = IvyFEM.Lapack.DoubleMatrix.Transpose(Bs);
                    var Kbeip = transBb * Cb * Bb;
                    var Kseip = transBs * Cs * Bs;
                    Kseip = IvyFEM.Lapack.DoubleMatrix.Scal(Kseip, kappa);

                    for (int i = 0; i < nodeCnt * dof; i++)
                    {
                        for (int j = 0; j < nodeCnt * dof; j++)
                        {
                            localKe[i, j] += detJWeightZPt * detJWeight * (Kbeip[i, j] + Kseip[i, j]);
                        }
                    }
                }
            }

            //------------------------------------
            for (int iNode = 0; iNode < 3; iNode++)
            {
                // fictitious stiffness value
                double f = 0.0;
                for (int k = 0; k < (uDof + aDof + bDof); k++)
                {
                    double diagVal = localKe[iNode * dof + k, iNode * dof + k];
                    System.Diagnostics.Debug.Assert(diagVal >= IvyFEM.Constants.PrecisionLowerLimit);
                    if (diagVal > f)
                    {
                        f = diagVal;
                    }
                }
                f *= 1.0e-3;
                localKe[iNode * dof + zOffset, iNode * dof + zOffset] = f;
            }

            //-------------------------------------------
            // global
            var T = CalcMITCLinearTransferMatrix(dof, nodeCnt, Vns, V1s, V2s);
            var transT = IvyFEM.Lapack.DoubleMatrix.Transpose(T);
            var Ke = transT * localKe * T;

            return Ke;
        }

        public static IvyFEM.Lapack.DoubleMatrix CalcMITCLinearPlateMe(
            TriangleFE dTriFE, double h, OpenTK.Vector3d[] xPts, double rho)
        {
            int nodeCnt = xPts.Length;
            System.Diagnostics.Debug.Assert(nodeCnt == 3);

            OpenTK.Vector3d[] Vns;
            OpenTK.Vector3d[] V1s;
            OpenTK.Vector3d[] V2s;
            CalcMITCLinearDirectorVectors(xPts, out Vns, out V1s, out V2s);

            int uDof = 3;
            int aDof = 1;
            int bDof = 1;
            int zDof = 1;
            int dof = uDof + aDof + bDof + zDof;
            int aOffset = uDof;
            int bOffset = aOffset + aDof;
            int zOffset = bOffset + bDof;
            System.Diagnostics.Debug.Assert(dof == 6);

            IvyFEM.Lapack.DoubleMatrix localMe = new IvyFEM.Lapack.DoubleMatrix(dof * 3, dof * 3);

            IntegrationPoints ipZ = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point2);
            IntegrationPoints ipXY = TriangleFE.GetIntegrationPoints(TriangleIntegrationPointCount.Point4);
            System.Diagnostics.Debug.Assert(ipZ.Ls.Length == 2);
            System.Diagnostics.Debug.Assert(ipXY.Ls.Length == 4);
            for (int ipZPt = 0; ipZPt < ipZ.PointCount; ipZPt++)
            {
                double[] LZPt = ipZ.Ls[ipZPt];
                double lineLenZPt = h;
                double weightZPt = ipZ.Weights[ipZPt];
                double detJWeightZPt = (lineLenZPt / 2.0) * weightZPt;

                double r3 = 2.0 * LZPt[1] - 1.0;

                for (int ipXYPt = 0; ipXYPt < ipXY.PointCount; ipXYPt++)
                {
                    double[] L = ipXY.Ls[ipXYPt];
                    double detJ = dTriFE.GetDetJacobian(L);
                    double weight = ipXY.Weights[ipXYPt];
                    double detJWeight = (1.0 / 2.0) * weight * detJ;

                    double r1 = L[1];
                    double r2 = L[2];
                    double[] r = new double[3] { r1, r2, r3 };

                    OpenTK.Vector3d[] coeffBs = CalcMITCLinearUCoeffVectors(dof, nodeCnt, h, Vns, V1s, V2s, r);
                    IvyFEM.Lapack.DoubleMatrix B = new IvyFEM.Lapack.DoubleMatrix(3, nodeCnt * dof);
                    for (int k = 0; k < nodeCnt * dof; k++)
                    {
                        B[0, k] = coeffBs[k].X;
                        B[1, k] = coeffBs[k].Y;
                        B[2, k] = coeffBs[k].Z;
                    }
                    var transB = IvyFEM.Lapack.DoubleMatrix.Transpose(B);
                    var Meip = transB * B;
                    Meip = IvyFEM.Lapack.DoubleMatrix.Scal(Meip, rho);

                    for (int i = 0; i < nodeCnt * dof; i++)
                    {
                        for (int j = 0; j < nodeCnt * dof; j++)
                        {
                            localMe[i, j] += detJWeightZPt * detJWeight * Meip[i, j];
                        }
                    }
                }
            }

            //-------------------------------------------
            // global
            var T = CalcMITCLinearTransferMatrix(dof, nodeCnt, Vns, V1s, V2s);
            var transT = IvyFEM.Lapack.DoubleMatrix.Transpose(T);
            var Me = transT * localMe * T;

            return Me;
        }
    }
}
