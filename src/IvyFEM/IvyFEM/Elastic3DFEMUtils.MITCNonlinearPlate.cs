using OpenTK.Audio.OpenAL;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Animation;

namespace IvyFEM
{
    public partial class Elastic3DFEMUtils
    {
        protected static void CalcMITCNonlinearInitialDirectorVectors(
            OpenTK.Vector3d[] xPt0s, out OpenTK.Vector3d[] Vn0s, out OpenTK.Vector3d[] V10s, out OpenTK.Vector3d[] V20s)
        {
            int nodeCnt = xPt0s.Length;
            System.Diagnostics.Debug.Assert(nodeCnt == 3);
            OpenTK.Vector3d normal = CadUtils3D.TriNormal(xPt0s[0], xPt0s[1], xPt0s[2]);
            OpenTK.Vector3d e1 = new OpenTK.Vector3d(1.0, 0.0, 0.0);
            OpenTK.Vector3d e2 = new OpenTK.Vector3d(0.0, 1.0, 0.0);
            OpenTK.Vector3d e3 = new OpenTK.Vector3d(0.0, 0.0, 1.0);
            // director vectors etc
            Vn0s = new OpenTK.Vector3d[nodeCnt];
            V10s = new OpenTK.Vector3d[nodeCnt];
            V20s = new OpenTK.Vector3d[nodeCnt];
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
                Vn0s[i] = vn;
                V10s[i] = v1;
                V20s[i] = v2;
            }
        }

        protected static void UpdateMITCNonlinearDirectorVectors(
            int nodeCnt, double[] alphas, double[] betas,
            OpenTK.Vector3d[] prevVns, OpenTK.Vector3d[] prevV1s, OpenTK.Vector3d[] prevV2s,
            out OpenTK.Vector3d[] Vns, out OpenTK.Vector3d[] V1s, out OpenTK.Vector3d[] V2s)
        {
            Vns = new OpenTK.Vector3d[nodeCnt];
            V1s = new OpenTK.Vector3d[nodeCnt];
            V2s = new OpenTK.Vector3d[nodeCnt];
            for (int iNode = 0; iNode < nodeCnt; iNode++)
            {
                double alpha = alphas[iNode];
                double beta = betas[iNode];
                double theta = Math.Sqrt(alpha * alpha + beta * beta);
                var S = new IvyFEM.Lapack.DoubleMatrix(3, 3);
                S[0, 0] = 0.0;
                S[0, 1] = 0.0;
                S[0, 2] = beta;
                S[1, 0] = 0.0;
                S[1, 1] = 0.0;
                S[1, 2] = -alpha;
                S[2, 0] = -beta;
                S[2, 1] = alpha;
                S[2, 2] = 0.0;
                var I = new IvyFEM.Lapack.DoubleMatrix(3, 3);
                I.Identity();
                var work1 = S;
                var work2 = IvyFEM.Lapack.DoubleMatrix.Scal(S * S, (1.0 / 2.0));
                var R = new IvyFEM.Lapack.DoubleMatrix(3, 3);
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        R[i, j] = I[i, j] + work1[i, j] + work2[i, j];
                    }
                }

                OpenTK.Vector3d prevVn = prevVns[iNode];
                OpenTK.Vector3d prevV1 = prevV1s[iNode];
                OpenTK.Vector3d prevV2 = prevV2s[iNode];
                double[] prevVnVec = { prevVn.X, prevVn.Y, prevVn.Z };
                double[] prevV1Vec = { prevV1.X, prevV1.Y, prevV1.Z };
                double[] prevV2Vec = { prevV2.X, prevV2.Y, prevV2.Z };
                double[] VnVec = R * prevVnVec;
                double[] V1Vec = R * prevV1Vec;
                double[] V2Vec = R * prevV2Vec;
                OpenTK.Vector3d Vn = new OpenTK.Vector3d(VnVec[0], VnVec[1], VnVec[2]);
                OpenTK.Vector3d V1 = new OpenTK.Vector3d(V1Vec[0], V1Vec[1], V1Vec[2]);
                OpenTK.Vector3d V2 = new OpenTK.Vector3d(V2Vec[0], V2Vec[1], V2Vec[2]);
                Vns[iNode] = Vn;
                V1s[iNode] = V1;
                V2s[iNode] = V2;
            }
        }

        protected static OpenTK.Vector3d[] CalcMITCNonlinearG0s(
            OpenTK.Vector3d[] xPt0s, double h, double[] r)
        {
            int nodeCnt = xPt0s.Length;
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

            // initial director vectors
            OpenTK.Vector3d[] Vn0s;
            OpenTK.Vector3d[] V10s;
            OpenTK.Vector3d[] V20s;
            CalcMITCNonlinearInitialDirectorVectors(
                xPt0s, out Vn0s, out V10s, out V20s);

            // g0
            OpenTK.Vector3d[] g0s = new OpenTK.Vector3d[3];
            for (int i = 0; i < 3; i++)
            {
                double deltai3 = i == 2 ? 1.0 : 0.0;
                g0s[i] = new OpenTK.Vector3d();
                for (int k = 0; k < nodeCnt; k++)
                {
                    g0s[i] += Nrs[k][i] * xPt0s[k] + (r3 / 2.0) * Nrs[k][i] * h * Vn0s[k] +
                        (1.0 / 2.0) * deltai3 * Ns[k] * h * Vn0s[k];
                }
            }
            return g0s;
        }

        protected static OpenTK.Vector3d[] CalcMITCNonlinearCurGs(
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

            // gt
            OpenTK.Vector3d[] curGs = new OpenTK.Vector3d[3];
            for (int i = 0; i < 3; i++)
            {
                double deltai3 = i == 2 ? 1.0 : 0.0;
                curGs[i] = new OpenTK.Vector3d();
                for (int k = 0; k < nodeCnt; k++)
                {
                    curGs[i] += Nrs[k][i] * xPts[k] + (r3 / 2.0) * Nrs[k][i] * h * Vns[k] +
                        (1.0 / 2.0) * deltai3 * Ns[k] * h * Vns[k];
                }
            }
            return curGs;
        }

        protected static IvyFEM.Lapack.DoubleMatrix CalcMITCNonliearCurStrain(
            OpenTK.Vector3d[] g0s, OpenTK.Vector3d[] curGs)
        {
            IvyFEM.Lapack.DoubleMatrix E = new IvyFEM.Lapack.DoubleMatrix(3, 3);
            for (int ig = 0; ig < 3; ig++)
            {
                for (int ih = 0; ih < 3; ih++)
                {
                    E[ig, ih] = (1.0 / 2.0) *( 
                        OpenTK.Vector3d.Dot(curGs[ig], curGs[ih]) -
                        OpenTK.Vector3d.Dot(g0s[ig], g0s[ih]));
                }
            }
            return E;
        }

        protected static void CalcMITCNonlinearUrgCoeffVectors(
            int dof, int nodeCnt, 
            double h, OpenTK.Vector3d[] Vns, OpenTK.Vector3d[] V1s, OpenTK.Vector3d[] V2s, double[] r,
            out OpenTK.Vector3d[][] coeffA1s, out OpenTK.Vector3d[][,] coeffA2s)
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

            coeffA1s = new OpenTK.Vector3d[3][];
            for (int ig = 0; ig < 3; ig++)
            {
                coeffA1s[ig] = new OpenTK.Vector3d[nodeCnt * dof];
                double deltaig3 = ig == 2 ? 1.0 : 0.0;
                for (int k = 0; k < nodeCnt; k++)
                {
                    for (int d = 0; d < 3; d++)
                    {
                        coeffA1s[ig][k * dof + d] = Nrs[k][ig] * es[d];
                    }
                    coeffA1s[ig][k * dof + 3] = -(r3 / 2.0) * Nrs[k][ig] * h * V2s[k] -
                        (1.0 / 2.0) * deltaig3 * Ns[k] * h * V2s[k];
                    coeffA1s[ig][k * dof + 4] = (r3 / 2.0) * Nrs[k][ig] * h * V1s[k] +
                        (1.0 / 2.0) * deltaig3 * Ns[k] * h * V1s[k];
                    coeffA1s[ig][k * dof + 5] = new OpenTK.Vector3d();
                }
            }

            coeffA2s = new OpenTK.Vector3d[3][,];
            for (int ig = 0; ig < 3; ig++)
            {
                coeffA2s[ig] = new OpenTK.Vector3d[nodeCnt * dof, nodeCnt * dof];
                double deltaig3 = ig == 2 ? 1.0 : 0.0;
                // only diagonal
                for (int k = 0; k < nodeCnt; k++)
                {
                    for (int d = 0; d < 3; d++)
                    {
                        coeffA2s[ig][k * dof + d, k * dof + d] = new OpenTK.Vector3d();
                    }
                    coeffA2s[ig][k * dof + 3, k * dof + 3] =
                        -(r3 / 4.0) * Nrs[k][ig] * h * Vns[k] - (1.0 / 4.0) * deltaig3 * Ns[k] * h * Vns[k];
                    coeffA2s[ig][k * dof + 4, k * dof + 4] =
                        -(r3 / 4.0) * Nrs[k][ig] * h * Vns[k] - (1.0 / 4.0) * deltaig3 * Ns[k] * h * Vns[k];
                    coeffA2s[ig][k * dof + 5, k * dof + 5] = new OpenTK.Vector3d();
                }
            }
        }

        protected static void CalcMITCNonlinearStrainCoeff(
            int dof, int nodeCnt,
            OpenTK.Vector3d[] curGs, OpenTK.Vector3d[][] coeffA1s, OpenTK.Vector3d[][,] coeffA2s,
            out double[,][] B1s, out IvyFEM.Lapack.DoubleMatrix[,] B2s)
        {
            System.Diagnostics.Debug.Assert(dof == 6);
            System.Diagnostics.Debug.Assert(nodeCnt == 3);

            B1s = new double[3, 3][];
            for (int ig = 0; ig < 3; ig++)
            {
                for (int ih = 0; ih < 3; ih++)
                {
                    B1s[ig, ih] = new double[nodeCnt* dof];
                    for (int k = 0; k < nodeCnt * dof; k++)
                    {
                        B1s[ig, ih][k] = (1.0 / 2.0) * (
                            OpenTK.Vector3d.Dot(curGs[ih], coeffA1s[ig][k]) +
                            OpenTK.Vector3d.Dot(curGs[ig], coeffA1s[ih][k]));
                    }
                }
            }

            B2s = new IvyFEM.Lapack.DoubleMatrix[3, 3];
            for (int ig = 0; ig < 3; ig++)
            {
                for (int ih = 0; ih < 3; ih++)
                {
                    B2s[ig, ih] = new IvyFEM.Lapack.DoubleMatrix(nodeCnt * dof, nodeCnt * dof);
                    for (int k = 0; k < nodeCnt * dof; k++)
                    {
                        B2s[ig, ih][k, k] +=
                            OpenTK.Vector3d.Dot(curGs[ih], coeffA2s[ig][k, k]) +
                            OpenTK.Vector3d.Dot(curGs[ig], coeffA2s[ih][k, k]);
                    }
                    for (int m = 0; m < nodeCnt * dof; m++)
                    {
                        for (int n = 0; n < nodeCnt * dof; n++)
                        {
                            B2s[ig, ih][m, n] +=
                                OpenTK.Vector3d.Dot(coeffA1s[ig][m], coeffA1s[ih][n]);  
                        }
                    }
                }
            }
        }

        protected static OpenTK.Vector3d[] CalcMITCNonlinearContravariantBasisVectors(OpenTK.Vector3d[] gs)
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

        // St.Venant-Kirchhoff
        protected static double[,,,] CalcMITCNonlinearStVenantC4Tensor(
            double lambda, double mu, OpenTK.Vector3d[] h0s)
        {
            double[,,,] C4 = new double[3, 3, 3, 3];
            for (int ig = 0; ig < 3; ig++)
            {
                for (int ih = 0; ih < 3; ih++)
                {
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            OpenTK.Vector3d h0i = h0s[i];
                            OpenTK.Vector3d h0j = h0s[j];
                            OpenTK.Vector3d h0g = h0s[ig];
                            OpenTK.Vector3d h0h = h0s[ih];
                            double h0ij = OpenTK.Vector3d.Dot(h0i, h0j);
                            double h0gh = OpenTK.Vector3d.Dot(h0g, h0h);
                            double h0ig = OpenTK.Vector3d.Dot(h0i, h0g);
                            double h0jh = OpenTK.Vector3d.Dot(h0j, h0h);
                            C4[ig, ih, i, j] = lambda * h0ij * h0gh + 2.0 * mu * h0ig * h0jh;
                        }
                    }
                }
            }
            return C4;
        }

        protected static IvyFEM.Lapack.DoubleMatrix CalcMITCNonlinearStress(double[,,,] C4, IvyFEM.Lapack.DoubleMatrix E)
        {
            IvyFEM.Lapack.DoubleMatrix S = new IvyFEM.Lapack.DoubleMatrix(3, 3);
            for (int ig = 0; ig < 3; ig++)
            {
                for (int ih = 0; ih < 3; ih++)
                {
                    double tmp = 0.0;
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            tmp += C4[ig, ih, i, j] * E[i, j];
                        }
                    }
                    S[ig, ih] = tmp;
                }
            }
            return S;
        }

        protected static void CalcMITCNonlinearUCoeffVectors(
            int dof, int nodeCnt,
            double h, OpenTK.Vector3d[] Vns, OpenTK.Vector3d[] V1s, OpenTK.Vector3d[] V2s, double[] r,
            out OpenTK.Vector3d[] coeffB1s, out OpenTK.Vector3d[,] coeffB2s)
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

            coeffB1s = new OpenTK.Vector3d[nodeCnt * dof];
            for (int k = 0; k < nodeCnt; k++)
            {
                for (int d = 0; d < 3; d++)
                {
                    coeffB1s[k * dof + d] = Ns[k] * es[d];
                }
                coeffB1s[k * dof + 3] = -(r3 / 2.0) * Ns[k] * h * V2s[k];
                coeffB1s[k * dof + 4] = (r3 / 2.0) * Ns[k] * h * V1s[k];
                coeffB1s[k * dof + 5] = new OpenTK.Vector3d();
            }

            coeffB2s = new OpenTK.Vector3d[nodeCnt * dof, nodeCnt * dof];
            {
                // only diagonal
                for (int k = 0; k < nodeCnt; k++)
                {
                    for (int d = 0; d < 3; d++)
                    {
                        coeffB2s[k * dof + d, k * dof + d] = new OpenTK.Vector3d();
                    }
                    coeffB2s[k * dof + 3, k * dof + 3] = -(r3 / 4.0) * Ns[k] * h * Vns[k];
                    coeffB2s[k * dof + 4, k * dof + 4] = -(r3 / 4.0) * Ns[k] * h * Vns[k];
                    coeffB2s[k * dof + 5, k * dof + 5] = new OpenTK.Vector3d();
                }
            }
        }

        public static IvyFEM.Lapack.DoubleMatrix CalcMITCNonlinearTransferMatrix(
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

        public static void InitMITCNonlinearPlateNodeValues(int coId, FieldValue uFV)
        {
            //-------------------------
            // displacement
            int uDof = 3;
            System.Diagnostics.Debug.Assert(uFV.Dof == uDof);
            System.Diagnostics.Debug.Assert(uFV.NodeType == FieldValueNodeType.Node);
            double[] uValues = uFV.GetDoubleValues(FieldDerivativeType.Value);
            for (int iDof = 0; iDof < uDof; iDof++)
            {
                uValues[coId * uDof + iDof] = 0.0;
            }
            //-------------------------
        }

        public static void InitMITCNonlinearPlateElementNodeValues(
            uint feId, OpenTK.Vector3d[] xPt0s, FieldValue vnFV, FieldValue v1FV, FieldValue v2FV)
        {
            int nodeCnt = xPt0s.Length;

            //----------------------------
            // director vectors
            OpenTK.Vector3d[] V0ns;
            OpenTK.Vector3d[] V01s;
            OpenTK.Vector3d[] V02s;
            CalcMITCNonlinearInitialDirectorVectors(xPt0s, out V0ns, out V01s, out V02s);

            int vDof = 3;
            System.Diagnostics.Debug.Assert(vnFV.Dof == vDof);
            System.Diagnostics.Debug.Assert(v1FV.Dof == vDof);
            System.Diagnostics.Debug.Assert(v2FV.Dof == vDof);
            System.Diagnostics.Debug.Assert(vnFV.NodeType == FieldValueNodeType.ElementNode);
            System.Diagnostics.Debug.Assert(v1FV.NodeType == FieldValueNodeType.ElementNode);
            System.Diagnostics.Debug.Assert(v2FV.NodeType == FieldValueNodeType.ElementNode);
            double[] vnValues = vnFV.GetDoubleValues(FieldDerivativeType.Value);
            double[] v1Values = v1FV.GetDoubleValues(FieldDerivativeType.Value);
            double[] v2Values = v2FV.GetDoubleValues(FieldDerivativeType.Value);
            System.Diagnostics.Debug.Assert(feId > 0);
            for (int iNode = 0; iNode < nodeCnt; iNode++)
            {
                int id = (int)((feId - 1) * nodeCnt + iNode);
                OpenTK.Vector3d Vn = V0ns[iNode];
                OpenTK.Vector3d V1 = V01s[iNode];
                OpenTK.Vector3d V2 = V02s[iNode];
                double[] doubleVn = { Vn.X, Vn.Y, Vn.Z };
                double[] doubleV1 = { V1.X, V1.Y, V1.Z };
                double[] doubleV2 = { V2.X, V2.Y, V2.Z };
                for (int iDof = 0; iDof < vDof; iDof++)
                {
                    vnValues[id * vDof + iDof] = doubleVn[iDof];
                    v1Values[id * vDof + iDof] = doubleV1[iDof];
                    v2Values[id * vDof + iDof] = doubleV2[iDof];
                }
            }
            //----------------------------
        }

        public static void UpdateMITCNonlinearPlateNodeValues(int coId, double[] curNodeUg, FieldValue uFV)
        {
            // uは変換しなくてよい
            OpenTK.Vector3d u = new OpenTK.Vector3d(curNodeUg[0], curNodeUg[1], curNodeUg[2]);

            //-------------------------
            // displacement
            int uDof = 3;
            System.Diagnostics.Debug.Assert(uFV.Dof == uDof);
            System.Diagnostics.Debug.Assert(uFV.NodeType == FieldValueNodeType.Node);
            double[] uValues = uFV.GetDoubleValues(FieldDerivativeType.Value);
            double[] doubleU = { u.X, u.Y, u.Z };
            for (int iDof = 0; iDof < uDof; iDof++)
            {
                uValues[coId * uDof + iDof] += doubleU[iDof];
            }
            //-------------------------
        }

        public static void UpdateMITCNonlinearPlateElementNodeValues(
            uint feId, int nodeCnt, double[] curUg, FieldValue vnFV, FieldValue v1FV, FieldValue v2FV)
        {
            //----------------------------
            // director vectors
            int vDof = 3;
            System.Diagnostics.Debug.Assert(vnFV.Dof == vDof);
            System.Diagnostics.Debug.Assert(v1FV.Dof == vDof);
            System.Diagnostics.Debug.Assert(v2FV.Dof == vDof);
            System.Diagnostics.Debug.Assert(vnFV.NodeType == FieldValueNodeType.ElementNode);
            System.Diagnostics.Debug.Assert(v1FV.NodeType == FieldValueNodeType.ElementNode);
            System.Diagnostics.Debug.Assert(v2FV.NodeType == FieldValueNodeType.ElementNode);
            double[] vnValues = vnFV.GetDoubleValues(FieldDerivativeType.Value);
            double[] v1Values = v1FV.GetDoubleValues(FieldDerivativeType.Value);
            double[] v2Values = v2FV.GetDoubleValues(FieldDerivativeType.Value);

            OpenTK.Vector3d[] prevVns = new OpenTK.Vector3d[nodeCnt];
            OpenTK.Vector3d[] prevV1s = new OpenTK.Vector3d[nodeCnt];
            OpenTK.Vector3d[] prevV2s = new OpenTK.Vector3d[nodeCnt];
            System.Diagnostics.Debug.Assert(feId > 0);
            for (int iNode = 0; iNode < nodeCnt; iNode++)
            {
                int id = (int)((feId - 1) * nodeCnt + iNode);
                prevVns[iNode] = new OpenTK.Vector3d(
                    vnValues[id * vDof], vnValues[id * vDof + 1], vnValues[id * vDof + 2]);
                prevV1s[iNode] = new OpenTK.Vector3d(
                    v1Values[id * vDof], v1Values[id * vDof + 1], v1Values[id * vDof + 2]);
                prevV2s[iNode] = new OpenTK.Vector3d(
                    v2Values[id * vDof], v2Values[id * vDof + 1], v2Values[id * vDof + 2]);
            }

            int uDof = 3;
            int aDof = 1;
            int bDof = 1;
            int zDof = 1;
            int dof = uDof + aDof + bDof + zDof;
            int aOffset = uDof;
            int bOffset = aOffset + aDof;
            int zOffset = bOffset + bDof;
            System.Diagnostics.Debug.Assert(dof == 6);

            // transfer matrix
            var T = CalcMITCNonlinearTransferMatrix(dof, nodeCnt, prevVns, prevV1s, prevV2s);
            var transT = IvyFEM.Lapack.DoubleMatrix.Transpose(T);
            // to local
            double[] curUl = T * curUg;
            OpenTK.Vector3d[] us = new OpenTK.Vector3d[nodeCnt];
            double[] alphas = new double[nodeCnt];
            double[] betas = new double[nodeCnt];
            for (int iNode = 0; iNode < nodeCnt; iNode++)
            {
                us[iNode] = new OpenTK.Vector3d(curUl[iNode * dof], curUl[iNode * dof + 1], curUl[iNode * dof + 2]);
                alphas[iNode] = curUl[iNode * dof + aOffset];
                betas[iNode] = curUl[iNode * dof + bOffset];
            }

            OpenTK.Vector3d[] Vns;
            OpenTK.Vector3d[] V1s;
            OpenTK.Vector3d[] V2s;
            UpdateMITCNonlinearDirectorVectors(
                nodeCnt, alphas, betas, prevVns, prevV1s, prevV2s, out Vns, out V1s, out V2s);

            for (int iNode = 0; iNode < nodeCnt; iNode++)
            {
                int id = (int)((feId - 1) * nodeCnt + iNode);
                OpenTK.Vector3d Vn = Vns[iNode];
                OpenTK.Vector3d V1 = V1s[iNode];
                OpenTK.Vector3d V2 = V2s[iNode];
                double[] doubleVn = { Vn.X, Vn.Y, Vn.Z };
                double[] doubleV1 = { V1.X, V1.Y, V1.Z };
                double[] doubleV2 = { V2.X, V2.Y, V2.Z };
                for (int iDof = 0; iDof < vDof; iDof++)
                {
                    vnValues[id * vDof + iDof] = doubleVn[iDof];
                    v1Values[id * vDof + iDof] = doubleV1[iDof];
                    v2Values[id * vDof + iDof] = doubleV2[iDof];
                }
            }
            //----------------------------
        }

        public static void CalcMITCStVenantPateKe(
            double lambda, double mu,
            uint feId, TriangleFE dTriFE, IList<int> dCoIds, double h, OpenTK.Vector3d[] xPt0s,
            FieldValue uFV, FieldValue vnFV, FieldValue v1FV, FieldValue v2FV,
            out double[] Qe, out IvyFEM.Lapack.DoubleMatrix Ke)
        {
            int nodeCnt = xPt0s.Length;
            System.Diagnostics.Debug.Assert(nodeCnt == 3);

            OpenTK.Vector3d[] xPts = new OpenTK.Vector3d[nodeCnt];
            for (int iNode = 0; iNode < nodeCnt; iNode++)
            {
                int coId = dCoIds[iNode];
                double[] prevU = uFV.GetDoubleValue(coId, FieldDerivativeType.Value);
                xPts[iNode] = xPt0s[iNode] + new OpenTK.Vector3d(prevU[0], prevU[1], prevU[2]);
            }
            OpenTK.Vector3d[] Vns = new OpenTK.Vector3d[nodeCnt];
            OpenTK.Vector3d[] V1s = new OpenTK.Vector3d[nodeCnt];
            OpenTK.Vector3d[] V2s = new OpenTK.Vector3d[nodeCnt];
            System.Diagnostics.Debug.Assert(feId > 0);
            for (int iNode = 0; iNode < nodeCnt; iNode++)
            {
                int id = (int)((feId - 1) * nodeCnt + iNode);
                double[] doubleVn = vnFV.GetDoubleValue(id, FieldDerivativeType.Value);
                double[] doubleV1 = v1FV.GetDoubleValue(id, FieldDerivativeType.Value);
                double[] doubleV2 = v2FV.GetDoubleValue(id, FieldDerivativeType.Value);
                Vns[iNode] = new OpenTK.Vector3d(doubleVn[0], doubleVn[1], doubleVn[2]);
                V1s[iNode] = new OpenTK.Vector3d(doubleV1[0], doubleV1[1], doubleV1[2]);
                V2s[iNode] = new OpenTK.Vector3d(doubleV2[0], doubleV2[1], doubleV2[2]);
            }

            int uDof = 3;
            int aDof = 1;
            int bDof = 1;
            int zDof = 1;
            int dof = uDof + aDof + bDof + zDof;
            int aOffset = uDof;
            int bOffset = aOffset + aDof;
            int zOffset = bOffset + bDof;
            System.Diagnostics.Debug.Assert(dof == 6);

            // transfer matrix
            var T = CalcMITCNonlinearTransferMatrix(dof, nodeCnt, Vns, V1s, V2s);
            var transT = IvyFEM.Lapack.DoubleMatrix.Transpose(T);

            double[] localQe = new double[nodeCnt * dof];
            IvyFEM.Lapack.DoubleMatrix localKe = new IvyFEM.Lapack.DoubleMatrix(nodeCnt * dof, nodeCnt * dof);

            //------------------------------------
            // sampling for MITC
            double[][] sampleRs = new double[3][];
            sampleRs[0] = new double[3] { 1.0 / 2.0, 0.0, 0.0 };
            sampleRs[1] = new double[3] { 0.0, 1.0 / 2.0, 0.0 };
            sampleRs[2] = new double[3] { 1.0 / 2.0, 1.0 / 2.0, 0.0 };
            double[][,][] sampleB1ss = new double[3][,][];
            IvyFEM.Lapack.DoubleMatrix[][,] sampleB2ss = new IvyFEM.Lapack.DoubleMatrix[3][,];
            IvyFEM.Lapack.DoubleMatrix[] sampleEs = new IvyFEM.Lapack.DoubleMatrix[3];
            for (int iSample = 0; iSample < 3; iSample++)
            {
                double[] sampleR = sampleRs[iSample];
                OpenTK.Vector3d[] g0s = CalcMITCNonlinearG0s(
                    xPt0s, h, sampleR);
                OpenTK.Vector3d[] curGs = CalcMITCNonlinearCurGs(
                    xPts, h, Vns, sampleR);
                OpenTK.Vector3d[][] coeffA1s;
                OpenTK.Vector3d[][,] coeffA2s;
                CalcMITCNonlinearUrgCoeffVectors(
                    dof, nodeCnt,
                    h, Vns, V1s, V2s, sampleR,
                    out coeffA1s, out coeffA2s);
                double[,][] B1s;
                IvyFEM.Lapack.DoubleMatrix[,] B2s;
                CalcMITCNonlinearStrainCoeff(
                    dof, nodeCnt,
                    curGs, coeffA1s, coeffA2s,
                    out B1s, out B2s);
                sampleB1ss[iSample] = B1s;
                sampleB2ss[iSample] = B2s;
                IvyFEM.Lapack.DoubleMatrix E = CalcMITCNonliearCurStrain(g0s, curGs);
                sampleEs[iSample] = E;
            }
            double[] sampleB113Pt1 = sampleB1ss[0][0, 2];
            double[] sampleB123Pt2 = sampleB1ss[1][1, 2];
            double[] sampleB113Pt3 = sampleB1ss[2][0, 2];
            double[] sampleB123Pt3 = sampleB1ss[2][1, 2];

            IvyFEM.Lapack.DoubleMatrix sampleB213Pt1 = sampleB2ss[0][0, 2];
            IvyFEM.Lapack.DoubleMatrix sampleB223Pt2 = sampleB2ss[1][1, 2];
            IvyFEM.Lapack.DoubleMatrix sampleB213Pt3 = sampleB2ss[2][0, 2];
            IvyFEM.Lapack.DoubleMatrix sampleB223Pt3 = sampleB2ss[2][1, 2];

            double sampleE13Pt1 = sampleEs[0][0, 2];
            double sampleE23Pt2 = sampleEs[1][1, 2];
            double sampleE13Pt3 = sampleEs[2][0, 2];
            double sampleE23Pt3 = sampleEs[2][1, 2];

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
                    double detJ = dTriFE.GetDetJacobian(L);
                    double weight = ipXY.Weights[ipXYPt];
                    double detJWeight = (1.0 / 2.0) * weight * detJ;

                    double r1 = L[1];
                    double r2 = L[2];
                    double[] r = new double[3] { r1, r2, r3 };

                    OpenTK.Vector3d[] g0s = CalcMITCNonlinearG0s(
                        xPt0s, h, r);
                    OpenTK.Vector3d[] curGs = CalcMITCNonlinearCurGs(
                        xPts, h, Vns, r);
                    OpenTK.Vector3d[][] coeffA1s;
                    OpenTK.Vector3d[][,] coeffA2s;
                    CalcMITCNonlinearUrgCoeffVectors(
                        dof, nodeCnt,
                        h, Vns, V1s, V2s, r,
                        out coeffA1s, out coeffA2s);
                    double[,][] B1s;
                    IvyFEM.Lapack.DoubleMatrix[,] B2s;
                    CalcMITCNonlinearStrainCoeff(
                        dof, nodeCnt,
                        curGs, coeffA1s, coeffA2s,
                        out B1s, out B2s);
                    IvyFEM.Lapack.DoubleMatrix E = CalcMITCNonliearCurStrain(g0s, curGs);

                    // MTIC
                    for (int k = 0; k < nodeCnt * dof; k++)
                    {
                        double coeffC = sampleB113Pt3[k] - sampleB113Pt1[k] -
                            (sampleB123Pt3[k] - sampleB123Pt2[k]);
                        // rt
                        double hB113 = sampleB113Pt1[k] + coeffC * r2;
                        double hB123 = sampleB123Pt2[k] - coeffC * r1;

                        // rt
                        B1s[0, 2][k] = hB113;
                        B1s[2, 0][k] = hB113;
                        // st
                        B1s[1, 2][k] = hB123;
                        B1s[2, 1][k] = hB123;
                    }
                    for (int p = 0; p < nodeCnt * dof; p++)
                    {
                        for (int q = 0; q < nodeCnt * dof; q++)
                        {
                            double coeffC = sampleB213Pt3[p, q] - sampleB213Pt1[p, q] -
                                (sampleB223Pt3[p, q] - sampleB223Pt2[p, q]);
                            // rt
                            double hB213 = sampleB213Pt1[p, q] + coeffC * r2;
                            double hB223 = sampleB223Pt2[p, q] - coeffC * r1;

                            // rt
                            B2s[0, 2][p, q] = hB213;
                            B2s[2, 0][p, q] = hB213;
                            // st
                            B2s[1, 2][p, q] = hB223;
                            B2s[2, 1][p, q] = hB223;

                        }
                    }
                    {
                        double coeffC = sampleE13Pt3 - sampleE13Pt1 -
                            (sampleE23Pt3 - sampleE23Pt2);
                        // rt
                        double hE13 = sampleE13Pt1 + coeffC * r2;
                        double hE23 = sampleE23Pt2 - coeffC * r1;

                        // rt
                        E[0, 2] = hE13;
                        E[2, 0] = hE13;
                        // st
                        E[1, 2] = hE23;
                        E[2, 1] = hE23;
                    }
                    //-------------------------------------

                    //-------------------------------------
                    OpenTK.Vector3d[] h0s = CalcMITCNonlinearContravariantBasisVectors(g0s);
                    double[,,,] C4 = CalcMITCNonlinearStVenantC4Tensor(
                        lambda, mu, h0s);
                    IvyFEM.Lapack.DoubleMatrix curS = CalcMITCNonlinearStress(C4, E);

                    //-------------------------------------
                    // Q
                    for (int p = 0; p < nodeCnt; p++)
                    {
                        for (int iDof = 0; iDof < dof; iDof++)
                        {
                            double tmp = 0.0;
                            for (int ig = 0; ig < 3; ig++)
                            {
                                for (int ih = 0; ih < 3; ih++)
                                {
                                    tmp += B1s[ig, ih][p * dof + iDof] * curS[ig, ih];
                                }
                            }
                            localQe[p * dof + iDof] += detJWeightZPt * detJWeight * tmp;
                        }
                    }
                    //-------------------------------------

                    //-------------------------------------
                    // K
                    for (int p = 0; p < nodeCnt; p++)
                    {
                        for (int q = 0; q < nodeCnt; q++)
                        {
                            for (int iDof = 0; iDof < dof; iDof++)
                            {
                                for (int jDof = 0; jDof < dof; jDof++)
                                {
                                    double tmp1 = 0.0;
                                    for (int ig = 0; ig < 3; ig++)
                                    {
                                        for (int ih = 0; ih < 3; ih++)
                                        {
                                            for (int e = 0; e < 3; e++)
                                            {
                                                for (int f = 0; f < 3; f++)
                                                {
                                                    tmp1 += B1s[ig, ih][p * dof + iDof] *
                                                        C4[ig, ih, e, f] * B1s[e, f][q * dof + jDof];
                                                }
                                            }
                                        }
                                    }
                                    double tmp2 = 0.0;
                                    for (int ig = 0; ig < 3; ig++)
                                    {
                                        for (int ih = 0; ih < 3; ih++)
                                        {
                                            tmp2 += curS[ig, ih] *
                                                B2s[ig, ih][p * dof + iDof, q * dof + jDof];
                                        }
                                    }
                                    double tmp = tmp1 + tmp2;
                                    localKe[p * dof + iDof, q * dof + jDof] += detJWeightZPt * detJWeight * tmp;
                                }
                            }
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
            Qe = transT * localQe;
            Ke = transT * localKe * T;
        }

        public static IvyFEM.Lapack.DoubleMatrix CalcMITCNonlinearPlateMe(
            uint feId, TriangleFE dTriFE, IList<int> dCoIds, double h, OpenTK.Vector3d[] xPt0s, double rho,
            FieldValue uFV, FieldValue vnFV, FieldValue v1FV, FieldValue v2FV)
        {
            int nodeCnt = xPt0s.Length;
            System.Diagnostics.Debug.Assert(nodeCnt == 3);

            OpenTK.Vector3d[] xPts = new OpenTK.Vector3d[nodeCnt];
            for (int iNode = 0; iNode < nodeCnt; iNode++)
            {
                int coId = dCoIds[iNode];
                double[] prevU = uFV.GetDoubleValue(coId, FieldDerivativeType.Value);
                xPts[iNode] = xPt0s[iNode] + new OpenTK.Vector3d(prevU[0], prevU[1], prevU[2]);
            }
            OpenTK.Vector3d[] Vns = new OpenTK.Vector3d[nodeCnt];
            OpenTK.Vector3d[] V1s = new OpenTK.Vector3d[nodeCnt];
            OpenTK.Vector3d[] V2s = new OpenTK.Vector3d[nodeCnt];
            System.Diagnostics.Debug.Assert(feId > 0);
            for (int iNode = 0; iNode < nodeCnt; iNode++)
            {
                int id = (int)((feId - 1) * nodeCnt + iNode);
                double[] doubleVn = vnFV.GetDoubleValue(id, FieldDerivativeType.Value);
                double[] doubleV1 = v1FV.GetDoubleValue(id, FieldDerivativeType.Value);
                double[] doubleV2 = v2FV.GetDoubleValue(id, FieldDerivativeType.Value);
                Vns[iNode] = new OpenTK.Vector3d(doubleVn[0], doubleVn[1], doubleVn[2]);
                V1s[iNode] = new OpenTK.Vector3d(doubleV1[0], doubleV1[1], doubleV1[2]);
                V2s[iNode] = new OpenTK.Vector3d(doubleV2[0], doubleV2[1], doubleV2[2]);
            }

            int uDof = 3;
            int aDof = 1;
            int bDof = 1;
            int zDof = 1;
            int dof = uDof + aDof + bDof + zDof;
            int aOffset = uDof;
            int bOffset = aOffset + aDof;
            int zOffset = bOffset + bDof;
            System.Diagnostics.Debug.Assert(dof == 6);

            // transfer matrix
            var T = CalcMITCNonlinearTransferMatrix(dof, nodeCnt, Vns, V1s, V2s);
            var transT = IvyFEM.Lapack.DoubleMatrix.Transpose(T);

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

                    OpenTK.Vector3d[] coeffB1s;
                    OpenTK.Vector3d[,] coeffB2s;
                    CalcMITCNonlinearUCoeffVectors(dof, nodeCnt, h, Vns, V1s, V2s, r, out coeffB1s, out coeffB2s);
                    // TODO:
                    // NOTE: 運動エネルギーの線形項のみ考慮している。
                    //       B2sを考慮すると方程式の線形性がなくなるのでNewton-Raphson法を使用しなければならなくなる
                    IvyFEM.Lapack.DoubleMatrix B = new IvyFEM.Lapack.DoubleMatrix(3, nodeCnt * dof);
                    for (int k = 0; k < nodeCnt * dof; k++)
                    {
                        B[0, k] = coeffB1s[k].X;
                        B[1, k] = coeffB1s[k].Y;
                        B[2, k] = coeffB1s[k].Z;
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
            var Me = transT * localMe * T;

            return Me;
        }
    }
}
