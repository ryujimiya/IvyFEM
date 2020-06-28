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
        protected static void CalcMITCThicknessStretchInitialDirectorVectors(
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

        protected static void UpdateMITCThicknessStretchDirectorVectors(
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

        protected static OpenTK.Vector3d[] CalcMITCThicknessStretchG0s(
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

            // initial strethcing
            double lambda0 = 1.0;
            double lambda1 = 0.0;

            // initial director vectors
            OpenTK.Vector3d[] Vn0s;
            OpenTK.Vector3d[] V10s;
            OpenTK.Vector3d[] V20s;
            CalcMITCThicknessStretchInitialDirectorVectors(
                xPt0s, out Vn0s, out V10s, out V20s);

            OpenTK.Vector3d ptVn0 = new OpenTK.Vector3d();
            for (int k = 0; k < nodeCnt; k++)
            {
                ptVn0 += Ns[k] * Vn0s[k];
            }
            double absVn0 = ptVn0.Length;

            // g0
            OpenTK.Vector3d[] g0s = new OpenTK.Vector3d[3];
            for (int ig = 0; ig < 3; ig++)
            {
                double deltaig3 = ig == 2 ? 1.0 : 0.0;
                g0s[ig] = new OpenTK.Vector3d();
                for (int k = 0; k < nodeCnt; k++)
                {
                    double r3Nkrg = r3 * Nrs[k][ig] + deltaig3 * Ns[k];
                    double r32Nkrg = r3 * r3 * Nrs[k][ig] + 2.0 * r3 * deltaig3 * Ns[k];
                    g0s[ig] += Nrs[k][ig] * xPt0s[k] +
                        (1.0 / 2.0) * (r3Nkrg * lambda0 + r32Nkrg * lambda1) * (h / absVn0) * Vn0s[k];
                }
            }
            return g0s;
        }

        protected static OpenTK.Vector3d[] CalcMITCThicknessStretchCurGs(
            OpenTK.Vector3d[] xPts, double h, double lambda0, double lambda1, OpenTK.Vector3d[] Vns, double[] r)
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

            OpenTK.Vector3d ptVn = new OpenTK.Vector3d();
            for (int k = 0; k < nodeCnt; k++)
            {
                ptVn += Ns[k] * Vns[k];
            }
            double absVn = ptVn.Length;

            // gt
            OpenTK.Vector3d[] curGs = new OpenTK.Vector3d[3];
            for (int ig = 0; ig < 3; ig++)
            {
                double deltaig3 = ig == 2 ? 1.0 : 0.0;
                curGs[ig] = new OpenTK.Vector3d();
                for (int k = 0; k < nodeCnt; k++)
                {
                    double r3Nkrg = r3 * Nrs[k][ig] + deltaig3 * Ns[k];
                    double r32Nkrg = r3 * r3 * Nrs[k][ig] + 2.0 * r3 * deltaig3 * Ns[k];
                    curGs[ig] += Nrs[k][ig] * xPts[k] +
                        (1.0 / 2.0) * (r3Nkrg * lambda0 + r32Nkrg * lambda1) * (h / absVn) * Vns[k];
                }
            }
            return curGs;
        }

        protected static void CalcMITCThicknessStretchUrg(
            int dof, int additionalDof, int nodeCnt,
            OpenTK.Vector3d[] us, double[] alphas, double[] betas, double deltaLambda0, double deltaLambda1,
            double h, double lambda0, double lambda1,
            OpenTK.Vector3d[] Vns, OpenTK.Vector3d[] V1s, OpenTK.Vector3d[] V2s, double[] r,
            out OpenTK.Vector3d[] u1rg, out OpenTK.Vector3d[] u2rg, out OpenTK.Vector3d[] u3rg)
        {
            System.Diagnostics.Debug.Assert(dof == 6);
            System.Diagnostics.Debug.Assert(additionalDof == 2);
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

            OpenTK.Vector3d ptVn = new OpenTK.Vector3d();
            for (int k = 0; k < nodeCnt; k++)
            {
                ptVn += Ns[k] * Vns[k];
            }
            double absVn = ptVn.Length;

            u1rg = new OpenTK.Vector3d[3];
            u2rg = new OpenTK.Vector3d[3];
            u3rg = new OpenTK.Vector3d[3];
            for (int ig = 0; ig < 3; ig++)
            {
                double deltaig3 = ig == 2 ? 1.0 : 0.0;

                for (int k = 0; k < nodeCnt; k++)
                {
                    double r3Nkrg = r3 * Nrs[k][ig] + deltaig3 * Ns[k];
                    double r32Nkrg = r3 * r3 * Nrs[k][ig] + 2.0 * r3 * deltaig3 * Ns[k];

                    u1rg[ig] += Nrs[k][ig] * us[k] +
                        (1.0 / 2.0) * (r3Nkrg * lambda0 + r32Nkrg * lambda1) *
                        (h / absVn) * (-alphas[k] * V2s[k] + betas[k] * V1s[k]) +
                        (1.0 / 2.0) * (r3Nkrg * deltaLambda0 + r32Nkrg * deltaLambda1) *
                        (h / absVn) * Vns[k];
                    u2rg[ig] +=
                        -(1.0 / 4.0) * (r3Nkrg * lambda0 + r32Nkrg * lambda1) *
                        (alphas[k] * alphas[k] + betas[k] * betas[k]) * (h / absVn) * Vns[k] +
                        (1.0 / 2.0) * (r3Nkrg * deltaLambda0 + r32Nkrg * deltaLambda1) *
                        (h / absVn) * (-alphas[k] * V2s[k] + betas[k] * V1s[k]);
                    u3rg[ig] +=
                        -(1.0 / 4.0) * (r3Nkrg * deltaLambda0 + r32Nkrg * deltaLambda1) *
                        (alphas[k] * alphas[k] + betas[k] * betas[k]) * (h / absVn) * Vns[k];
                }
            }
        }

        protected static void CalcMITCThicknessStretchUrgCoeffVectors(
            int dof, int additionalDof, int nodeCnt,
            double[] alphas, double[] betas, double deltaLambda0, double deltaLambda1, 
            double h, double lambda0, double lambda1,
            OpenTK.Vector3d[] Vns, OpenTK.Vector3d[] V1s, OpenTK.Vector3d[] V2s, double[] r,
            out OpenTK.Vector3d[][] coeffA1s, out OpenTK.Vector3d[][] coeffA2s, out OpenTK.Vector3d[][] coeffA3s)
        {
            System.Diagnostics.Debug.Assert(dof == 6);
            System.Diagnostics.Debug.Assert(additionalDof == 2);
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

            OpenTK.Vector3d ptVn = new OpenTK.Vector3d();
            for (int k = 0; k < nodeCnt; k++)
            {
                ptVn += Ns[k] * Vns[k];
            }
            double absVn = ptVn.Length;

            coeffA1s = new OpenTK.Vector3d[3][];
            for (int ig = 0; ig < 3; ig++)
            {
                coeffA1s[ig] = new OpenTK.Vector3d[nodeCnt * dof + additionalDof];
                double deltaig3 = ig == 2 ? 1.0 : 0.0;

                for (int k = 0; k < nodeCnt; k++)
                {
                    double r3Nkrg = r3 * Nrs[k][ig] + deltaig3 * Ns[k];
                    double r32Nkrg = r3 * r3 * Nrs[k][ig] + 2.0 * r3 * deltaig3 * Ns[k];
                    for (int d = 0; d < 3; d++)
                    {
                        coeffA1s[ig][k * dof + d] = Nrs[k][ig] * es[d];
                    }
                    coeffA1s[ig][k * dof + 3] = -(1.0 / 2.0) * (r3Nkrg * lambda0 + r32Nkrg * lambda1) *
                        (h / absVn) * V2s[k];
                    coeffA1s[ig][k * dof + 4] = (1.0 / 2.0) * (r3Nkrg * lambda0 + r32Nkrg * lambda1) *
                        (h / absVn) * V1s[k];
                    coeffA1s[ig][k * dof + 5] = new OpenTK.Vector3d();
                }
                {
                    OpenTK.Vector3d tmp1 = new OpenTK.Vector3d();
                    OpenTK.Vector3d tmp2 = new OpenTK.Vector3d();
                    for (int k = 0; k < nodeCnt; k++)
                    {
                        double r3Nkrg = r3 * Nrs[k][ig] + deltaig3 * Ns[k];
                        double r32Nkrg = r3 * r3 * Nrs[k][ig] + 2.0 * r3 * deltaig3 * Ns[k];
                        tmp1 += (1.0 / 2.0) * r3Nkrg * (h / absVn) * Vns[k];
                        tmp2 += (1.0 / 2.0) * r32Nkrg * (h / absVn) * Vns[k];
                    }
                    coeffA1s[ig][nodeCnt * dof] = tmp1;
                    coeffA1s[ig][nodeCnt * dof + 1] = tmp2;
                }
            }

            coeffA2s = new OpenTK.Vector3d[3][];
            for (int ig = 0; ig < 3; ig++)
            {
                coeffA2s[ig] = new OpenTK.Vector3d[nodeCnt * dof + additionalDof];
                double deltaig3 = ig == 2 ? 1.0 : 0.0;

                for (int k = 0; k < nodeCnt; k++)
                {
                    double r3Nkrg = r3 * Nrs[k][ig] + deltaig3 * Ns[k];
                    double r32Nkrg = r3 * r3 * Nrs[k][ig] + 2.0 * r3 * deltaig3 * Ns[k];
                    for (int d = 0; d < 3; d++)
                    {
                        coeffA2s[ig][k * dof + d] = new OpenTK.Vector3d();
                    }
                    coeffA2s[ig][k * dof + 3] =
                        -(1.0 / 2.0) * (r3Nkrg * lambda0 + r32Nkrg * lambda1) * alphas[k] * (h  / absVn) * Vns[k] -
                        (1.0 / 2.0) * (r3Nkrg * deltaLambda0 + r32Nkrg * deltaLambda1) * (h / absVn) * V2s[k];
                    coeffA2s[ig][k * dof + 4] =
                        -(1.0 / 2.0) * (r3Nkrg * lambda0 + r32Nkrg * lambda1) * betas[k] * (h / absVn) * Vns[k] +
                        (1.0 / 2.0) * (r3Nkrg * deltaLambda0 + r32Nkrg * deltaLambda1) * (h / absVn) * V1s[k];
                    coeffA2s[ig][k * dof + 5] = new OpenTK.Vector3d();
                }
                {
                    OpenTK.Vector3d tmp1 = new OpenTK.Vector3d();
                    OpenTK.Vector3d tmp2 = new OpenTK.Vector3d();
                    for (int k = 0; k < nodeCnt; k++)
                    {
                        double r3Nkrg = r3 * Nrs[k][ig] + deltaig3 * Ns[k];
                        double r32Nkrg = r3 * r3 * Nrs[k][ig] + 2.0 * r3 * deltaig3 * Ns[k];
                        tmp1 +=
                            (1.0 / 2.0) * r3Nkrg *
                            (h / absVn) * (-alphas[k] * V2s[k] + betas[k] * V1s[k]);
                        tmp2 +=
                            (1.0 / 2.0) * r32Nkrg *
                            (h / absVn) * (-alphas[k] * V2s[k] + betas[k] * V1s[k]);
                    }
                    coeffA2s[ig][nodeCnt * dof] = tmp1;
                    coeffA2s[ig][nodeCnt * dof + 1] = tmp2;
                }
            }

            coeffA3s = new OpenTK.Vector3d[3][];
            for (int ig = 0; ig < 3; ig++)
            {
                coeffA3s[ig] = new OpenTK.Vector3d[nodeCnt * dof + additionalDof];
                double deltaig3 = ig == 2 ? 1.0 : 0.0;

                for (int k = 0; k < nodeCnt; k++)
                {
                    double r3Nkrg = r3 * Nrs[k][ig] + deltaig3 * Ns[k];
                    double r32Nkrg = r3 * r3 * Nrs[k][ig] + 2.0 * r3 * deltaig3 * Ns[k];
                    for (int d = 0; d < 3; d++)
                    {
                        coeffA3s[ig][k * dof + d] = new OpenTK.Vector3d();
                    }
                    coeffA3s[ig][k * dof + 3] =
                        -(1.0 / 2.0) * (r3Nkrg * deltaLambda0 + r32Nkrg * deltaLambda1) *
                        alphas[k] * (h / absVn) * Vns[k];
                    coeffA3s[ig][k * dof + 4] =
                        -(1.0 / 2.0) * (r3Nkrg * deltaLambda0 + r32Nkrg * deltaLambda1) *
                        betas[k] * (h / absVn) * Vns[k];
                    coeffA3s[ig][k * dof + 5] = new OpenTK.Vector3d();
                }
                {
                    OpenTK.Vector3d tmp1 = new OpenTK.Vector3d();
                    OpenTK.Vector3d tmp2 = new OpenTK.Vector3d();
                    for (int k = 0; k < nodeCnt; k++)
                    {
                        double r3Nkrg = r3 * Nrs[k][ig] + deltaig3 * Ns[k];
                        double r32Nkrg = r3 * r3 * Nrs[k][ig] + 2.0 * r3 * deltaig3 * Ns[k];
                        tmp1 +=
                            -(1.0 / 4.0) * r3Nkrg * (alphas[k] * alphas[k] + betas[k] * betas[k]) * (h / absVn) * Vns[k];
                        tmp2 +=
                            -(1.0 / 4.0) * r32Nkrg * (alphas[k] * alphas[k] + betas[k] * betas[k]) * (h / absVn) * Vns[k];
                    }
                    coeffA3s[ig][nodeCnt * dof] = tmp1;
                    coeffA3s[ig][nodeCnt * dof + 1] = tmp2;
                }
            }
        }

        protected static void CalcMITCThicknessStretchStrainCoeff(
            int dof, int additionalDof, int nodeCnt,
            OpenTK.Vector3d[] u1rg, OpenTK.Vector3d[] u2rg, OpenTK.Vector3d[] u3rg,
            OpenTK.Vector3d[] curGs,
            OpenTK.Vector3d[][] coeffA1s, OpenTK.Vector3d[][] coeffA2s, OpenTK.Vector3d[][] coeffA3s,
            out double[,][] B1s, out double[,][] B2s)
        {
            System.Diagnostics.Debug.Assert(dof == 6);
            System.Diagnostics.Debug.Assert(additionalDof == 2);
            System.Diagnostics.Debug.Assert(nodeCnt == 3);

            B1s = new double[3, 3][];
            for (int ig = 0; ig < 3; ig++)
            {
                for (int ih = 0; ih < 3; ih++)
                {
                    B1s[ig, ih] = new double[nodeCnt* dof + additionalDof];
                    for (int k = 0; k < (nodeCnt * dof + additionalDof); k++)
                    {
                        B1s[ig, ih][k] = (1.0 / 2.0) * (
                            OpenTK.Vector3d.Dot(curGs[ih], coeffA1s[ig][k]) +
                            OpenTK.Vector3d.Dot(curGs[ig], coeffA1s[ih][k]));
                    }
                }
            }

            B2s = new double[3, 3][];
            for (int ig = 0; ig < 3; ig++)
            {
                for (int ih = 0; ih < 3; ih++)
                {
                    B2s[ig, ih] = new double[nodeCnt * dof + additionalDof];
                    for (int k = 0; k < (nodeCnt * dof + additionalDof); k++)
                    {
                        //
                        OpenTK.Vector3d a23g = coeffA2s[ig][k] + coeffA3s[ig][k];
                        OpenTK.Vector3d a123g = coeffA1s[ig][k] + coeffA2s[ig][k] + coeffA3s[ig][k];
                        OpenTK.Vector3d u123rg = u1rg[ig] + u2rg[ig] + u3rg[ig];
                        //
                        OpenTK.Vector3d a23h = coeffA2s[ih][k] + coeffA3s[ih][k];
                        OpenTK.Vector3d a123h = coeffA1s[ih][k] + coeffA2s[ih][k] + coeffA3s[ih][k];
                        OpenTK.Vector3d u123rh = u1rg[ih] + u2rg[ih] + u3rg[ih];

                        B2s[ig, ih][k] += (1.0 / 2.0) * (
                            OpenTK.Vector3d.Dot(a23g, curGs[ih]) +
                            OpenTK.Vector3d.Dot(a23h, curGs[ig]) +
                            OpenTK.Vector3d.Dot(a123g, u123rh) +
                            OpenTK.Vector3d.Dot(a123h, u123rg));
                    }
                }
            }
        }

        protected static void CalcMITCThicknessStretchUrgCoeffVectorsU(
            int dof, int additionalDof, int nodeCnt,
            double[] alphas, double[] betas, double deltaLambda0, double deltaLambda1,
            double h, double lambda0, double lambda1,
            OpenTK.Vector3d[] Vns, OpenTK.Vector3d[] V1s, OpenTK.Vector3d[] V2s, double[] r,
            out OpenTK.Vector3d[][,] coeffAu2s, out OpenTK.Vector3d[][,] coeffAu3s)
        {
            System.Diagnostics.Debug.Assert(dof == 6);
            System.Diagnostics.Debug.Assert(additionalDof == 2);
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

            OpenTK.Vector3d ptVn = new OpenTK.Vector3d();
            for (int k = 0; k < nodeCnt; k++)
            {
                ptVn += Ns[k] * Vns[k];
            }
            double absVn = ptVn.Length;

            //----------------------------------------------------------------------
            coeffAu2s = new OpenTK.Vector3d[3][,];
            for (int ig = 0; ig < 3; ig++)
            {
                coeffAu2s[ig] = new OpenTK.Vector3d[nodeCnt * dof + additionalDof, nodeCnt * dof + additionalDof];
                double deltaig3 = ig == 2 ? 1.0 : 0.0;

                // u

                // alpha
                for (int iNode = 0; iNode < nodeCnt; iNode++)
                {
                    // u

                    // alpha
                    double r3Nirg = r3 * Nrs[iNode][ig] + deltaig3 * Ns[iNode];
                    double r32Nirg = r3 * r3 * Nrs[iNode][ig] + 2.0 * r3 * deltaig3 * Ns[iNode];
                    coeffAu2s[ig][iNode * dof + 3, iNode * dof + 3] =
                        -(1.0 / 2.0) * (r3Nirg * lambda0 + r32Nirg * lambda1) * (h / absVn) * Vns[iNode];

                    // beta

                    // lambda0
                    {
                        OpenTK.Vector3d tmp = new OpenTK.Vector3d();
                        for (int kNode = 0; kNode < nodeCnt; kNode++)
                        {
                            double r3Nkrg = r3 * Nrs[kNode][ig] + deltaig3 * Ns[kNode];
                            double r32Nkrg = r3 * r3 * Nrs[kNode][ig] + 2.0 * r3 * deltaig3 * Ns[kNode];
                            tmp += -(1.0 / 2.0) * r3Nkrg * (h / absVn) * V2s[kNode];
                        }
                        coeffAu2s[ig][iNode * dof + 3, nodeCnt * dof] = tmp;
                    }

                    // lambda1
                    {
                        OpenTK.Vector3d tmp = new OpenTK.Vector3d();
                        for (int kNode = 0; kNode < nodeCnt; kNode++)
                        {
                            double r3Nkrg = r3 * Nrs[kNode][ig] + deltaig3 * Ns[kNode];
                            double r32Nkrg = r3 * r3 * Nrs[kNode][ig] + 2.0 * r3 * deltaig3 * Ns[kNode];
                            tmp += -(1.0 / 2.0) * r32Nkrg * (h / absVn) * V2s[kNode];
                        }
                        coeffAu2s[ig][iNode * dof + 3, nodeCnt * dof + 1] = tmp;
                    }
                }

                // betak
                for (int iNode = 0; iNode < nodeCnt; iNode++)
                {
                    // u

                    // alpha

                    // beta
                    double r3Nirg = r3 * Nrs[iNode][ig] + deltaig3 * Ns[iNode];
                    double r32Nirg = r3 * r3 * Nrs[iNode][ig] + 2.0 * r3 * deltaig3 * Ns[iNode];
                    coeffAu2s[ig][iNode * dof + 4, iNode * dof + 4] =
                        -(1.0 / 2.0) * (r3Nirg * lambda0 + r32Nirg * lambda1) * (h / absVn) * Vns[iNode];

                    // lambda0
                    {
                        OpenTK.Vector3d tmp = new OpenTK.Vector3d();
                        for (int kNode = 0; kNode < nodeCnt; kNode++)
                        {
                            double r3Nkrg = r3 * Nrs[kNode][ig] + deltaig3 * Ns[kNode];
                            double r32Nkrg = r3 * r3 * Nrs[kNode][ig] + 2.0 * r3 * deltaig3 * Ns[kNode];
                            tmp += (1.0 / 2.0) * r3Nkrg * (h / absVn) * V1s[kNode];
                        }
                        coeffAu2s[ig][iNode * dof + 4, nodeCnt * dof] = tmp;
                    }

                    // lambda1
                    {
                        OpenTK.Vector3d tmp = new OpenTK.Vector3d();
                        for (int kNode = 0; kNode < nodeCnt; kNode++)
                        {
                            double r3Nkrg = r3 * Nrs[kNode][ig] + deltaig3 * Ns[kNode];
                            double r32Nkrg = r3 * r3 * Nrs[kNode][ig] + 2.0 * r3 * deltaig3 * Ns[kNode];
                            tmp += (1.0 / 2.0) * r32Nkrg * (h / absVn) * V1s[kNode];
                        }
                        coeffAu2s[ig][iNode * dof + 4, nodeCnt * dof + 1] = tmp;
                    }
                }

                // lambdap0
                {
                    // u

                    // alpha
                    for (int jNode = 0; jNode < nodeCnt; jNode++)
                    {
                        OpenTK.Vector3d tmp = new OpenTK.Vector3d();
                        for (int kNode = 0; kNode < nodeCnt; kNode++)
                        {
                            double r3Nkrg = r3 * Nrs[kNode][ig] + deltaig3 * Ns[kNode];
                            double r32Nkrg = r3 * r3 * Nrs[kNode][ig] + 2.0 * r3 * deltaig3 * Ns[kNode];
                            tmp += -(1.0 / 2.0) * r3Nkrg * (h / absVn) * V2s[kNode];
                        }
                        coeffAu2s[ig][nodeCnt * dof, jNode * dof + 3] = tmp;
                    }

                    // beta
                    for (int jNode = 0; jNode < nodeCnt; jNode++)
                    {
                        OpenTK.Vector3d tmp = new OpenTK.Vector3d();
                        for (int kNode = 0; kNode < nodeCnt; kNode++)
                        {
                            double r3Nkrg = r3 * Nrs[kNode][ig] + deltaig3 * Ns[kNode];
                            double r32Nkrg = r3 * r3 * Nrs[kNode][ig] + 2.0 * r3 * deltaig3 * Ns[kNode];
                            tmp += (1.0 / 2.0) * r3Nkrg * (h / absVn) * V1s[kNode];
                        }
                        coeffAu2s[ig][nodeCnt * dof, jNode * dof + 4] = tmp;
                    }

                    // lambda0

                    // lambda1
                }

                // lambda1
                {
                    // u

                    // alpha
                    for (int jNode = 0; jNode < nodeCnt; jNode++)
                    {
                        OpenTK.Vector3d tmp = new OpenTK.Vector3d();
                        for (int kNode = 0; kNode < nodeCnt; kNode++)
                        {
                            double r3Nkrg = r3 * Nrs[kNode][ig] + deltaig3 * Ns[kNode];
                            double r32Nkrg = r3 * r3 * Nrs[kNode][ig] + 2.0 * r3 * deltaig3 * Ns[kNode];
                            tmp += -(1.0 / 2.0) * r32Nkrg * (h / absVn) * V2s[kNode];
                        }
                        coeffAu2s[ig][nodeCnt * dof + 1, jNode * dof + 3] = tmp;
                    }

                    // beta
                    for (int jNode = 0; jNode < nodeCnt; jNode++)
                    {
                        OpenTK.Vector3d tmp = new OpenTK.Vector3d();
                        for (int kNode = 0; kNode < nodeCnt; kNode++)
                        {
                            double r3Nkrg = r3 * Nrs[kNode][ig] + deltaig3 * Ns[kNode];
                            double r32Nkrg = r3 * r3 * Nrs[kNode][ig] + 2.0 * r3 * deltaig3 * Ns[kNode];
                            tmp += (1.0 / 2.0) * r32Nkrg * (h / absVn) * V1s[kNode];
                        }
                        coeffAu2s[ig][nodeCnt * dof + 1, jNode * dof + 4] = tmp;
                    }

                    // lambda0

                    // lambda1
                }
            }

            //----------------------------------------------------------------------
            coeffAu3s = new OpenTK.Vector3d[3][,];
            for (int ig = 0; ig < 3; ig++)
            {
                coeffAu3s[ig] = new OpenTK.Vector3d[nodeCnt * dof + additionalDof, nodeCnt * dof + additionalDof];
                double deltaig3 = ig == 2 ? 1.0 : 0.0;

                // u

                // alpha
                for (int iNode = 0; iNode < nodeCnt; iNode++)
                {
                    double r3Nirg = r3 * Nrs[iNode][ig] + deltaig3 * Ns[iNode];
                    double r32Nirg = r3 * r3 * Nrs[iNode][ig] + 2.0 * r3 * deltaig3 * Ns[iNode];
                    
                    // u

                    // alpha
                    coeffAu3s[ig][iNode * dof + 3, iNode * dof + 3] =
                        -(1.0 / 2.0) * (r3Nirg * deltaLambda0 + r32Nirg * deltaLambda1) * (h / absVn) * Vns[iNode];

                    // beta

                    // lambda0
                    {
                        OpenTK.Vector3d tmp = new OpenTK.Vector3d();
                        for (int kNode = 0; kNode < nodeCnt; kNode++)
                        {
                            double r3Nkrg = r3 * Nrs[kNode][ig] + deltaig3 * Ns[kNode];
                            double r32Nkrg = r3 * r3 * Nrs[kNode][ig] + 2.0 * r3 * deltaig3 * Ns[kNode];
                            tmp += -(1.0 / 2.0) * r3Nkrg * alphas[kNode] * (h / absVn) * Vns[kNode];
                        }
                        coeffAu3s[ig][iNode * dof + 3, nodeCnt * dof] = tmp;
                    }

                    // lambda1
                    {
                        OpenTK.Vector3d tmp = new OpenTK.Vector3d();
                        for (int kNode = 0; kNode < nodeCnt; kNode++)
                        {
                            double r3Nkrg = r3 * Nrs[kNode][ig] + deltaig3 * Ns[kNode];
                            double r32Nkrg = r3 * r3 * Nrs[kNode][ig] + 2.0 * r3 * deltaig3 * Ns[kNode];
                            tmp += -(1.0 / 2.0) * r32Nkrg * alphas[kNode] * (h / absVn) * Vns[kNode];
                        }
                        coeffAu3s[ig][iNode * dof + 3, nodeCnt * dof + 1] = tmp;
                    }
                }

                // beta
                for (int iNode = 0; iNode < nodeCnt; iNode++)
                {
                    // u

                    // alpha

                    // beta
                    double r3Nirg = r3 * Nrs[iNode][ig] + deltaig3 * Ns[iNode];
                    double r32Nirg = r3 * r3 * Nrs[iNode][ig] + 2.0 * r3 * deltaig3 * Ns[iNode];
                    coeffAu3s[ig][iNode * dof + 4, iNode * dof + 4] =
                        -(1.0 / 2.0) * (r3Nirg * deltaLambda0 + r32Nirg * deltaLambda1) * (h / absVn) * Vns[iNode];

                    // lambda0
                    {
                        OpenTK.Vector3d tmp = new OpenTK.Vector3d();
                        for (int kNode = 0; kNode < nodeCnt; kNode++)
                        {
                            double r3Nkrg = r3 * Nrs[kNode][ig] + deltaig3 * Ns[kNode];
                            double r32Nkrg = r3 * r3 * Nrs[kNode][ig] + 2.0 * r3 * deltaig3 * Ns[kNode];
                            tmp += -(1.0 / 2.0) * r3Nkrg * betas[kNode] * (h / absVn) * Vns[kNode];
                        }
                        coeffAu3s[ig][iNode * dof + 4, nodeCnt * dof] = tmp;
                    }

                    // lambda1
                    {
                        OpenTK.Vector3d tmp = new OpenTK.Vector3d();
                        for (int kNode = 0; kNode < nodeCnt; kNode++)
                        {
                            double r3Nkrg = r3 * Nrs[kNode][ig] + deltaig3 * Ns[kNode];
                            double r32Nkrg = r3 * r3 * Nrs[kNode][ig] + 2.0 * r3 * deltaig3 * Ns[kNode];
                            tmp += -(1.0 / 2.0) * r32Nkrg * betas[kNode] * (h / absVn) * Vns[kNode];
                        }
                        coeffAu3s[ig][iNode * dof + 4, nodeCnt * dof + 1] = tmp;
                    }
                }

                // lambda0
                {
                    // u

                    // alpha
                    for (int jNode = 0; jNode < nodeCnt; jNode++)
                    {
                        OpenTK.Vector3d tmp = new OpenTK.Vector3d();
                        for (int kNode = 0; kNode < nodeCnt; kNode++)
                        {
                            double r3Nkrg = r3 * Nrs[kNode][ig] + deltaig3 * Ns[kNode];
                            double r32Nkrg = r3 * r3 * Nrs[kNode][ig] + 2.0 * r3 * deltaig3 * Ns[kNode];
                            tmp += -(1.0 / 2.0) * r3Nkrg * alphas[kNode] * (h / absVn) * Vns[kNode];
                        }
                        coeffAu3s[ig][nodeCnt * dof, jNode * dof + 3] = tmp;
                    }

                    // beta
                    for (int jNode = 0; jNode < nodeCnt; jNode++)
                    {
                        OpenTK.Vector3d tmp = new OpenTK.Vector3d();
                        for (int kNode = 0; kNode < nodeCnt; kNode++)
                        {
                            double r3Nkrg = r3 * Nrs[kNode][ig] + deltaig3 * Ns[kNode];
                            double r32Nkrg = r3 * r3 * Nrs[kNode][ig] + 2.0 * r3 * deltaig3 * Ns[kNode];
                            tmp += -(1.0 / 2.0) * r3Nkrg * betas[kNode] * (h / absVn) * Vns[kNode];
                        }
                        coeffAu3s[ig][nodeCnt * dof, jNode * dof + 4] = tmp;
                    }

                    // lambda0

                    // lambda1
                }

                // lambda1
                {
                    // u

                    // alpha
                    for (int jNode = 0; jNode < nodeCnt; jNode++)
                    {
                        OpenTK.Vector3d tmp = new OpenTK.Vector3d();
                        for (int kNode = 0; kNode < nodeCnt; kNode++)
                        {
                            double r3Nkrg = r3 * Nrs[kNode][ig] + deltaig3 * Ns[kNode];
                            double r32Nkrg = r3 * r3 * Nrs[kNode][ig] + 2.0 * r3 * deltaig3 * Ns[kNode];
                            tmp += -(1.0 / 2.0) * r32Nkrg * alphas[kNode] * (h / absVn) * Vns[kNode];
                        }
                        coeffAu3s[ig][nodeCnt * dof + 1, jNode * dof + 3] = tmp;
                    }

                    // beta
                    for (int jNode = 0; jNode < nodeCnt; jNode++)
                    {
                        OpenTK.Vector3d tmp = new OpenTK.Vector3d();
                        for (int kNode = 0; kNode < nodeCnt; kNode++)
                        {
                            double r3Nkrg = r3 * Nrs[kNode][ig] + deltaig3 * Ns[kNode];
                            double r32Nkrg = r3 * r3 * Nrs[kNode][ig] + 2.0 * r3 * deltaig3 * Ns[kNode];
                            tmp += -(1.0 / 2.0) * r32Nkrg * betas[kNode] * (h / absVn) * Vns[kNode];
                        }
                        coeffAu3s[ig][nodeCnt * dof + 1, jNode * dof + 4] = tmp;
                    }

                    // lambda0

                    // lambda1
                }
            }
        }

        protected static void CalcMITCThicknessStretchStrainCoeffU(
            int dof, int additionalDof, int nodeCnt,
            OpenTK.Vector3d[] u1rg, OpenTK.Vector3d[] u2rg, OpenTK.Vector3d[] u3rg,
            OpenTK.Vector3d[] curGs,
            OpenTK.Vector3d[][] coeffA1s, OpenTK.Vector3d[][] coeffA2s, OpenTK.Vector3d[][] coeffA3s,
            OpenTK.Vector3d[][,] coeffAu2s, OpenTK.Vector3d[][,] coeffAu3s,
            out double[,][,] Bu2s)
        {
            System.Diagnostics.Debug.Assert(dof == 6);
            System.Diagnostics.Debug.Assert(additionalDof == 2);
            System.Diagnostics.Debug.Assert(nodeCnt == 3);

            Bu2s = new double[3, 3][,];
            for (int ig = 0; ig < 3; ig++)
            {
                for (int ih = 0; ih < 3; ih++)
                {
                    Bu2s[ig, ih] = new double[nodeCnt * dof + additionalDof, nodeCnt * dof + additionalDof];
                    for (int i = 0; i < (nodeCnt * dof + additionalDof); i++)
                    {
                        for (int j = 0; j < (nodeCnt * dof + additionalDof); j++)
                        {
                            //
                            OpenTK.Vector3d au23gij = coeffAu2s[ig][i, j] + coeffAu3s[ig][i, j];
                            //
                            OpenTK.Vector3d au23hij = coeffAu2s[ih][i, j] + coeffAu3s[ih][i, j];

                            //
                            OpenTK.Vector3d a123gi = coeffA1s[ig][i] + coeffA2s[ig][i] + coeffA3s[ig][i];
                            OpenTK.Vector3d a123gj = coeffA1s[ig][j] + coeffA2s[ig][j] + coeffA3s[ig][j];
                            OpenTK.Vector3d u123rg = u1rg[ig] + u2rg[ig] + u3rg[ig];
                            //
                            OpenTK.Vector3d a123hi = coeffA1s[ih][i] + coeffA2s[ih][i] + coeffA3s[ih][i];
                            OpenTK.Vector3d a123hj = coeffA1s[ih][j] + coeffA2s[ih][j] + coeffA3s[ih][j];
                            OpenTK.Vector3d u123rh = u1rg[ih] + u2rg[ih] + u3rg[ih];

                            Bu2s[ig, ih][i, j] = (1.0 / 2.0) * (
                                OpenTK.Vector3d.Dot(au23gij, curGs[ih]) +
                                OpenTK.Vector3d.Dot(au23hij, curGs[ig]) +
                                OpenTK.Vector3d.Dot(au23gij, u123rh) +
                                OpenTK.Vector3d.Dot(a123gi, a123hj) +
                                OpenTK.Vector3d.Dot(au23hij, u123rg) +
                                OpenTK.Vector3d.Dot(a123hi, a123gj));
                        }
                    }
                }
            }
        }

        protected static OpenTK.Vector3d[] CalcMITCThicknessStretchContravariantBasisVectors(OpenTK.Vector3d[] gs)
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
        protected static void CalcMITCThicknessStretchStVenantC4StressAndTensor(
            double lambda, double mu,
            OpenTK.Vector3d[] g0s, OpenTK.Vector3d[] h0s,
            IvyFEM.Lapack.DoubleMatrix cgh,
            out IvyFEM.Lapack.DoubleMatrix S, out double[,,,] C4)
        {
            IvyFEM.Lapack.DoubleMatrix E = new IvyFEM.Lapack.DoubleMatrix(3, 3);
            for (int ig = 0; ig < 3; ig++)
            {
                for (int ih = 0; ih < 3; ih++)
                {
                    E[ig, ih] = (1.0 / 2.0) * (
                        cgh[ig, ih] -
                        OpenTK.Vector3d.Dot(g0s[ig], g0s[ih]));
                }
            }

            C4 = new double[3, 3, 3, 3];
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

            S = new IvyFEM.Lapack.DoubleMatrix(3, 3);
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
        }

        public static IvyFEM.Lapack.DoubleMatrix CalcMITCThicknessStretchTransferMatrix(
            int dof, int additionalDof, int nodeCnt,
            OpenTK.Vector3d[] Vns, OpenTK.Vector3d[] V1s, OpenTK.Vector3d[] V2s)
        {
            IvyFEM.Lapack.DoubleMatrix T = new Lapack.DoubleMatrix(
                nodeCnt * dof + additionalDof, nodeCnt * dof + additionalDof);
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
            T[nodeCnt * dof, nodeCnt * dof] = 1.0;
            T[nodeCnt * dof + 1, nodeCnt * dof + 1] = 1.0;
            return T;
        }

        public static void InitMITCThicknessStretchPlateNodeValues(int coId, FieldValue uFV)
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

        public static void InitMITCThicknessStretchPlateElementNodeValues(
            uint feId, OpenTK.Vector3d[] xPt0s,
            FieldValue vnFV, FieldValue v1FV, FieldValue v2FV, FieldValue lambdaFV)
        {
            int nodeCnt = xPt0s.Length;

            //----------------------------
            // director vectors
            OpenTK.Vector3d[] V0ns;
            OpenTK.Vector3d[] V01s;
            OpenTK.Vector3d[] V02s;
            CalcMITCThicknessStretchInitialDirectorVectors(xPt0s, out V0ns, out V01s, out V02s);

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

            //----------------------------
            // λ0 λ1
            int lDof = 2;
            System.Diagnostics.Debug.Assert(lambdaFV.Dof == lDof);
            System.Diagnostics.Debug.Assert(lambdaFV.NodeType == FieldValueNodeType.Bubble);
            double[] lambdaValues = lambdaFV.GetDoubleValues(FieldDerivativeType.Value);
            System.Diagnostics.Debug.Assert(feId > 0);
            {
                int id = (int)((feId - 1) * 2);
                lambdaValues[id] = 1.0;
                lambdaValues[id + 1] = 0.0;
            }
        }

        public static void UpdateMITCThicknessStretchPlateNodeValues(int coId, double[] curNodeUg, FieldValue uFV)
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

        public static void UpdateMITCThicknessStretchPlateElementNodeValues(
            uint feId, int nodeCnt, double[] curUg,
            FieldValue vnFV, FieldValue v1FV, FieldValue v2FV, FieldValue lambdaFV)
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
            //----------------------------

            //----------------------------
            // λ0 λ1
            int lDof = 2;
            System.Diagnostics.Debug.Assert(lambdaFV.Dof == lDof);
            System.Diagnostics.Debug.Assert(lambdaFV.NodeType == FieldValueNodeType.Bubble);
            double[] lambdaValues = lambdaFV.GetDoubleValues(FieldDerivativeType.Value);
            System.Diagnostics.Debug.Assert(feId > 0);
            double prevLambda0;
            double prevLambda1;
            {
                int id = (int)((feId - 1) * 2);
                prevLambda0 = lambdaValues[id];
                prevLambda1 = lambdaValues[id + 1];
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
            int additionalDof = 2;

            // transfer matrix
            var T = CalcMITCThicknessStretchTransferMatrix(dof, additionalDof, nodeCnt, prevVns, prevV1s, prevV2s);
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
            double lambda0 = curUl[nodeCnt * dof] + prevLambda0;
            double lambda1 = curUl[nodeCnt * dof + 1] + prevLambda1;
            System.Diagnostics.Debug.Assert(lambda0 > Constants.PrecisionLowerLimit); // lamda0 > 0
            //System.Diagnostics.Debug.Assert(lambda1 >= 0); // lamda1 >= 0 とは限らない？

            //----------------------------
            // director vector
            OpenTK.Vector3d[] Vns;
            OpenTK.Vector3d[] V1s;
            OpenTK.Vector3d[] V2s;
            UpdateMITCThicknessStretchDirectorVectors(
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

            //----------------------------
            // stretch
            {
                int id = (int)((feId - 1) * 2);
                lambdaValues[id] = lambda0;
                lambdaValues[id + 1] = lambda1;
            }
        }

        public static void CalcMITCStVenantThicknessStretchPateKe(
            double lambda, double mu,
            uint feId, TriangleFE dTriFE, IList<int> dCoIds, double h, OpenTK.Vector3d[] xPt0s,
            FieldValue uFV, FieldValue vnFV, FieldValue v1FV, FieldValue v2FV, FieldValue lambdaFV,
            double[] curUg,
            out double[] Qe0, out double[] Qe2, out IvyFEM.Lapack.DoubleMatrix Ke1, out IvyFEM.Lapack.DoubleMatrix Ke2)
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
            double lambda0;
            double lambda1;
            {
                int id = (int)(feId - 1); // Note:ここは2倍しない(番号指定)
                double[] doubleLambda = lambdaFV.GetDoubleValue(id, FieldDerivativeType.Value);
                lambda0 = doubleLambda[0];
                lambda1 = doubleLambda[1];
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
            int additionalDof = 2;

            // transfer matrix
            var T = CalcMITCThicknessStretchTransferMatrix(dof, additionalDof, nodeCnt, Vns, V1s, V2s);
            var transT = IvyFEM.Lapack.DoubleMatrix.Transpose(T);

            //-------------------------------------------------
            // 現在の節点値
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
            double deltaLambda0 = curUl[nodeCnt * dof];
            double deltaLambda1 = curUl[nodeCnt * dof + 1];

            double[] localQe0 = new double[nodeCnt * dof + additionalDof];
            double[] localQe2 = new double[nodeCnt * dof + additionalDof];
            IvyFEM.Lapack.DoubleMatrix localKe1 = new IvyFEM.Lapack.DoubleMatrix(
                nodeCnt * dof + additionalDof, nodeCnt * dof + additionalDof);
            IvyFEM.Lapack.DoubleMatrix localKe2 = new IvyFEM.Lapack.DoubleMatrix(
                nodeCnt * dof + additionalDof, nodeCnt * dof + additionalDof);

            //------------------------------------
            // sampling for MITC
            double[][] sampleRs = new double[3][];
            sampleRs[0] = new double[3] { 1.0 / 2.0, 0.0, 0.0 };
            sampleRs[1] = new double[3] { 0.0, 1.0 / 2.0, 0.0 };
            sampleRs[2] = new double[3] { 1.0 / 2.0, 1.0 / 2.0, 0.0 };
            double[][,][] sampleB1ss = new double[3][,][];
            double[][,][] sampleB2ss = new double[3][,][];
            double[][,][,] sampleBu2ss = new double[3][,][,];
            IvyFEM.Lapack.DoubleMatrix[] sampleCghs = new IvyFEM.Lapack.DoubleMatrix[3];
            for (int iSample = 0; iSample < 3; iSample++)
            {
                double[] sampleR = sampleRs[iSample];
                OpenTK.Vector3d[] g0s = CalcMITCThicknessStretchG0s(
                    xPt0s, h, sampleR);
                OpenTK.Vector3d[] curGs = CalcMITCThicknessStretchCurGs(
                    xPts, h, lambda0, lambda1, Vns, sampleR);
                OpenTK.Vector3d[] u1rg;
                OpenTK.Vector3d[] u2rg;
                OpenTK.Vector3d[] u3rg;
                CalcMITCThicknessStretchUrg(
                    dof, additionalDof, nodeCnt,
                    us, alphas, betas, deltaLambda0, deltaLambda1,
                    h, lambda0, lambda1,
                    Vns, V1s, V2s, sampleR,
                    out u1rg, out u2rg, out u3rg);
                OpenTK.Vector3d[][] coeffA1s;
                OpenTK.Vector3d[][] coeffA2s;
                OpenTK.Vector3d[][] coeffA3s;
                CalcMITCThicknessStretchUrgCoeffVectors(
                    dof, additionalDof, nodeCnt,
                    alphas, betas, deltaLambda0, deltaLambda1,
                    h, lambda0, lambda1, Vns, V1s, V2s, sampleR,
                    out coeffA1s, out coeffA2s, out coeffA3s);
                double[,][] B1s;
                double[,][] B2s;
                CalcMITCThicknessStretchStrainCoeff(
                    dof, additionalDof, nodeCnt,
                    u1rg, u2rg, u3rg,
                    curGs, coeffA1s, coeffA2s, coeffA3s,
                    out B1s, out B2s);
                sampleB1ss[iSample] = B1s;
                sampleB2ss[iSample] = B2s;
                OpenTK.Vector3d[][,] coeffAu2s;
                OpenTK.Vector3d[][,] coeffAu3s;
                CalcMITCThicknessStretchUrgCoeffVectorsU(
                    dof, additionalDof, nodeCnt,
                    alphas, betas, deltaLambda0, deltaLambda1,
                    h, lambda0, lambda1,
                    Vns, V1s, V2s, sampleR,
                    out coeffAu2s, out coeffAu3s);
                double[,][,] Bu2s;
                CalcMITCThicknessStretchStrainCoeffU(
                    dof, additionalDof, nodeCnt,
                    u1rg, u2rg, u3rg,
                    curGs,
                    coeffA1s, coeffA2s, coeffA3s,
                    coeffAu2s, coeffAu3s,
                    out Bu2s);
                sampleBu2ss[iSample] = Bu2s;
                IvyFEM.Lapack.DoubleMatrix cgh = CalcMITCMooneyRivlinC(curGs);
                sampleCghs[iSample] = cgh;
            }
            double[] sampleB113Pt1 = sampleB1ss[0][0, 2];
            double[] sampleB123Pt2 = sampleB1ss[1][1, 2];
            double[] sampleB113Pt3 = sampleB1ss[2][0, 2];
            double[] sampleB123Pt3 = sampleB1ss[2][1, 2];

            double[] sampleB213Pt1 = sampleB2ss[0][0, 2];
            double[] sampleB223Pt2 = sampleB2ss[1][1, 2];
            double[] sampleB213Pt3 = sampleB2ss[2][0, 2];
            double[] sampleB223Pt3 = sampleB2ss[2][1, 2];

            double[,] sampleBu213Pt1 = sampleBu2ss[0][0, 2];
            double[,] sampleBu223Pt2 = sampleBu2ss[1][1, 2];
            double[,] sampleBu213Pt3 = sampleBu2ss[2][0, 2];
            double[,] sampleBu223Pt3 = sampleBu2ss[2][1, 2];

            double sampleCgh13Pt1 = sampleCghs[0][0, 2];
            double sampleCgh23Pt2 = sampleCghs[1][1, 2];
            double sampleCgh13Pt3 = sampleCghs[2][0, 2];
            double sampleCgh23Pt3 = sampleCghs[2][1, 2];
            //------------------------------------

            IntegrationPoints ipZ = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point3);
            IntegrationPoints ipXY = TriangleFE.GetIntegrationPoints(TriangleIntegrationPointCount.Point7);
            System.Diagnostics.Debug.Assert(ipZ.Ls.Length == 3);
            System.Diagnostics.Debug.Assert(ipXY.Ls.Length == 7);
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

                    OpenTK.Vector3d[] g0s = CalcMITCThicknessStretchG0s(
                        xPt0s, h, r);
                    OpenTK.Vector3d[] curGs = CalcMITCThicknessStretchCurGs(
                        xPts, h, lambda0, lambda1, Vns, r);
                    OpenTK.Vector3d[] u1rg;
                    OpenTK.Vector3d[] u2rg;
                    OpenTK.Vector3d[] u3rg;
                    CalcMITCThicknessStretchUrg(
                        dof, additionalDof, nodeCnt,
                        us, alphas, betas, deltaLambda0, deltaLambda1,
                        h, lambda0, lambda1,
                        Vns, V1s, V2s, r,
                        out u1rg, out u2rg, out u3rg);
                    OpenTK.Vector3d[][] coeffA1s;
                    OpenTK.Vector3d[][] coeffA2s;
                    OpenTK.Vector3d[][] coeffA3s;
                    CalcMITCThicknessStretchUrgCoeffVectors(
                        dof, additionalDof, nodeCnt,
                        alphas, betas, deltaLambda0, deltaLambda1,
                        h, lambda0, lambda1, Vns, V1s, V2s, r,
                        out coeffA1s, out coeffA2s, out coeffA3s);
                    double[,][] B1s;
                    double[,][] B2s;
                    CalcMITCThicknessStretchStrainCoeff(
                        dof, additionalDof, nodeCnt,
                        u1rg, u2rg, u3rg,
                        curGs, coeffA1s, coeffA2s, coeffA3s,
                        out B1s, out B2s);
                    OpenTK.Vector3d[][,] coeffAu2s;
                    OpenTK.Vector3d[][,] coeffAu3s;
                    CalcMITCThicknessStretchUrgCoeffVectorsU(
                        dof, additionalDof, nodeCnt,
                        alphas, betas, deltaLambda0, deltaLambda1,
                        h, lambda0, lambda1,
                        Vns, V1s, V2s, r,
                        out coeffAu2s, out coeffAu3s);
                    double[,][,] Bu2s;
                    CalcMITCThicknessStretchStrainCoeffU(
                        dof, additionalDof, nodeCnt,
                        u1rg, u2rg, u3rg,
                        curGs,
                        coeffA1s, coeffA2s, coeffA3s,
                        coeffAu2s, coeffAu3s,
                        out Bu2s);
                    IvyFEM.Lapack.DoubleMatrix cgh = CalcMITCMooneyRivlinC(curGs);

                    //-------------------------------------
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

                    for (int k = 0; k < nodeCnt * dof; k++)
                    {
                        double coeffC = sampleB213Pt3[k] - sampleB213Pt1[k] -
                            (sampleB223Pt3[k] - sampleB223Pt2[k]);
                        // rt
                        double hB213 = sampleB213Pt1[k] + coeffC * r2;
                        double hB223 = sampleB223Pt2[k] - coeffC * r1;

                        // rt
                        B2s[0, 2][k] = hB213;
                        B2s[2, 0][k] = hB213;
                        // st
                        B2s[1, 2][k] = hB223;
                        B2s[2, 1][k] = hB223;
                    }

                    for (int k = 0; k < nodeCnt * dof; k++)
                    {
                        for (int l = 0; l < nodeCnt * dof; l++)
                        {
                            double coeffC = sampleBu213Pt3[k, l] - sampleBu213Pt1[k, l] -
                                (sampleBu223Pt3[k, l] - sampleBu223Pt2[k, l]);
                            // rt
                            double hBu213 = sampleBu213Pt1[k, l] + coeffC * r2;
                            double hBu223 = sampleBu223Pt2[k, l] - coeffC * r1;

                            // rt
                            Bu2s[0, 2][k, l] = hBu213;
                            Bu2s[2, 0][k, l] = hBu213;
                            // st
                            Bu2s[1, 2][k, l] = hBu223;
                            Bu2s[2, 1][k, l] = hBu223;
                        }
                    }
                    {
                        double coeffC = sampleCgh13Pt3 - sampleCgh13Pt1 -
                            (sampleCgh23Pt3 - sampleCgh23Pt2);
                        // rt
                        double hCgh13 = sampleCgh13Pt1 + coeffC * r2;
                        double hCgh23 = sampleCgh23Pt2 - coeffC * r1;

                        // rt
                        cgh[0, 2] = hCgh13;
                        cgh[2, 0] = hCgh13;
                        // st
                        cgh[1, 2] = hCgh23;
                        cgh[2, 1] = hCgh23;
                    }
                    //-------------------------------------

                    //-------------------------------------
                    OpenTK.Vector3d[] h0s = CalcMITCThicknessStretchContravariantBasisVectors(g0s);
                    IvyFEM.Lapack.DoubleMatrix curS;
                    double[,,,] C4;
                    CalcMITCThicknessStretchStVenantC4StressAndTensor(
                        lambda, mu,
                        g0s, h0s,
                        cgh,
                        out curS, out C4);

                    //-------------------------------------
                    // Q
                    //-------------------------------------
                    // Q0
                    for (int p = 0; p < (nodeCnt * dof + additionalDof); p++)
                    {
                        double tmp = 0.0;
                        for (int ig = 0; ig < 3; ig++)
                        {
                            for (int ih = 0; ih < 3; ih++)
                            {
                                tmp += B1s[ig, ih][p] * curS[ig, ih];
                            }
                        }
                        localQe0[p] += detJWeightZPt * detJWeight * tmp;
                    }
                    //-------------------------------------
                    // Q2
                    for (int p = 0; p < (nodeCnt * dof + additionalDof); p++)
                    {
                        double tmp = 0.0;
                        for (int ig = 0; ig < 3; ig++)
                        {
                            for (int ih = 0; ih < 3; ih++)
                            {
                                tmp += curS[ig, ih] * B2s[ig, ih][p];
                            }
                        }
                        localQe2[p] += detJWeightZPt * detJWeight * tmp;
                    }
                    //-------------------------------------

                    //-------------------------------------
                    // K
                    //-------------------------------------
                    // K1
                    for (int p = 0; p < (nodeCnt * dof + additionalDof); p++)
                    {
                        for (int q = 0; q < (nodeCnt * dof + additionalDof); q++)
                        {
                            double tmp = 0.0;
                            for (int ig = 0; ig < 3; ig++)
                            {
                                for (int ih = 0; ih < 3; ih++)
                                {
                                    for (int e = 0; e < 3; e++)
                                    {
                                        for (int f = 0; f < 3; f++)
                                        {
                                            tmp += B1s[ig, ih][p] *
                                                C4[ig, ih, e, f] * B1s[e, f][q];
                                        }
                                    }
                                }
                            }
                            localKe1[p, q] += detJWeightZPt * detJWeight * tmp;
                        }
                    }

                    // K2
                    for (int p = 0; p < (nodeCnt * dof + additionalDof); p++)
                    {
                        for (int q = 0; q < (nodeCnt * dof + additionalDof); q++)
                        {
                            double tmp = 0.0;
                            for (int ig = 0; ig < 3; ig++)
                            {
                                for (int ih = 0; ih < 3; ih++)
                                {
                                    tmp += curS[ig, ih] * Bu2s[ig, ih][p, q];
                                }
                            }
                            localKe2[p, q] += detJWeightZPt * detJWeight * tmp;
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
                    double diagVal = localKe1[iNode * dof + k, iNode * dof + k];
                    //double diagVal = Math.Abs(localKe1[iNode * dof + k, iNode * dof + k]); //DEBUG
                    System.Diagnostics.Debug.Assert(diagVal >= IvyFEM.Constants.PrecisionLowerLimit);
                    if (diagVal > f)
                    {
                        f = diagVal;
                    }
                }
                f *= 1.0e-3;
                localKe1[iNode * dof + zOffset, iNode * dof + zOffset] = f;
            }
            for (int iNode = 0; iNode < 3; iNode++)
            {
                // fictitious stiffness value
                double f = 0.0;
                for (int k = 0; k < (uDof + aDof + bDof); k++)
                {
                    double diagVal = localKe2[iNode * dof + k, iNode * dof + k];
                    //double diagVal = Math.Abs(localKe2[iNode * dof + k, iNode * dof + k]); //DEBUG
                    //System.Diagnostics.Debug.Assert(diagVal >= IvyFEM.Constants.PrecisionLowerLimit);
                    if (diagVal > f)
                    {
                        f = diagVal;
                    }
                }
                f *= 1.0e-3;
                localKe2[iNode * dof + zOffset, iNode * dof + zOffset] = f;
            }

            //DEBUG
            for (int i = 0; i < (nodeCnt * dof + additionalDof); i++)
            {
                System.Diagnostics.Debug.Assert(Math.Abs(localKe1[i, i]) >= Constants.PrecisionLowerLimit);
                //System.Diagnostics.Debug.Assert(Math.Abs(localKe2[i, i]) >= Constants.PrecisionLowerLimit); // lambdaの対角項は0
            }

            //-------------------------------------------
            // global
            Qe0 = transT * localQe0;
            Qe2 = transT * localQe2;
            Ke1 = transT * localKe1 * T;
            Ke2 = transT * localKe2 * T;

            //DEBUG
            for (int i = 0; i < (nodeCnt * dof + additionalDof); i++)
            {
                System.Diagnostics.Debug.Assert(Math.Abs(Ke1[i, i]) >= Constants.PrecisionLowerLimit);
                //System.Diagnostics.Debug.Assert(Math.Abs(Ke2[i, i]) >= Constants.PrecisionLowerLimit); // lambdaの対角項は0
            }
        }
    }
}
