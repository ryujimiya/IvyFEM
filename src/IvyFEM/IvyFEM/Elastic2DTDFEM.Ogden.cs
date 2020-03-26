using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DTDFEM
    {
        protected void CalcOgdenRivlinHyperelasticElementAB(
            uint feId, IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            {
                uint quantityId0 = 0;
                TriangleFE workTriFE = World.GetTriangleFE(quantityId0, feId);
                Material workMa0 = World.GetMaterial(workTriFE.MaterialId);
                if (!(workMa0 is OgdenHyperelasticMaterial))
                {
                    return;
                }
            }

            uint uQuantityId = 0;
            uint lQuantityId = 1;
            System.Diagnostics.Debug.Assert(World.GetDof(uQuantityId) == 2);
            System.Diagnostics.Debug.Assert(World.GetDof(lQuantityId) == 1);
            int uDof = 2;
            int lDof = 1;
            int uNodeCnt = (int)World.GetNodeCount(uQuantityId);
            int lNodeCnt = (int)World.GetNodeCount(lQuantityId);
            //System.Diagnostics.Debug.Assert(uNodeCnt * uDof + lNodeCnt * lDof == A.RowLength);
            int offset = GetOffset(lQuantityId);
            System.Diagnostics.Debug.Assert(offset == uNodeCnt * uDof);

            double dt = TimeStep;
            double beta = NewmarkBeta;
            double gamma = NewmarkGamma;
            var uFV = World.GetFieldValue(UValueId);

            TriangleFE uTriFE = World.GetTriangleFE(uQuantityId, feId);
            TriangleFE lTriFE = World.GetTriangleFE(lQuantityId, feId);
            Material ma0 = World.GetMaterial(uTriFE.MaterialId);
            System.Diagnostics.Debug.Assert(ma0 is OgdenHyperelasticMaterial);

            uint vertexCnt = uTriFE.VertexCount;
            for (int iVertex = 0; iVertex < vertexCnt; iVertex++)
            {
                System.Diagnostics.Debug.Assert(uTriFE.VertexCoordIds[iVertex] == lTriFE.VertexCoordIds[iVertex]);
            }
            int[] uCoIds = uTriFE.NodeCoordIds;
            uint uElemNodeCnt = uTriFE.NodeCount;
            int[] uNodes = new int[uElemNodeCnt];
            for (int iNode = 0; iNode < uElemNodeCnt; iNode++)
            {
                int coId = uCoIds[iNode];
                int nodeId = World.Coord2Node(uQuantityId, coId);
                uNodes[iNode] = nodeId;
            }
            int[] lCoIds = lTriFE.NodeCoordIds;
            uint lElemNodeCnt = lTriFE.NodeCount;
            int[] lNodes = new int[lElemNodeCnt];
            for (int iNode = 0; iNode < lElemNodeCnt; iNode++)
            {
                int coId = lCoIds[iNode];
                int nodeId = World.Coord2Node(lQuantityId, coId);
                lNodes[iNode] = nodeId;
            }

            var ma = ma0 as OgdenHyperelasticMaterial;
            bool isCompressible = ma.IsCompressible;
            double d1 = ma.D1;
            int oOrder = ma.Order;
            double[] oMus = ma.Mus;
            double[] oAlphas = ma.Alphas;
            double rho = ma.MassDensity;
            double[] g = { ma.GravityX, ma.GravityY };

            double[] uSN = uTriFE.CalcSN();
            double[,] uSNN = uTriFE.CalcSNN();
            double[,] lSNN = lTriFE.CalcSNN();
            System.Diagnostics.Debug.Assert((int)World.TriIntegrationPointCount >= 3);
            IntegrationPoints ip = TriangleFE.GetIntegrationPoints(World.TriIntegrationPointCount);
            System.Diagnostics.Debug.Assert(ip.Ls.Length == (int)World.TriIntegrationPointCount);

            double[,] qu = new double[uElemNodeCnt, uDof];
            double[] ql = new double[lElemNodeCnt];

            for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
            {
                double[] L = ip.Ls[ipPt];
                double[] uN = uTriFE.CalcN(L);
                double[] lN = lTriFE.CalcN(L);
                double[][] uNu = uTriFE.CalcNu(L);
                double detJ = uTriFE.GetDetJacobian(L);
                double weight = ip.Weights[ipPt];
                double detJWeight = (1.0 / 2.0) * weight * detJ;

                // 変位の微分
                double[,] uu = new double[uDof, uDof];
                for (int iDof = 0; iDof < uDof; iDof++)
                {
                    for (int jDof = 0; jDof < uDof; jDof++)
                    {
                        for (int iNode = 0; iNode < uElemNodeCnt; iNode++)
                        {
                            int iNodeId = uNodes[iNode];
                            if (iNodeId == -1)
                            {
                                continue;
                            }
                            uu[iDof, jDof] += U[iNodeId * uDof + iDof] * uNu[jDof][iNode];
                        }
                    }
                }

                // ラグランジュの未定乗数
                double l = 0;
                for (int iNode = 0; iNode < lElemNodeCnt; iNode++)
                {
                    int iNodeId = lNodes[iNode];
                    if (iNodeId == -1)
                    {
                        continue;
                    }
                    l += U[offset + iNodeId] * lN[iNode];
                }

                // Green-Lagrangeのひずみ
                double[,] e = new double[uDof, uDof];
                for (int iDof = 0; iDof < uDof; iDof++)
                {
                    for (int jDof = 0; jDof < uDof; jDof++)
                    {
                        e[iDof, jDof] = (1.0 / 2.0) * (uu[iDof, jDof] + uu[jDof, iDof]);
                        for (int kDof = 0; kDof < uDof; kDof++)
                        {
                            e[iDof, jDof] += (1.0 / 2.0) * uu[kDof, iDof] * uu[kDof, jDof];
                        }
                    }
                }

                // 変形勾配テンソル
                double[,] f = new double[uDof, uDof];
                for (int iDof = 0; iDof < uDof; iDof++)
                {
                    for (int jDof = 0; jDof < uDof; jDof++)
                    {
                        f[iDof, jDof] = uu[iDof, jDof];
                    }
                    f[iDof, iDof] += 1.0;
                }

                // 右Cauchy-Green変形テンソル
                double[,] c = new double[uDof, uDof];
                for (int iDof = 0; iDof < uDof; iDof++)
                {
                    for (int jDof = 0; jDof < uDof; jDof++)
                    {
                        c[iDof, jDof] = uu[iDof, jDof] + uu[jDof, iDof];
                        for (int kDof = 0; kDof < uDof; kDof++)
                        {
                            c[iDof, jDof] += uu[kDof, iDof] * uu[kDof, jDof];
                        }
                    }
                    c[iDof, iDof] += 1.0;
                }

                // Cのテンソル不変量
                double I1 = c[0, 0] + c[1, 1] + 1.0;
                double I2 = c[0, 0] * c[1, 1] + c[0, 0] + c[1, 1] - c[0, 1] * c[1, 0];
                double I3 = c[0, 0] * c[1, 1] - c[0, 1] * c[1, 0];
                double inv13I3 = 1.0 / Math.Pow(I3, 1.0 / 3.0);
                double inv23I3 = 1.0 / Math.Pow(I3, 2.0 / 3.0);
                double[,] invC = new double[uDof, uDof];
                {
                    double invI3 = 1.0 / I3;
                    invC[0, 0] = invI3 * c[1, 1];
                    invC[0, 1] = -invI3 * c[0, 1];
                    invC[1, 0] = -invI3 * c[1, 0];
                    invC[1, 1] = invI3 * c[0, 0];
                }

                // Cの主値の1/2乗と主軸
                // Note: これらは3次元
                int dim3 = 3;
                System.Numerics.Complex[] fLambdas;
                System.Numerics.Complex[][] cNormals;
                SolvePrincipalValues(c, out fLambdas, out cNormals);
                double inv16I3 = 1.0 / Math.Pow(I3, 1.0 / 6.0);
                // 修正Cの主値
                System.Numerics.Complex[] barFLambdas = fLambdas.Select(a => a * inv16I3).ToArray();

                // 主第2Piola-Kirchhoff応力
                System.Numerics.Complex[] principalS = new System.Numerics.Complex[uDof];
                /*
                for (int kDof = 0; kDof < uDof; kDof++)
                {
                    for (int r = 0; r < oOrder; r++)
                    {
                        double mu = oMus[r];
                        double alpha = oAlphas[r];
                        principalS[kDof] += (mu / (fLambdas[kDof] * fLambdas[kDof])) * (
                            System.Numerics.Complex.Pow(barFLambdas[kDof], alpha) -
                            (1.0 / 3.0) * (
                            System.Numerics.Complex.Pow(barFLambdas[0], alpha) +
                            System.Numerics.Complex.Pow(barFLambdas[1], alpha) +
                            System.Numerics.Complex.Pow(barFLambdas[2], alpha)
                            ));
                    }
                }
                */
                for (int kDof = 0; kDof < uDof; kDof++)
                {
                    for (int r = 0; r < oOrder; r++)
                    {
                        double mu = oMus[r];
                        double alpha = oAlphas[r];
                        principalS[kDof] += (mu / (fLambdas[kDof] * fLambdas[kDof])) * (
                            System.Numerics.Complex.Pow(barFLambdas[kDof], alpha) -
                            (1.0 / 3.0) * (
                            System.Numerics.Complex.Pow(barFLambdas[0], alpha) +
                            System.Numerics.Complex.Pow(barFLambdas[1], alpha) +
                            Math.Pow(I3, -alpha / 6.0)
                            ));
                    }
                }

                // 第2Piola-Kirchhoff応力
                double[,] s = new double[uDof, uDof];
                for (int iDof = 0; iDof < uDof; iDof++)
                {
                    for (int jDof = 0; jDof < uDof; jDof++)
                    {
                        for (int kDof = 0; kDof < uDof; kDof++)
                        {
                            System.Numerics.Complex sij =
                                principalS[kDof] * cNormals[kDof][iDof] * cNormals[kDof][jDof];
                            s[iDof, jDof] += sij.Real;
                            //System.Diagnostics.Debug.Assert(
                            //    Math.Abs(sij.Imaginary) < IvyFEM.Constants.PrecisionLowerLimit);
                        }
                    }
                }
                // Lagrangeの未定乗数の寄与項
                {
                    double tmp = 2.0 * l * I3;
                    for (int iDof = 0; iDof < uDof; iDof++)
                    {
                        for (int jDof = 0; jDof < uDof; jDof++)
                        {
                            s[iDof, jDof] += tmp * invC[iDof, jDof];
                        }
                    }
                }

                // 主軸構成則テンソル(オリジナルの方)
                System.Numerics.Complex[,] principalC21st = new System.Numerics.Complex[uDof, uDof];
                System.Numerics.Complex[,] principalC22nd = new System.Numerics.Complex[uDof, uDof];
                /*
                for (int pDof = 0; pDof < uDof; pDof++)
                {
                    for (int qDof = 0; qDof < uDof; qDof++)
                    {
                        for (int r = 0; r < oOrder; r++)
                        {
                            double mu = oMus[r];
                            double alpha = oAlphas[r];
                            System.Numerics.Complex tmp1 = 0;
                            if (pDof == qDof)
                            {
                                tmp1 = -2.0 * (
                                    System.Numerics.Complex.Pow(barFLambdas[qDof], alpha) -
                                    (1.0 / 3.0) * (
                                    System.Numerics.Complex.Pow(barFLambdas[0], alpha) +
                                    System.Numerics.Complex.Pow(barFLambdas[1], alpha) +
                                    System.Numerics.Complex.Pow(barFLambdas[2], alpha)
                                    ));
                            }
                            System.Numerics.Complex tmp2 = 0;
                            for (int mDof = 0; mDof < dim3; mDof++)
                            {
                                double dmp = mDof == pDof ? 1 : 0;
                                double dmq = mDof == qDof ? 1 : 0;
                                tmp2 += alpha *
                                    System.Numerics.Complex.Pow(barFLambdas[mDof], alpha) *
                                    (dmp - 1.0 / 3.0) *
                                    (dmq - 1.0 / 3.0);
                            }
                            principalC21st[pDof, qDof] +=
                                (mu / (fLambdas[pDof] * fLambdas[pDof] * fLambdas[qDof] * fLambdas[qDof])) *
                                (tmp1 + tmp2);
                        }
                    }
                }
                for (int pDof = 0; pDof < uDof; pDof++)
                {
                    for (int qDof = 0; qDof < uDof; qDof++)
                    {
                        for (int r = 0; r < oOrder; r++)
                        {
                            double mu = oMus[r];
                            double alpha = oAlphas[r];
                            System.Numerics.Complex tmp = 0;
                            if (pDof != qDof)
                            {
                                System.Numerics.Complex diffFLambda = fLambdas[pDof] - fLambdas[qDof];
                                if (diffFLambda.Magnitude >= IvyFEM.Constants.PrecisionLowerLimit)
                                {
                                    // λp != λq
                                    System.Numerics.Complex squareFLambdaP = fLambdas[pDof] * fLambdas[pDof];
                                    System.Numerics.Complex squareFLambdaQ = fLambdas[qDof] * fLambdas[qDof];
                                    tmp = (mu / (squareFLambdaP * squareFLambdaQ)) * (
                                        (
                                        squareFLambdaQ * System.Numerics.Complex.Pow(barFLambdas[pDof], alpha) -
                                        squareFLambdaP * System.Numerics.Complex.Pow(barFLambdas[qDof], alpha)
                                        ) / (squareFLambdaP - squareFLambdaQ) +
                                        (1.0 / 3.0) *
                                        (
                                        System.Numerics.Complex.Pow(barFLambdas[0], alpha) +
                                        System.Numerics.Complex.Pow(barFLambdas[1], alpha) +
                                        System.Numerics.Complex.Pow(barFLambdas[2], alpha)
                                        ));
                                }
                                else
                                {
                                    // λp == λq
                                    tmp = (mu /
                                        (fLambdas[pDof] * fLambdas[pDof] * fLambdas[pDof] * fLambdas[pDof])) * (
                                        -((2.0 - alpha) / 2.0) *
                                        System.Numerics.Complex.Pow(barFLambdas[pDof], alpha) +
                                        (1.0 / 3.0) * (
                                        System.Numerics.Complex.Pow(barFLambdas[0], alpha) +
                                        System.Numerics.Complex.Pow(barFLambdas[1], alpha) +
                                        System.Numerics.Complex.Pow(barFLambdas[2], alpha)
                                        ));
                                }
                                principalC22nd[pDof, qDof] += tmp;
                            }
                        }
                    }
                }
                */
                uint dim2 = 2;
                for (int pDof = 0; pDof < uDof; pDof++)
                {
                    for (int qDof = 0; qDof < uDof; qDof++)
                    {
                        for (int r = 0; r < oOrder; r++)
                        {
                            double mu = oMus[r];
                            double alpha = oAlphas[r];
                            System.Numerics.Complex tmp1 = 0;
                            if (pDof == qDof)
                            {
                                tmp1 = -2.0 * (
                                    System.Numerics.Complex.Pow(barFLambdas[qDof], alpha) -
                                    (1.0 / 3.0) * (
                                    System.Numerics.Complex.Pow(barFLambdas[0], alpha) +
                                    System.Numerics.Complex.Pow(barFLambdas[1], alpha) +
                                    Math.Pow(I3, -alpha / 6.0)
                                    ));
                            }
                            System.Numerics.Complex tmp2 = 0;
                            for (int mDof = 0; mDof < dim2; mDof++)
                            {
                                double dmp = mDof == pDof ? 1 : 0;
                                double dmq = mDof == qDof ? 1 : 0;
                                tmp2 += alpha *
                                    System.Numerics.Complex.Pow(barFLambdas[mDof], alpha) *
                                    (dmp - 1.0 / 3.0) *
                                    (dmq - 1.0 / 3.0);
                            }
                            tmp2 += (alpha / 9.0) * Math.Pow(I3, -alpha / 6.0);
                            principalC21st[pDof, qDof] +=
                                (mu / (fLambdas[pDof] * fLambdas[pDof] * fLambdas[qDof] * fLambdas[qDof])) *
                                (tmp1 + tmp2);
                        }
                    }
                }
                for (int pDof = 0; pDof < uDof; pDof++)
                {
                    for (int qDof = 0; qDof < uDof; qDof++)
                    {
                        for (int r = 0; r < oOrder; r++)
                        {
                            double mu = oMus[r];
                            double alpha = oAlphas[r];
                            if (pDof != qDof)
                            {
                                System.Numerics.Complex tmp = 0;
                                System.Numerics.Complex diffFLambda = fLambdas[pDof] - fLambdas[qDof];
                                if (diffFLambda.Magnitude >= IvyFEM.Constants.PrecisionLowerLimit)
                                {
                                    // λp != λq
                                    System.Numerics.Complex squareFLambdaP = fLambdas[pDof] * fLambdas[pDof];
                                    System.Numerics.Complex squareFLambdaQ = fLambdas[qDof] * fLambdas[qDof];
                                    tmp = (mu / (squareFLambdaP * squareFLambdaQ)) * (
                                        (
                                        squareFLambdaQ * System.Numerics.Complex.Pow(barFLambdas[pDof], alpha) -
                                        squareFLambdaP * System.Numerics.Complex.Pow(barFLambdas[qDof], alpha)
                                        ) / (squareFLambdaP - squareFLambdaQ) +
                                        (1.0 / 3.0) *
                                        (
                                        System.Numerics.Complex.Pow(barFLambdas[0], alpha) +
                                        System.Numerics.Complex.Pow(barFLambdas[1], alpha) +
                                        System.Numerics.Complex.Pow(I3, -alpha / 6.0)
                                        ));
                                }
                                else
                                {
                                    // λp == λq
                                    tmp = (mu /
                                        (fLambdas[pDof] * fLambdas[pDof] * fLambdas[pDof] * fLambdas[pDof])) * (
                                        -((2.0 - alpha) / 2.0) *
                                        System.Numerics.Complex.Pow(barFLambdas[pDof], alpha) +
                                        (1.0 / 3.0) * (
                                        System.Numerics.Complex.Pow(barFLambdas[0], alpha) +
                                        System.Numerics.Complex.Pow(barFLambdas[1], alpha) +
                                        System.Numerics.Complex.Pow(I3, -alpha / 6.0)
                                        ));
                                }
                                principalC22nd[pDof, qDof] += tmp;
                            }
                        }
                    }
                }

                // 構成則テンソル（オリジナルの方）
                double[,,,] c4 = new double[uDof, uDof, uDof, uDof];
                for (int gDof = 0; gDof < uDof; gDof++)
                {
                    for (int hDof = 0; hDof < uDof; hDof++)
                    {
                        for (int eDof = 0; eDof < uDof; eDof++)
                        {
                            for (int fDof = 0; fDof < uDof; fDof++)
                            {
                                for (int pDof = 0; pDof < uDof; pDof++)
                                {
                                    for (int qDof = 0; qDof < uDof; qDof++)
                                    {
                                        System.Numerics.Complex c4ghef =
                                            principalC21st[pDof, qDof] *
                                            cNormals[pDof][gDof] * cNormals[pDof][hDof] *
                                            cNormals[qDof][eDof] * cNormals[qDof][fDof] +
                                            principalC22nd[pDof, qDof] *
                                            cNormals[pDof][gDof] * cNormals[qDof][hDof] *
                                            (
                                            cNormals[pDof][eDof] * cNormals[qDof][fDof] +
                                            cNormals[qDof][eDof] * cNormals[pDof][fDof]
                                            );
                                        c4[gDof, hDof, eDof, fDof] += c4ghef.Real;
                                        //System.Diagnostics.Debug.Assert(
                                        //    Math.Abs(c4ghef.Imaginary) < IvyFEM.Constants.PrecisionLowerLimit);
                                    }
                                }
                            }
                        }
                    }
                }
                // 構成則テンソル
                double[,,,] barC4 = new double[uDof, uDof, uDof, uDof];
                // 圧力構成則テンソル
                for (int gDof = 0; gDof < uDof; gDof++)
                {
                    for (int hDof = 0; hDof < uDof; hDof++)
                    {
                        for (int eDof = 0; eDof < uDof; eDof++)
                        {
                            for (int fDof = 0; fDof < uDof; fDof++)
                            {
                                barC4[gDof, hDof, eDof, fDof] +=
                                    4.0 * l * I3 * invC[gDof, hDof] * invC[eDof, fDof] -
                                    2.0 * l * I3 * (
                                    invC[gDof, eDof] * invC[fDof, hDof] +
                                    invC[gDof, fDof] * invC[eDof, hDof]);
                            }
                        }
                    }
                }
                // 変位構成則テンソル
                for (int gDof = 0; gDof < uDof; gDof++)
                {
                    for (int hDof = 0; hDof < uDof; hDof++)
                    {
                        for (int eDof = 0; eDof < uDof; eDof++)
                        {
                            for (int fDof = 0; fDof < uDof; fDof++)
                            {
                                barC4[gDof, hDof, eDof, fDof] +=
                                    (1.0 / 2.0) * (c4[gDof, hDof, eDof, fDof] + c4[gDof, hDof, fDof, eDof]);
                            }
                        }
                    }
                }

                double[,,,] b = new double[uElemNodeCnt, uDof, uDof, uDof];
                {
                    for (int iNode = 0; iNode < uElemNodeCnt; iNode++)
                    {
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            for (int gDof = 0; gDof < uDof; gDof++)
                            {
                                for (int hDof = 0; hDof < uDof; hDof++)
                                {
                                    b[iNode, iDof, gDof, hDof] = uNu[hDof][iNode] * f[iDof, gDof];
                                }
                            }
                        }
                    }
                }

                for (int row = 0; row < uElemNodeCnt; row++)
                {
                    int rowNodeId = uNodes[row];
                    if (rowNodeId == -1)
                    {
                        continue;
                    }
                    for (int col = 0; col < uElemNodeCnt; col++)
                    {
                        int colNodeId = uNodes[col];
                        if (colNodeId == -1)
                        {
                            continue;
                        }

                        double[,] kuu = new double[uDof, uDof];
                        for (int rowDof = 0; rowDof < uDof; rowDof++)
                        {
                            for (int colDof = 0; colDof < uDof; colDof++)
                            {
                                double tmp1 = 0.0;
                                for (int gDof = 0; gDof < uDof; gDof++)
                                {
                                    for (int hDof = 0; hDof < uDof; hDof++)
                                    {
                                        for (int eDof = 0; eDof < uDof; eDof++)
                                        {
                                            for (int fDof = 0; fDof < uDof; fDof++)
                                            {
                                                tmp1 += barC4[gDof, hDof, eDof, fDof] *
                                                    b[row, rowDof, eDof, fDof] *
                                                    b[col, colDof, gDof, hDof];
                                            }
                                        }
                                    }
                                }
                                kuu[rowDof, colDof] += detJWeight * tmp1;
                            }
                        }
                        {
                            double tmp = 0.0;
                            for (int gDof = 0; gDof < uDof; gDof++)
                            {
                                for (int hDof = 0; hDof < uDof; hDof++)
                                {
                                    tmp += s[gDof, hDof] * uNu[gDof][row] * uNu[hDof][col];
                                }
                            }
                            for (int rowDof = 0; rowDof < uDof; rowDof++)
                            {
                                kuu[rowDof, rowDof] += detJWeight * tmp;
                            }
                        }

                        for (int rowDof = 0; rowDof < uDof; rowDof++)
                        {
                            for (int colDof = 0; colDof < uDof; colDof++)
                            {
                                A[rowNodeId * uDof + rowDof, colNodeId * uDof + colDof] +=
                                    kuu[rowDof, colDof];
                                B[rowNodeId * uDof + rowDof] +=
                                    kuu[rowDof, colDof] * U[colNodeId * uDof + colDof];
                            }
                        }
                    }
                }
                for (int row = 0; row < uElemNodeCnt; row++)
                {
                    int rowNodeId = uNodes[row];
                    if (rowNodeId == -1)
                    {
                        continue;
                    }
                    for (int col = 0; col < lElemNodeCnt; col++)
                    {
                        int colNodeId = lNodes[col];
                        if (colNodeId == -1)
                        {
                            continue;
                        }

                        double[,] kul = new double[uDof, lDof];
                        double[,] klu = new double[lDof, uDof];
                        for (int rowDof = 0; rowDof < uDof; rowDof++)
                        {
                            double tmp = 0.0;
                            for (int gDof = 0; gDof < uDof; gDof++)
                            {
                                for (int hDof = 0; hDof < uDof; hDof++)
                                {
                                    tmp += invC[gDof, hDof] * b[row, rowDof, gDof, hDof];
                                }
                            }

                            kul[rowDof, 0] +=
                                detJWeight * tmp * 2.0 * lN[col] * I3;
                            klu[0, rowDof] +=
                                detJWeight * tmp * 2.0 * lN[col] * I3;
                        }

                        for (int rowDof = 0; rowDof < uDof; rowDof++)
                        {
                            A[rowNodeId * uDof + rowDof, offset + colNodeId] += kul[rowDof, 0];
                            A[offset + colNodeId, rowNodeId * uDof + rowDof] += klu[0, rowDof];
                            B[rowNodeId * uDof + rowDof] +=
                                kul[rowDof, 0] * U[offset + colNodeId];
                            B[offset + colNodeId] +=
                                klu[0, rowDof] * U[rowNodeId * uDof + rowDof];
                        }
                    }
                }

                // Note: kllは数値積分しなくて求められる

                for (int iNode = 0; iNode < uElemNodeCnt; iNode++)
                {
                    for (int iDof = 0; iDof < uDof; iDof++)
                    {
                        for (int gDof = 0; gDof < uDof; gDof++)
                        {
                            for (int hDof = 0; hDof < uDof; hDof++)
                            {
                                qu[iNode, iDof] += detJWeight * s[gDof, hDof] * b[iNode, iDof, gDof, hDof];
                            }
                        }
                    }
                }
                for (int iNode = 0; iNode < lElemNodeCnt; iNode++)
                {
                    if (isCompressible)
                    {
                        ql[iNode] += detJWeight * lN[iNode] * ((I3 - 1) - l / d1);
                    }
                    else
                    {
                        ql[iNode] += detJWeight * lN[iNode] * (I3 - 1);
                    }
                }
            }

            for (int row = 0; row < uElemNodeCnt; row++)
            {
                int rowNodeId = uNodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                for (int col = 0; col < uElemNodeCnt; col++)
                {
                    int colNodeId = uNodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }

                    int colCoId = uCoIds[col];
                    double[] u = uFV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                    double[] velU = uFV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                    double[] accU = uFV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);
                    double[,] m = new double[2, 2];

                    m[0, 0] = rho * uSNN[row, col];
                    m[1, 0] = 0.0;
                    m[0, 1] = 0.0;
                    m[1, 1] = rho * uSNN[row, col];
                    for (int rowDof = 0; rowDof < uDof; rowDof++)
                    {
                        for (int colDof = 0; colDof < uDof; colDof++)
                        {
                            A[rowNodeId * uDof + rowDof, colNodeId * uDof + colDof] +=
                                (1.0 / (beta * dt * dt)) * m[rowDof, colDof];
                            B[rowNodeId * uDof + rowDof] +=
                                m[rowDof, colDof] * (
                                (1.0 / (beta * dt * dt)) * u[colDof] +
                                (1.0 / (beta * dt)) * velU[colDof] +
                                (1.0 / (2.0 * beta) - 1.0) * accU[colDof]);
                        }
                    }
                }
            }

            if (isCompressible)
            {
                for (int row = 0; row < lElemNodeCnt; row++)
                {
                    int rowNodeId = lNodes[row];
                    if (rowNodeId == -1)
                    {
                        continue;
                    }
                    for (int col = 0; col < lElemNodeCnt; col++)
                    {
                        int colNodeId = lNodes[col];
                        if (colNodeId == -1)
                        {
                            continue;
                        }

                        double kll = -(1.0 / d1) * lSNN[row, col];
                        A[offset + rowNodeId, offset + colNodeId] += kll;
                        B[offset + rowNodeId] += kll * U[offset + colNodeId];
                    }
                }
            }

            double[,] fg = new double[uElemNodeCnt, uDof];
            for (int iNode = 0; iNode < uElemNodeCnt; iNode++)
            {
                for (int iDof = 0; iDof < uDof; iDof++)
                {
                    fg[iNode, iDof] = rho * g[iDof] * uSN[iNode];
                }
            }

            for (int row = 0; row < uElemNodeCnt; row++)
            {
                int rowNodeId = uNodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                for (int rowDof = 0; rowDof < uDof; rowDof++)
                {
                    B[rowNodeId * uDof + rowDof] += fg[row, rowDof] - qu[row, rowDof];
                }
            }
            for (int row = 0; row < lElemNodeCnt; row++)
            {
                int rowNodeId = lNodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                B[offset + rowNodeId] += -ql[row];
            }
        }
    }
}
