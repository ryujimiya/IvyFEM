using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Fluid2DFEM
    {
        private void CalcSUPGPressurePoissonWithBellAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            int quantityCnt = World.GetQuantityCount();
            System.Diagnostics.Debug.Assert(quantityCnt == 7);
            uint vQuantityId = 0;
            uint pQuantityId = 1;
            uint pxQuantityId = 2;
            uint pyQuantityId = 3;
            uint pxxQuantityId = 4;
            uint pxyQuantityId = 5;
            uint pyyQuantityId = 6;
            int vDof = 2;
            int pDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int pxNodeCnt = (int)World.GetNodeCount(pxQuantityId);
            int pyNodeCnt = (int)World.GetNodeCount(pyQuantityId);
            int pxxNodeCnt = (int)World.GetNodeCount(pxxQuantityId);
            int pxyNodeCnt = (int)World.GetNodeCount(pxyQuantityId);
            int pyyNodeCnt = (int)World.GetNodeCount(pyyQuantityId);
            int offsetp = vNodeCnt * vDof;
            int offsetpx = offsetp + pNodeCnt;
            int offsetpy = offsetpx + pxNodeCnt;
            int offsetpxx = offsetpy + pyNodeCnt;
            int offsetpxy = offsetpxx + pxxNodeCnt;
            int offsetpyy = offsetpxy + pxyNodeCnt;

            IList<uint> feIds = World.GetTriangleFEIds(vQuantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE vTriFE = World.GetTriangleFE(vQuantityId, feId);
                TriangleFE pTriFE = World.GetTriangleFE(pQuantityId, feId);
                TriangleFE pxTriFE = World.GetTriangleFE(pxQuantityId, feId);
                TriangleFE pyTriFE = World.GetTriangleFE(pyQuantityId, feId);
                TriangleFE pxxTriFE = World.GetTriangleFE(pxxQuantityId, feId);
                TriangleFE pxyTriFE = World.GetTriangleFE(pxyQuantityId, feId);
                TriangleFE pyyTriFE = World.GetTriangleFE(pyyQuantityId, feId);
                uint vertexCnt = vTriFE.VertexCount;
                for (int iVertex = 0; iVertex < vertexCnt; iVertex++)
                {
                    System.Diagnostics.Debug.Assert(vTriFE.VertexCoordIds[iVertex] == pTriFE.VertexCoordIds[iVertex]);
                }

                int[] vCoIds = vTriFE.NodeCoordIds;
                uint vElemNodeCnt = vTriFE.NodeCount;
                int[] vNodes = new int[vElemNodeCnt];
                for (int iNode = 0; iNode < vElemNodeCnt; iNode++)
                {
                    int coId = vCoIds[iNode];
                    int nodeId = World.Coord2Node(vQuantityId, coId);
                    vNodes[iNode] = nodeId;
                }
                int[] pCoIds = pTriFE.NodeCoordIds;
                uint pElemNodeCnt = pTriFE.NodeCount;
                int[] pNodes = new int[pElemNodeCnt];
                for (int iNode = 0; iNode < pElemNodeCnt; iNode++)
                {
                    int coId = pCoIds[iNode];
                    int nodeId = World.Coord2Node(pQuantityId, coId);
                    pNodes[iNode] = nodeId;
                }
                int[] pxCoIds = pxTriFE.NodeCoordIds;
                uint pxElemNodeCnt = pxTriFE.NodeCount;
                int[] pxNodes = new int[pxElemNodeCnt];
                for (int iNode = 0; iNode < pxElemNodeCnt; iNode++)
                {
                    int coId = pxCoIds[iNode];
                    int nodeId = World.Coord2Node(pxQuantityId, coId);
                    pxNodes[iNode] = nodeId;
                }
                int[] pyCoIds = pyTriFE.NodeCoordIds;
                uint pyElemNodeCnt = pyTriFE.NodeCount;
                int[] pyNodes = new int[pyElemNodeCnt];
                for (int iNode = 0; iNode < pyElemNodeCnt; iNode++)
                {
                    int coId = pyCoIds[iNode];
                    int nodeId = World.Coord2Node(pyQuantityId, coId);
                    pyNodes[iNode] = nodeId;
                }
                int[] pxxCoIds = pxxTriFE.NodeCoordIds;
                uint pxxElemNodeCnt = pxxTriFE.NodeCount;
                int[] pxxNodes = new int[pxxElemNodeCnt];
                for (int iNode = 0; iNode < pxxElemNodeCnt; iNode++)
                {
                    int coId = pxxCoIds[iNode];
                    int nodeId = World.Coord2Node(pxxQuantityId, coId);
                    pxxNodes[iNode] = nodeId;
                }
                int[] pxyCoIds = pxyTriFE.NodeCoordIds;
                uint pxyElemNodeCnt = pxyTriFE.NodeCount;
                int[] pxyNodes = new int[pxyElemNodeCnt];
                for (int iNode = 0; iNode < pxyElemNodeCnt; iNode++)
                {
                    int coId = pxyCoIds[iNode];
                    int nodeId = World.Coord2Node(pxyQuantityId, coId);
                    pxyNodes[iNode] = nodeId;
                }
                int[] pyyCoIds = pyyTriFE.NodeCoordIds;
                uint pyyElemNodeCnt = pyyTriFE.NodeCount;
                int[] pyyNodes = new int[pyyElemNodeCnt];
                for (int iNode = 0; iNode < pyyElemNodeCnt; iNode++)
                {
                    int coId = pyyCoIds[iNode];
                    int nodeId = World.Coord2Node(pyyQuantityId, coId);
                    pyyNodes[iNode] = nodeId;
                }

                Material ma0 = World.GetMaterial(vTriFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is NewtonFluidMaterial);
                var ma = ma0 as NewtonFluidMaterial;
                double rho = ma.MassDensity;
                double mu = ma.Mu;
                double nu = mu / rho;
                double[] g = { ma.GravityX, ma.GravityY };

                double[][] velos = new double[vElemNodeCnt][];
                for (int iNode = 0; iNode < vElemNodeCnt; iNode++)
                {
                    double[] velo = new double[vDof];
                    int nodeId = vNodes[iNode];
                    if (nodeId == -1)
                    {
                        // 0
                    }
                    else
                    {
                        for (int iDof = 0; iDof < vDof; iDof++)
                        {
                            velo[iDof] = U[nodeId * vDof + iDof];
                        }
                    }
                    velos[iNode] = velo;
                }

                double taum = 0;
                double tauc = 0;
                {
                    double[] aveVelo = {
                        (velos[0][0] + velos[1][0] + velos[2][0]) / 3.0,
                        (velos[0][1] + velos[1][1] + velos[2][1]) / 3.0
                    };
                    double veloNorm = Math.Sqrt(aveVelo[0] * aveVelo[0] + aveVelo[1] * aveVelo[1]);
                    double[][] Lu = new double[vDof][];
                    {
                        double[] a;
                        double[] b;
                        double[] c;
                        vTriFE.CalcTransMatrix(out a, out b, out c);
                        // Lx
                        Lu[0] = b;
                        // Ly
                        Lu[1] = c;
                    }
                    IvyFEM.Lapack.DoubleMatrix GMat = new IvyFEM.Lapack.DoubleMatrix(vDof, vDof);
                    double[] gVec = new double[vDof];
                    for (int iDof = 0; iDof < vDof; iDof++)
                    {
                        for (int jDof = 0; jDof < vDof; jDof++)
                        {
                            for (int kDof = 0; kDof < vDof; kDof++)
                            {
                                GMat[iDof, jDof] += Lu[iDof][kDof] * Lu[jDof][kDof];
                            }
                        }
                    }
                    for (int iDof = 0; iDof < vDof; iDof++)
                    {
                        for (int kDof = 0; kDof < vDof; kDof++)
                        {
                            gVec[iDof] += Lu[iDof][kDof];
                        }
                    }

                    double sqinvtaum1 = 0;
                    double sqinvtaum2 = 0;
                    {
                        double[] tmpVec = GMat * aveVelo;
                        sqinvtaum2 = IvyFEM.Lapack.Functions.ddot(aveVelo, tmpVec);
                    }
                    double sqinvtaum3 = 0;
                    {
                        IvyFEM.Lapack.DoubleMatrix GMatT = new Lapack.DoubleMatrix(GMat);
                        GMatT.Transpose();
                        double GMatDoubleDot = IvyFEM.Lapack.DoubleMatrix.DoubleDot(GMat, GMatT);
                        sqinvtaum3 = 30.0 * nu * nu * GMatDoubleDot;
                    }
                    double sqinvtaum = sqinvtaum1 + sqinvtaum2 + sqinvtaum3;
                    taum = 1.0 / Math.Sqrt(sqinvtaum);

                    double gDot = IvyFEM.Lapack.Functions.ddot(gVec, gVec);
                    tauc = 1.0 / (taum * gDot);
                }

                double[] vSN = vTriFE.CalcSN();
                IntegrationPoints ip = TriangleFE.GetIntegrationPoints(World.TriIntegrationPointCount);//Point7
                for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                {
                    double[] L = ip.Ls[ipPt];
                    double[] vN = vTriFE.CalcN(L);
                    double[][] vNu = vTriFE.CalcNu(L);
                    double[] vNx = vNu[0];
                    double[] vNy = vNu[1];
                    double[,][] vNuv = vTriFE.CalcNuv(L);
                    double[] pN = pTriFE.CalcN(L);
                    double[] pxN = pxTriFE.CalcN(L);
                    double[] pyN = pyTriFE.CalcN(L);
                    double[] pxxN = pxxTriFE.CalcN(L);
                    double[] pxyN = pxyTriFE.CalcN(L);
                    double[] pyyN = pyyTriFE.CalcN(L);
                    double[][] pNu = pTriFE.CalcNu(L);
                    double[] pNx = pNu[0];
                    double[] pNy = pNu[1];
                    double[][] pxNu = pxTriFE.CalcNu(L);
                    double[] pxNx = pxNu[0];
                    double[] pxNy = pxNu[1];
                    double[][] pyNu = pyTriFE.CalcNu(L);
                    double[] pyNx = pyNu[0];
                    double[] pyNy = pyNu[1];
                    double[][] pxxNu = pxxTriFE.CalcNu(L);
                    double[] pxxNx = pxxNu[0];
                    double[] pxxNy = pxxNu[1];
                    double[][] pxyNu = pxyTriFE.CalcNu(L);
                    double[] pxyNx = pxyNu[0];
                    double[] pxyNy = pxyNu[1];
                    double[][] pyyNu = pyyTriFE.CalcNu(L);
                    double[] pyyNx = pyyNu[0];
                    double[] pyyNy = pyyNu[1];
                    double[,][] pNuv = pTriFE.CalcNuv(L);
                    double[,][] pxNuv = pxTriFE.CalcNuv(L);
                    double[,][] pyNuv = pyTriFE.CalcNuv(L);
                    double[,][] pxxNuv = pxxTriFE.CalcNuv(L);
                    double[,][] pxyNuv = pxyTriFE.CalcNuv(L);
                    double[,][] pyyNuv = pyyTriFE.CalcNuv(L);
                    int[] offsetps = {
                        offsetp, offsetpx, offsetpy,
                        offsetpxx, offsetpxy, offsetpyy
                    };
                    uint[] pElemNodeCnts = {
                        pElemNodeCnt, pxElemNodeCnt, pyElemNodeCnt,
                        pxxElemNodeCnt, pxyElemNodeCnt, pyyElemNodeCnt
                    };
                    int[][] pNodess = {
                        pNodes, pxNodes, pyNodes,
                        pxxNodes, pxyNodes, pyyNodes
                    };
                    double[][] pNs = {
                        pN, pxN, pyN,
                        pxxN, pxyN, pyyN
                    };
                    double[][] pNxs = {
                        pNx, pxNx, pyNx,
                        pxxNx, pxyNx, pyyNx
                    };
                    double[][] pNys = {
                        pNy, pxNy, pyNy,
                        pxxNy, pxyNy, pyyNy
                    };

                    double detJ = vTriFE.GetDetJacobian(L);
                    double weight = ip.Weights[ipPt];
                    double detJWeight = (1.0 / 2.0) * weight * detJ;

                    double[] v = new double[vDof];
                    double[] vx = new double[vDof];
                    double[] vy = new double[vDof];
                    double[] vxx = new double[vDof];
                    double[] vxy = new double[vDof];
                    double[] vyx = new double[vDof];
                    double[] vyy = new double[vDof];
                    double p = 0;
                    double px = 0;
                    double py = 0;
                    double pxx = 0;
                    double pxy = 0;
                    double pyy = 0;
                    for (int iNode = 0; iNode < vElemNodeCnt; iNode++)
                    {
                        int nodeId = vNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }
                        for (int iDof = 0; iDof < vDof; iDof++)
                        {
                            double vValue = U[nodeId * vDof + iDof];
                            v[iDof] += vValue * vN[iNode];
                            vx[iDof] += vValue * vNx[iNode];
                            vy[iDof] += vValue * vNy[iNode];
                            vxx[iDof] += vValue * vNuv[0, 0][iNode];
                            vxy[iDof] += vValue * vNuv[0, 1][iNode];
                            vyx[iDof] += vValue * vNuv[1, 0][iNode];
                            vyy[iDof] += vValue * vNuv[1, 1][iNode];
                        }
                    }
                    // 6変数 x 3頂点 = 21 shape functions
                    for (int iNode = 0; iNode < pElemNodeCnt; iNode++)
                    {
                        int nodeId = pNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }
                        double pValue1 = U[offsetp + nodeId];
                        p += pValue1 * pN[iNode];
                        px += pValue1 * pNx[iNode];
                        py += pValue1 * pNy[iNode];
                        pxx += pValue1 * pNuv[0, 0][iNode];
                        pxy += pValue1 * pNuv[0, 1][iNode];
                        pyy += pValue1 * pNuv[1, 1][iNode];
                    }
                    for (int iNode = 0; iNode < pxElemNodeCnt; iNode++)
                    {
                        int nodeId = pxNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }
                        double pValue2 = U[offsetpx + nodeId];
                        p += pValue2 * pxN[iNode];
                        px += pValue2 * pxNx[iNode];
                        py += pValue2 * pxNy[iNode];
                        pxx += pValue2 * pxNuv[0, 0][iNode];
                        pxy += pValue2 * pxNuv[0, 1][iNode];
                        pyy += pValue2 * pxNuv[1, 1][iNode];
                    }
                    for (int iNode = 0; iNode < pyElemNodeCnt; iNode++)
                    {
                        int nodeId = pyNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }
                        double pValue3 = U[offsetpy + nodeId];
                        p += pValue3 * pyN[iNode];
                        px += pValue3 * pyNx[iNode];
                        py += pValue3 * pyNy[iNode];
                        pxx += pValue3 * pyNuv[0, 0][iNode];
                        pxy += pValue3 * pyNuv[0, 1][iNode];
                        pyy += pValue3 * pyNuv[1, 1][iNode];
                    }
                    for (int iNode = 0; iNode < pxxElemNodeCnt; iNode++)
                    {
                        int nodeId = pxxNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }
                        double pValue4 = U[offsetpxx + nodeId];
                        p += pValue4 * pxxN[iNode];
                        px += pValue4 * pxxNx[iNode];
                        py += pValue4 * pxxNy[iNode];
                        pxx += pValue4 * pxxNuv[0, 0][iNode];
                        pxy += pValue4 * pxxNuv[0, 1][iNode];
                        pyy += pValue4 * pxxNuv[1, 1][iNode];
                    }
                    for (int iNode = 0; iNode < pxyElemNodeCnt; iNode++)
                    {
                        int nodeId = pxyNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }
                        double pValue5 = U[offsetpxy + nodeId];
                        p += pValue5 * pxyN[iNode];
                        px += pValue5 * pxyNx[iNode];
                        py += pValue5 * pxyNy[iNode];
                        pxx += pValue5 * pxyNuv[0, 0][iNode];
                        pxy += pValue5 * pxyNuv[0, 1][iNode];
                        pyy += pValue5 * pxyNuv[1, 1][iNode];
                    }
                    for (int iNode = 0; iNode < pyyElemNodeCnt; iNode++)
                    {
                        int nodeId = pyyNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }
                        double pValue6 = U[offsetpyy + nodeId];
                        p += pValue6 * pyyN[iNode];
                        px += pValue6 * pyyNx[iNode];
                        py += pValue6 * pyyNy[iNode];
                        pxx += pValue6 * pyyNuv[0, 0][iNode];
                        pxy += pValue6 * pyyNuv[0, 1][iNode];
                        pyy += pValue6 * pyyNuv[1, 1][iNode];
                    }
                    double[][] vu = new double[vDof][];
                    vu[0] = vx;
                    vu[1] = vy;
                    double[,][] vuv = new double[vDof, vDof][];
                    vuv[0, 0] = vxx;
                    vuv[0, 1] = vxy;
                    vuv[1, 0] = vyx;
                    vuv[1, 1] = vyy;
                    double[] pu = new double[vDof];
                    pu[0] = px;
                    pu[1] = py;
                    double[,] puv = new double[vDof, vDof];
                    puv[0, 0] = pxx;
                    puv[0, 1] = pxy;
                    puv[1, 0] = pxy;
                    puv[1, 1] = pyy;
                    double[] gx = new double[vDof];
                    double[] gy = new double[vDof];
                    for (int iNode = 0; iNode < vElemNodeCnt; iNode++)
                    {
                        int nodeId = vNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }
                        for (int iDof = 0; iDof < vDof; iDof++)
                        {
                            double gValue = g[iDof];
                            gx[iDof] += gValue * vNx[iNode];
                            gy[iDof] += gValue * vNy[iNode];
                        }
                    }

                    for (int row = 0; row < vElemNodeCnt; row++)
                    {
                        int rowNodeId = vNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < vElemNodeCnt; col++)
                        {
                            int colNodeId = vNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double[,] kvv1 = new double[vDof, vDof];
                            kvv1[0, 0] = detJWeight * mu * (vNx[row] * vNx[col] +
                                vNx[row] * vNx[col] + vNy[row] * vNy[col]);
                            kvv1[0, 1] = detJWeight * mu * vNy[row] * vNx[col];
                            kvv1[1, 0] = detJWeight * mu * vNx[row] * vNy[col];
                            kvv1[1, 1] = detJWeight * mu * (vNy[row] * vNy[col] +
                                vNx[row] * vNx[col] + vNy[row] * vNy[col]);

                            double[,] kvv2 = new double[vDof, vDof];
                            kvv2[0, 0] = detJWeight * rho * vN[row] * (v[0] * vNx[col] + v[1] * vNy[col]);
                            kvv2[0, 1] = 0;
                            kvv2[1, 0] = 0;
                            kvv2[1, 1] = detJWeight * rho * vN[row] * (v[0] * vNx[col] + v[1] * vNy[col]);

                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                for (int colDof = 0; colDof < vDof; colDof++)
                                {
                                    A[rowNodeId * vDof + rowDof, colNodeId * vDof + colDof] +=
                                        kvv1[rowDof, colDof] + kvv2[rowDof, colDof];
                                }
                            }
                        }
                    }

                    for (int row = 0; row < vElemNodeCnt; row++)
                    {
                        int rowNodeId = vNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }

                        for (uint colQuantityId = pQuantityId; colQuantityId < quantityCnt; colQuantityId++)
                        {
                            int coloffsetp = offsetps[colQuantityId - pQuantityId];
                            uint colpElemNodeCnt = pElemNodeCnts[colQuantityId - pQuantityId];
                            int[] colpNodes = pNodess[colQuantityId - pQuantityId];
                            double[] colpN = pNs[colQuantityId - pQuantityId];
                            double[] colpNx = pNxs[colQuantityId - pQuantityId];
                            double[] colpNy = pNys[colQuantityId - pQuantityId];

                            for (int col = 0; col < colpElemNodeCnt; col++)
                            {
                                int colNodeId = colpNodes[col];
                                if (colNodeId == -1)
                                {
                                    continue;
                                }

                                double[,] kvp = new double[vDof, pDof];
                                kvp[0, 0] = -detJWeight * vNx[row] * colpN[col];
                                kvp[1, 0] = -detJWeight * vNy[row] * colpN[col];

                                for (int rowDof = 0; rowDof < vDof; rowDof++)
                                {
                                    A[rowNodeId * vDof + rowDof, coloffsetp + colNodeId] += kvp[rowDof, 0];
                                }
                            }
                        }
                    }

                    // Pressure Poisson
                    // kpv
                    for (uint rowQuantityId = pQuantityId; rowQuantityId < quantityCnt; rowQuantityId++)
                    {
                        int rowoffsetp = offsetps[rowQuantityId - pQuantityId];
                        uint rowpElemNodeCnt = pElemNodeCnts[rowQuantityId - pQuantityId];
                        int[] rowpNodes = pNodess[rowQuantityId - pQuantityId];
                        double[] rowpN = pNs[rowQuantityId - pQuantityId];
                        double[] rowpNx = pNxs[rowQuantityId - pQuantityId];
                        double[] rowpNy = pNys[rowQuantityId - pQuantityId];

                        for (int row = 0; row < rowpElemNodeCnt; row++)
                        {
                            int rowNodeId = rowpNodes[row];
                            if (rowNodeId == -1)
                            {
                                continue;
                            }
                            for (int col = 0; col < vElemNodeCnt; col++)
                            {
                                int colNodeId = vNodes[col];
                                if (colNodeId == -1)
                                {
                                    continue;
                                }

                                double[,] kpv = new double[pDof, vDof];
                                kpv[0, 0] = detJWeight * rho * rowpNx[row] * (v[0] * vNx[col] + v[1] * vNy[col]);
                                kpv[0, 1] = detJWeight * rho * rowpNy[row] * (v[0] * vNx[col] + v[1] * vNy[col]);

                                for (int colDof = 0; colDof < vDof; colDof++)
                                {
                                    A[rowoffsetp + rowNodeId, colNodeId * vDof + colDof] += kpv[0, colDof];
                                }
                            }
                        }
                    }
                    // kpp
                    for (uint rowQuantityId = pQuantityId; rowQuantityId < quantityCnt; rowQuantityId++)
                    {
                        int rowoffsetp = offsetps[rowQuantityId - pQuantityId];
                        uint rowpElemNodeCnt = pElemNodeCnts[rowQuantityId - pQuantityId];
                        int[] rowpNodes = pNodess[rowQuantityId - pQuantityId];
                        double[] rowpN = pNs[rowQuantityId - pQuantityId];
                        double[] rowpNx = pNxs[rowQuantityId - pQuantityId];
                        double[] rowpNy = pNys[rowQuantityId - pQuantityId];

                        for (int row = 0; row < rowpElemNodeCnt; row++)
                        {
                            int rowNodeId = rowpNodes[row];
                            if (rowNodeId == -1)
                            {
                                continue;
                            }

                            for (uint colQuantityId = pQuantityId; colQuantityId < quantityCnt; colQuantityId++)
                            {
                                int coloffsetp = offsetps[colQuantityId - pQuantityId];
                                uint colpElemNodeCnt = pElemNodeCnts[colQuantityId - pQuantityId];
                                int[] colpNodes = pNodess[colQuantityId - pQuantityId];
                                double[] colpN = pNs[colQuantityId - pQuantityId];
                                double[] colpNx = pNxs[colQuantityId - pQuantityId];
                                double[] colpNy = pNys[colQuantityId - pQuantityId];

                                for (int col = 0; col < colpElemNodeCnt; col++)
                                {
                                    int colNodeId = colpNodes[col];
                                    if (colNodeId == -1)
                                    {
                                        continue;
                                    }

                                    double[,] kpp = new double[pDof, pDof];
                                    kpp[0, 0] = detJWeight * (rowpNx[row] * colpNx[col] + rowpNy[row] * colpNy[col]);

                                    A[rowoffsetp + rowNodeId, coloffsetp + colNodeId] += kpp[0, 0];
                                }
                            }
                        }
                    }
                    // fp
                    for (uint rowQuantityId = pQuantityId; rowQuantityId < quantityCnt; rowQuantityId++)
                    {
                        int rowoffsetp = offsetps[rowQuantityId - pQuantityId];
                        uint rowpElemNodeCnt = pElemNodeCnts[rowQuantityId - pQuantityId];
                        int[] rowpNodes = pNodess[rowQuantityId - pQuantityId];
                        double[] rowpN = pNs[rowQuantityId - pQuantityId];
                        double[] rowpNx = pNxs[rowQuantityId - pQuantityId];
                        double[] rowpNy = pNys[rowQuantityId - pQuantityId];

                        for (int row = 0; row < rowpElemNodeCnt; row++)
                        {
                            int rowNodeId = rowpNodes[row];
                            if (rowNodeId == -1)
                            {
                                continue;
                            }
                            double fp = detJWeight * rho * (g[0] * rowpNx[row] + g[1] * rowpNy[row]);
                            double qp = 0;
                            B[rowoffsetp + rowNodeId] +=
                                fp - qp;
                        }
                    }

                    //////////////////////////////////////////////////////////////
                    // SUPG
                    double[] rmi = new double[vDof];
                    double[,,] rmivj = new double[vDof, vDof, vElemNodeCnt];
                    uint allpElemNodeCnt =
                        pElemNodeCnt +
                        pxElemNodeCnt + pyElemNodeCnt +
                        pxxElemNodeCnt + pxyElemNodeCnt + pyyElemNodeCnt;
                    uint[] pElemNodeOffsets = new uint[6] {
                        0,
                        pElemNodeCnt,
                        pElemNodeCnt + pxElemNodeCnt,
                        pElemNodeCnt + pxElemNodeCnt + pyElemNodeCnt,
                        pElemNodeCnt + pxElemNodeCnt + pyElemNodeCnt + pxxElemNodeCnt,
                        pElemNodeCnt + pxElemNodeCnt + pyElemNodeCnt + pxxElemNodeCnt + pxyElemNodeCnt
                    };
                    double[,] rmip = new double[vDof, allpElemNodeCnt];
                    double rc = 0;
                    double[,] rcvj = new double[vDof, vElemNodeCnt];
                    double[] rcp = new double[allpElemNodeCnt];
                    for (int iDof = 0; iDof < vDof; iDof++)
                    {
                        rmi[iDof] =
                            -mu * (vuv[0, 0][iDof] + vuv[1, 1][iDof]) +
                            rho * (v[0] * vx[iDof] + v[1] * vy[iDof]) +
                            pu[iDof] - rho * g[iDof];

                        for (int jDof = 0; jDof < vDof; jDof++)
                        {
                            for (int jNode = 0; jNode < vElemNodeCnt; jNode++)
                            {
                                int jNodeId = vNodes[jNode];
                                if (jNodeId == -1)
                                {
                                    continue;
                                }

                                rmivj[iDof, jDof, jNode] = 0;
                                if (iDof == jDof)
                                {
                                    rmivj[iDof, jDof, jNode] +=
                                        -mu * (vNuv[0, 0][jNode] + vNuv[1, 1][jNode]);
                                }
                                if (iDof == jDof)
                                {
                                    rmivj[iDof, jDof, jNode] +=
                                        rho * (v[0] * vNu[0][jNode] + v[1] * vNu[1][jNode]);
                                }
                            }
                        }

                        // p
                        for (int jNode = 0; jNode < pElemNodeCnt; jNode++)
                        {
                            int jNodeId = pNodes[jNode];
                            if (jNodeId == -1)
                            {
                                continue;
                            }
                            rmip[iDof, jNode + pElemNodeOffsets[0]] = pNu[iDof][jNode];
                        }
                        // px
                        for (int jNode = 0; jNode < pxElemNodeCnt; jNode++)
                        {
                            int jNodeId = pxNodes[jNode];
                            if (jNodeId == -1)
                            {
                                continue;
                            }
                            rmip[iDof, jNode + pElemNodeOffsets[1]] = pxNu[iDof][jNode];
                        }
                        // py
                        for (int jNode = 0; jNode < pyElemNodeCnt; jNode++)
                        {
                            int jNodeId = pyNodes[jNode];
                            if (jNodeId == -1)
                            {
                                continue;
                            }
                            rmip[iDof, jNode + pElemNodeOffsets[2]] = pyNu[iDof][jNode];
                        }
                        // pxx
                        for (int jNode = 0; jNode < pxxElemNodeCnt; jNode++)
                        {
                            int jNodeId = pxxNodes[jNode];
                            if (jNodeId == -1)
                            {
                                continue;
                            }
                            rmip[iDof, jNode + pElemNodeOffsets[3]] = pxxNu[iDof][jNode];
                        }
                        // pxy
                        for (int jNode = 0; jNode < pxyElemNodeCnt; jNode++)
                        {
                            int jNodeId = pxyNodes[jNode];
                            if (jNodeId == -1)
                            {
                                continue;
                            }
                            rmip[iDof, jNode + pElemNodeOffsets[4]] = pxyNu[iDof][jNode];
                        }
                        // pyy
                        for (int jNode = 0; jNode < pyyElemNodeCnt; jNode++)
                        {
                            int jNodeId = pyyNodes[jNode];
                            if (jNodeId == -1)
                            {
                                continue;
                            }
                            rmip[iDof, jNode + pElemNodeOffsets[5]] = pyyNu[iDof][jNode];
                        }
                    }
                    // for Navier-Stokes
                    {
                        rc = vx[0] + vy[1];
                        for (int jDof = 0; jDof < vDof; jDof++)
                        {
                            for (int jNode = 0; jNode < vElemNodeCnt; jNode++)
                            {
                                int jNodeId = vNodes[jNode];
                                if (jNodeId == -1)
                                {
                                    continue;
                                }
                                rcvj[jDof, jNode] = vNu[jDof][jNode];
                            }
                        }
                    }
                    /*NG
                    //-----------------
                    // for Pressure Poisson
                    {
                        rc = puv[0, 0] + puv[1, 1] + rho * (
                            v[0] * (vuv[0, 0][0] + vuv[0, 1][1]) +
                            v[1] * (vuv[1, 0][0] + vuv[1, 1][1])) -
                            rho * (gx[0] + gy[1]);
                        for (int jDof = 0; jDof < vDof; jDof++)
                        {
                            for (int jNode = 0; jNode < vElemNodeCnt; jNode++)
                            {
                                int jNodeId = vNodes[jNode];
                                if (jNodeId == -1)
                                {
                                    continue;
                                }
                                rcvj[jDof, jNode] = rho * (
                                    v[0] * vNuv[0, jDof][jNode] +
                                    v[1] * vNuv[1, jDof][jNode]);
                            }
                        }
                        for (int jNode = 0; jNode < pElemNodeCnt; jNode++)
                        {
                            int jNodeId = pNodes[jNode];
                            if (jNodeId == -1)
                            {
                                continue;
                            }
                            rcp[jNode + pElemNodeOffsets[0]] = pNuv[0, 0][jNode] + pNuv[1, 1][jNode];
                        }
                        for (int jNode = 0; jNode < pxElemNodeCnt; jNode++)
                        {
                            int jNodeId = pxNodes[jNode];
                            if (jNodeId == -1)
                            {
                                continue;
                            }
                            rcp[jNode + pElemNodeOffsets[1]] = pxNuv[0, 0][jNode] + pxNuv[1, 1][jNode];
                        }
                        for (int jNode = 0; jNode < pyElemNodeCnt; jNode++)
                        {
                            int jNodeId = pyNodes[jNode];
                            if (jNodeId == -1)
                            {
                                continue;
                            }
                            rcp[jNode + pElemNodeOffsets[2]] = pyNuv[0, 0][jNode] + pyNuv[1, 1][jNode];
                        }
                        for (int jNode = 0; jNode < pxxElemNodeCnt; jNode++)
                        {
                            int jNodeId = pxxNodes[jNode];
                            if (jNodeId == -1)
                            {
                                continue;
                            }
                            rcp[jNode + pElemNodeOffsets[3]] = pxxNuv[0, 0][jNode] + pxxNuv[1, 1][jNode];
                        }
                        for (int jNode = 0; jNode < pxyElemNodeCnt; jNode++)
                        {
                            int jNodeId = pxyNodes[jNode];
                            if (jNodeId == -1)
                            {
                                continue;
                            }
                            rcp[jNode + pElemNodeOffsets[4]] = pxyNuv[0, 0][jNode] + pxyNuv[1, 1][jNode];
                        }
                        for (int jNode = 0; jNode < pyyElemNodeCnt; jNode++)
                        {
                            int jNodeId = pyyNodes[jNode];
                            if (jNodeId == -1)
                            {
                                continue;
                            }
                            rcp[jNode + pElemNodeOffsets[5]] = pyyNuv[0, 0][jNode] + pyyNuv[1, 1][jNode];
                        }
                    }
                    */
                    //-----------------
                    // kvv
                    for (int row = 0; row < vElemNodeCnt; row++)
                    {
                        int rowNodeId = vNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < vElemNodeCnt; col++)
                        {
                            int colNodeId = vNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double[,] kvv1 = new double[vDof, vDof];
                            double[,] kvv2 = new double[vDof, vDof];
                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                for (int colDof = 0; colDof < vDof; colDof++)
                                {
                                    kvv1[rowDof, colDof] =
                                        detJWeight * taum * (v[0] * vNu[0][row] + v[1] * vNu[1][row]) *
                                        rmivj[rowDof, colDof, col];
                                    // for Navier-Stokes
                                    kvv2[rowDof, colDof] =
                                        detJWeight * tauc * rho * vNu[rowDof][row] * rcvj[colDof, col];
                                    /*NG
                                    // for Pressure Poisson
                                    kvv2[rowDof, colDof] =
                                        detJWeight * tauc * rho * vNu[rowDof][row] * rcvj[colDof, col];
                                    */
                                }
                            }

                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                for (int colDof = 0; colDof < vDof; colDof++)
                                {
                                    A[rowNodeId * vDof + rowDof, colNodeId * vDof + colDof] +=
                                        kvv1[rowDof, colDof] + kvv2[rowDof, colDof];
                                }
                            }
                        }
                    }
                    // kvp
                    for (int row = 0; row < vElemNodeCnt; row++)
                    {
                        int rowNodeId = vNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (uint colQuantityId = pQuantityId; colQuantityId < quantityCnt; colQuantityId++)
                        {
                            int coloffsetp = offsetps[colQuantityId - pQuantityId];
                            uint colpElemNodeCnt = pElemNodeCnts[colQuantityId - pQuantityId];
                            int[] colpNodes = pNodess[colQuantityId - pQuantityId];
                            double[] colpN = pNs[colQuantityId - pQuantityId];
                            double[] colpNx = pNxs[colQuantityId - pQuantityId];
                            double[] colpNy = pNys[colQuantityId - pQuantityId];

                            for (int col = 0; col < colpElemNodeCnt; col++)
                            {
                                int colNodeId = colpNodes[col];
                                if (colNodeId == -1)
                                {
                                    continue;
                                }

                                double[,] kv2p = new double[vDof, pDof];
                                double[,] kvp = new double[vDof, pDof];
                                uint workcolOffset = pElemNodeOffsets[colQuantityId - pQuantityId];
                                for (int rowDof = 0; rowDof < vDof; rowDof++)
                                {
                                    /*NG
                                    //----------------
                                    // for Pressure Poisson
                                    kv2p[rowDof, 0] =
                                        detJWeight * tauc * rho * vNu[rowDof][row] * rcp[col + workcolOffset];
                                    //----------------
                                    */
                                    kvp[rowDof, 0] =
                                        detJWeight * taum *
                                        (v[0] * vNu[0][row] + v[1] * vNu[1][row]) * rmip[rowDof, col + workcolOffset];
                                }

                                for (int rowDof = 0; rowDof < vDof; rowDof++)
                                {
                                    // for Navier-Stokes
                                    A[rowNodeId * vDof + rowDof, coloffsetp + colNodeId] += kvp[rowDof, 0];
                                    /*
                                    // for Pressure Poisson
                                    A[rowNodeId * vDof + rowDof, coloffsetp + colNodeId] +=
                                        kvp[rowDof, 0] + kv2p[rowDof, 0];
                                    */
                                }
                            }
                        }
                    }
                    ////////
                    // kpv
                    for (uint rowQuantityId = pQuantityId; rowQuantityId < quantityCnt; rowQuantityId++)
                    {
                        int rowoffsetp = offsetps[rowQuantityId - pQuantityId];
                        uint rowpElemNodeCnt = pElemNodeCnts[rowQuantityId - pQuantityId];
                        int[] rowpNodes = pNodess[rowQuantityId - pQuantityId];
                        double[] rowpN = pNs[rowQuantityId - pQuantityId];
                        double[] rowpNx = pNxs[rowQuantityId - pQuantityId];
                        double[] rowpNy = pNys[rowQuantityId - pQuantityId];
                        double[][] rowpNu = { rowpNx, rowpNy };
                        uint workrowOffset = pElemNodeOffsets[rowQuantityId - pQuantityId];

                        for (int row = 0; row < rowpElemNodeCnt; row++)
                        {
                            int rowNodeId = rowpNodes[row];
                            if (rowNodeId == -1)
                            {
                                continue;
                            }
                            for (int col = 0; col < vElemNodeCnt; col++)
                            {
                                int colNodeId = vNodes[col];
                                if (colNodeId == -1)
                                {
                                    continue;
                                }

                                double[,] kpv = new double[pDof, vDof];
                                for (int colDof = 0; colDof < vDof; colDof++)
                                {
                                    kpv[0, colDof] =
                                        detJWeight * (1.0 / rho) * taum *
                                        (rowpNu[0][row] * rmivj[0, colDof, col] + rowpNu[1][row] * rmivj[1, colDof, col]);
                                }

                                for (int colDof = 0; colDof < vDof; colDof++)
                                {
                                    A[rowoffsetp + rowNodeId, colNodeId * vDof + colDof] += kpv[0, colDof];
                                }
                            }
                        }
                    }
                    // kpp
                    for (uint rowQuantityId = pQuantityId; rowQuantityId < quantityCnt; rowQuantityId++)
                    {
                        int rowoffsetp = offsetps[rowQuantityId - pQuantityId];
                        uint rowpElemNodeCnt = pElemNodeCnts[rowQuantityId - pQuantityId];
                        int[] rowpNodes = pNodess[rowQuantityId - pQuantityId];
                        double[] rowpN = pNs[rowQuantityId - pQuantityId];
                        double[] rowpNx = pNxs[rowQuantityId - pQuantityId];
                        double[] rowpNy = pNys[rowQuantityId - pQuantityId];
                        double[][] rowpNu = { rowpNx, rowpNy };
                        uint workrowOffset = pElemNodeOffsets[rowQuantityId - pQuantityId];

                        for (int row = 0; row < rowpElemNodeCnt; row++)
                        {
                            int rowNodeId = rowpNodes[row];
                            if (rowNodeId == -1)
                            {
                                continue;
                            }

                            for (uint colQuantityId = pQuantityId; colQuantityId < quantityCnt; colQuantityId++)
                            {
                                int coloffsetp = offsetps[colQuantityId - pQuantityId];
                                uint colpElemNodeCnt = pElemNodeCnts[colQuantityId - pQuantityId];
                                int[] colpNodes = pNodess[colQuantityId - pQuantityId];
                                double[] colpN = pNs[colQuantityId - pQuantityId];
                                double[] colpNx = pNxs[colQuantityId - pQuantityId];
                                double[] colpNy = pNys[colQuantityId - pQuantityId];
                                double[][] colpNu = { colpNx, colpNy };
                                uint workcolOffset = pElemNodeOffsets[colQuantityId - pQuantityId];

                                for (int col = 0; col < colpElemNodeCnt; col++)
                                {
                                    int colNodeId = colpNodes[col];
                                    if (colNodeId == -1)
                                    {
                                        continue;
                                    }

                                    double[,] kpp = new double[pDof, pDof];
                                    kpp[0, 0] =
                                        detJWeight * (1.0 / rho) * taum *
                                        (rowpNu[0][row] * rmip[0, col + workcolOffset] +
                                        rowpNu[1][row] * rmip[1, col + workcolOffset]);

                                    A[rowoffsetp + rowNodeId, coloffsetp + colNodeId] += kpp[0, 0];
                                }
                            }
                        }
                    }
                    /////
                }

                for (int row = 0; row < vElemNodeCnt; row++)
                {
                    int rowNodeId = vNodes[row];
                    if (rowNodeId == -1)
                    {
                        continue;
                    }
                    for (int rowDof = 0; rowDof < vDof; rowDof++)
                    {
                        B[rowNodeId * vDof + rowDof] += rho * g[rowDof] * vSN[row];
                    }
                }
            }
        }
    }
}
