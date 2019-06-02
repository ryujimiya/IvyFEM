using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Fluid2DFEM
    {
        private void CalcStdGVorticityAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint wQuantityId = 0;
            uint pQuantityId = 1;
            int wNodeCnt = (int)World.GetNodeCount(wQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int offset = wNodeCnt;

            IList<uint> feIds = World.GetTriangleFEIds(wQuantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE wTriFE = World.GetTriangleFE(wQuantityId, feId);
                TriangleFE pTriFE = World.GetTriangleFE(pQuantityId, feId);
                uint vertexCnt = wTriFE.VertexCount;
                for (int iVertex = 0; iVertex < vertexCnt; iVertex++)
                {
                    System.Diagnostics.Debug.Assert(wTriFE.VertexCoordIds[iVertex] == pTriFE.VertexCoordIds[iVertex]);
                }

                int[] wCoIds = wTriFE.NodeCoordIds;
                uint wElemNodeCnt = wTriFE.NodeCount;
                int[] wNodes = new int[wElemNodeCnt];
                for (int iNode = 0; iNode < wElemNodeCnt; iNode++)
                {
                    int coId = wCoIds[iNode];
                    int nodeId = World.Coord2Node(wQuantityId, coId);
                    wNodes[iNode] = nodeId;
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

                Material ma0 = World.GetMaterial(wTriFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is NewtonFluidMaterial);
                var ma = ma0 as NewtonFluidMaterial;
                double rho = ma.MassDensity;
                double mu = ma.Mu;
                double[] g = { ma.GravityX, ma.GravityY };

                IntegrationPoints ip = TriangleFE.GetIntegrationPoints(World.TriIntegrationPointCount);//Point7
                for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                {
                    double[] L = ip.Ls[ipPt];
                    double[] wN = wTriFE.CalcN(L);
                    double[][] wNu = wTriFE.CalcNu(L);
                    double[] wNx = wNu[0];
                    double[] wNy = wNu[1];
                    double[] pN = pTriFE.CalcN(L);
                    double[][] pNu = pTriFE.CalcNu(L);
                    double[] pNx = pNu[0];
                    double[] pNy = pNu[1];

                    double detJ = wTriFE.GetDetJacobian(L);
                    double weight = ip.Weights[ipPt];
                    double detJWeight = (1.0 / 2.0) * weight * detJ;

                    double w = 0;
                    double wx = 0;
                    double wy = 0;
                    for (int iNode = 0; iNode < wElemNodeCnt; iNode++)
                    {
                        int nodeId = wNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }
                        double wValue = U[nodeId];
                        w += wValue * wN[iNode];
                        wx += wValue * wNx[iNode];
                        wy += wValue * wNy[iNode];
                    }
                    double p = 0;
                    double px = 0;
                    double py = 0;
                    for (int iNode = 0; iNode < pElemNodeCnt; iNode++)
                    {
                        int nodeId = pNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }
                        double pValue = U[offset + nodeId];
                        p += pValue * pN[iNode];
                        px += pValue * pNx[iNode];
                        py += pValue * pNy[iNode];
                    }
                    double[] v = new double[2];
                    v[0] = py;
                    v[1] = -px;

                    for (int row = 0; row < wElemNodeCnt; row++)
                    {
                        int rowNodeId = wNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < wElemNodeCnt; col++)
                        {
                            int colNodeId = wNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double kww1 = detJWeight * mu * (wNx[row] * wNx[col] + wNy[row] * wNy[col]);
                            double kww2 = detJWeight * rho * wN[row] * (v[0] * wNx[col] + v[1] * wNy[col]);
                            A[rowNodeId, colNodeId] +=
                                kww1 + kww2;
                        }
                    }

                    for (int row = 0; row < pElemNodeCnt; row++)
                    {
                        int rowNodeId = pNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < wElemNodeCnt; col++)
                        {
                            int colNodeId = wNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double kpw = -detJWeight * pN[row] * wN[col];
                            A[offset + rowNodeId, colNodeId] += kpw;
                        }
                    }

                    for (int row = 0; row < pElemNodeCnt; row++)
                    {
                        int rowNodeId = pNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < pElemNodeCnt; col++)
                        {
                            int colNodeId = pNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double kpp = detJWeight * (pNx[row] * pNx[col] + pNy[row] * pNy[col]);
                            A[offset + rowNodeId, offset + colNodeId] += kpp;
                        }
                    }

                    for (int row = 0; row < wElemNodeCnt; row++)
                    {
                        int rowNodeId = wNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        double f = 0;
                        for (int kNode = 0; kNode < wElemNodeCnt; kNode++)
                        {
                            int kNodeId = wNodes[kNode];
                            if (kNodeId == -1)
                            {
                                continue;
                            }
                            f += detJWeight * rho * (wNx[kNode] * g[1] - wNy[kNode] * g[0]);
                        }
                        B[rowNodeId] += f;
                    }
                }
            }
        }
    }
}
