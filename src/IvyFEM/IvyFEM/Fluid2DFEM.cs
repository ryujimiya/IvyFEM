using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Fluid2DFEM : FEM
    {
        // output
        public double[] U { get; private set; } = null;

        public Fluid2DFEM(FEWorld world)
        {
            World = world;
        }

        private void CalcAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint vQuantityId = 0;
            uint pQuantityId = 1;
            int vDof = 2;
            int pDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int offset = vNodeCnt * vDof;

            IList<uint> feIds = World.GetTriangleFEIds(vQuantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE vTriFE = World.GetTriangleFE(vQuantityId, feId);
                TriangleFE pTriFE = World.GetTriangleFE(pQuantityId, feId);
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

                Material ma0 = World.GetMaterial(vTriFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is NewtonFluidMaterial);
                var ma = ma0 as NewtonFluidMaterial;
                double rho = ma.MassDensity;
                double mu = ma.Mu;
                double[] g = { ma.GravityX, ma.GravityY };

                double[] vSN = vTriFE.CalcSN();
                IntegrationPoints ip = TriangleFE.GetIntegrationPoints(World.TriIntegrationPointCount);//Point7
                for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                {
                    double[] L = ip.Ls[ipPt];
                    double[] vN = vTriFE.CalcN(L);
                    double[][] vNu = vTriFE.CalcNu(L);
                    double[] vNx = vNu[0];
                    double[] vNy = vNu[1];
                    double[] pN = pTriFE.CalcN(L);
                    double[][] pNu = pTriFE.CalcNu(L);
                    double[] pNx = pNu[0];
                    double[] pNy = pNu[1];

                    double detJ = vTriFE.GetDetJacobian(L);
                    double weight = ip.Weights[ipPt];
                    double detJWeight = (1.0 / 2.0) * weight * detJ;

                    double[] v = new double[vDof];
                    double[] vx = new double[vDof];
                    double[] vy = new double[vDof];
                    for (int iNode = 0; iNode < vElemNodeCnt; iNode++)
                    {
                        int nodeId = vNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }
                        for (int iDof = 0; iDof < vDof; iDof++)
                        {
                            v[iDof] += U[nodeId * vDof + iDof] * vN[iNode];
                            vx[iDof] += U[nodeId * vDof + iDof] * vNx[iNode];
                            vy[iDof] += U[nodeId * vDof + iDof] * vNy[iNode];
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
                            kvv1[0, 0] = detJWeight * mu * (vNx[row] * vNx[col] + vNy[row] * vNy[col]);
                            kvv1[1, 1] = detJWeight * mu * (vNx[row] * vNx[col] + vNy[row] * vNy[col]);

                            double[,] kvv2 = new double[vDof, vDof];
                            kvv2[0, 0] = detJWeight * rho * vN[row] * (
                                vN[col] * vx[0] + v[0] * vNx[col] + v[1] * vNy[col]);
                            kvv2[0, 1] = detJWeight * rho * vN[row] * vN[col] * vy[0];
                            kvv2[1, 0] = detJWeight * rho * vN[row] * vN[col] * vx[1];
                            kvv2[1, 1] = detJWeight * rho * vN[row] * (
                                vN[col] * vy[1] + v[0] * vNx[col] + v[1] * vNy[col]);

                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                for (int colDof = 0; colDof < vDof; colDof++)
                                {
                                    A[rowNodeId * vDof + rowDof, colNodeId * vDof + colDof] +=
                                        kvv1[rowDof, colDof] + kvv2[rowDof, colDof];

                                    B[rowNodeId * vDof + rowDof] +=
                                        kvv2[rowDof, colDof] * U[colNodeId * vDof + colDof];
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
                        for (int col = 0; col < pElemNodeCnt; col++)
                        {
                            int colNodeId = pNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double[,] kvp = new double[vDof, pDof];
                            kvp[0, 0] = -detJWeight * vNx[row] * pN[col];
                            kvp[1, 0] = -detJWeight * vNy[row] * pN[col];

                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                A[rowNodeId * vDof + rowDof, offset + colNodeId] += kvp[rowDof, 0];
                                A[offset + colNodeId, rowNodeId * vDof + rowDof] += kvp[rowDof, 0];
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
                    for (int rowDof = 0; rowDof < vDof; rowDof++)
                    {
                        B[rowNodeId * vDof + rowDof] += rho * g[rowDof] * vSN[row];
                    }
                }
            }
        }

        public override void Solve()
        {
            // Newton Raphson
            double sqNorm = 0;
            double sqInvNorm0 = 0;
            double convRatio = 1.0e+2 * IvyFEM.Linear.Constants.ConvRatioTolerance;
            double tolerance = convRatio;
            const int maxIter = IvyFEM.Linear.Constants.MaxIter;
            int iter = 0;

            uint vQuantityId = 0;
            int vDof = 2;
            uint pQuantityId = 1;
            int pDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int nodeCnt = vNodeCnt * vDof + pNodeCnt * pDof;
            U = new double[nodeCnt];

            for (iter = 0; iter < maxIter; iter++)
            {
                var A = new IvyFEM.Linear.DoubleSparseMatrix(nodeCnt, nodeCnt);
                var B = new double[nodeCnt];

                CalcAB(A, B);

                DoubleSetFixedCadsCondtion(A, B, new int[] { vNodeCnt, pNodeCnt }, new int[] { vDof, pDof });

                double[] AU = A * U;
                double[] R = IvyFEM.Lapack.Functions.daxpy(-1.0, B, AU);
                sqNorm = IvyFEM.Lapack.Functions.ddot(R, R);
                if (iter == 0)
                {
                    if (sqNorm < IvyFEM.Constants.PrecisionLowerLimit)
                    {
                        convRatio = 0;
                        break;
                    }
                    sqInvNorm0 = 1.0 / sqNorm;
                }
                else
                {
                    convRatio = Math.Sqrt(sqNorm * sqInvNorm0);
                    if (sqNorm * sqInvNorm0 < tolerance * tolerance)
                    {
                        break;
                    }
                }

                //---------------------------------------------------
                double[] X;
                Solver.DoubleSolve(out X, A, B);
                U = X;
                //---------------------------------------------------
            }
            System.Diagnostics.Debug.WriteLine("Newton Raphson iter = " + iter + " norm = " + convRatio);
            System.Diagnostics.Debug.Assert(iter < maxIter);
        }
    }
}
