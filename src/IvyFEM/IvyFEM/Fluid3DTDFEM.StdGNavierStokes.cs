using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Fluid3DTDFEM
    {
        private void CalcStdGNavierStokesAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint vQuantityId = 0;
            uint pQuantityId = 1;
            int vDof = 3;
            int pDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int offset = vNodeCnt * vDof;

            double dt = TimeStep;
            double beta = NewmarkBeta;
            double gamma = NewmarkGamma;
            var FV = World.GetFieldValue(ValueId);
            IList<uint> feIds = World.GetTetrahedronFEIds(vQuantityId);
            foreach (uint feId in feIds)
            {
                TetrahedronFE vTetFE = World.GetTetrahedronFE(vQuantityId, feId);
                TetrahedronFE pTetFE = World.GetTetrahedronFE(pQuantityId, feId);
                uint vertexCnt = vTetFE.VertexCount;
                for (int iVertex = 0; iVertex < vertexCnt; iVertex++)
                {
                    System.Diagnostics.Debug.Assert(vTetFE.VertexCoordIds[iVertex] == pTetFE.VertexCoordIds[iVertex]);
                }

                int[] vCoIds = vTetFE.NodeCoordIds;
                uint vElemNodeCnt = vTetFE.NodeCount;
                int[] vNodes = new int[vElemNodeCnt];
                for (int iNode = 0; iNode < vElemNodeCnt; iNode++)
                {
                    int coId = vCoIds[iNode];
                    int nodeId = World.Coord2Node(vQuantityId, coId);
                    vNodes[iNode] = nodeId;
                }
                int[] pCoIds = pTetFE.NodeCoordIds;
                uint pElemNodeCnt = pTetFE.NodeCount;
                int[] pNodes = new int[pElemNodeCnt];
                for (int iNode = 0; iNode < pElemNodeCnt; iNode++)
                {
                    int coId = pCoIds[iNode];
                    int nodeId = World.Coord2Node(pQuantityId, coId);
                    pNodes[iNode] = nodeId;
                }

                Material ma0 = World.GetMaterial(vTetFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is NewtonFluidMaterial);
                var ma = ma0 as NewtonFluidMaterial;
                double rho = ma.MassDensity;
                double mu = ma.Mu;
                double[] g = { ma.GravityX, ma.GravityY, ma.GravityZ };

                double[] vSN = vTetFE.CalcSN();
                IntegrationPoints ip = TetrahedronFE.GetIntegrationPoints(World.TetIntegrationPointCount);
                for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                {
                    double[] L = ip.Ls[ipPt];
                    double[] vN = vTetFE.CalcN(L);
                    double[][] vNu = vTetFE.CalcNu(L);
                    double[] vNx = vNu[0];
                    double[] vNy = vNu[1];
                    double[] vNz = vNu[2];
                    double[] pN = pTetFE.CalcN(L);
                    double[][] pNu = pTetFE.CalcNu(L);
                    double[] pNx = pNu[0];
                    double[] pNy = pNu[1];
                    double[] pNz = pNu[2];

                    double detJ = vTetFE.GetDetJacobian(L);
                    double weight = ip.Weights[ipPt];
                    double detJWeight = (1.0 / 6.0) * weight * detJ;

                    double[] v = new double[vDof];
                    double[] vx = new double[vDof];
                    double[] vy = new double[vDof];
                    double[] vz = new double[vDof];
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
                            vz[iDof] += vValue * vNz[iNode];
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
                            int colCoId = vCoIds[col];
                            double[] u = FV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                            double[] vel = FV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                            double[] acc = FV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);

                            double[,] kvv1 = new double[vDof, vDof];
                            kvv1[0, 0] = detJWeight * mu * (vNx[row] * vNx[col] +
                                vNx[row] * vNx[col] + vNy[row] * vNy[col] + vNz[row] * vNz[col]);
                            kvv1[0, 1] = detJWeight * mu * vNy[row] * vNx[col];
                            kvv1[0, 2] = detJWeight * mu * vNz[row] * vNx[col];
                            kvv1[1, 0] = detJWeight * mu * vNx[row] * vNy[col];
                            kvv1[1, 1] = detJWeight * mu * (vNy[row] * vNy[col] +
                                vNx[row] * vNx[col] + vNy[row] * vNy[col] + vNz[row] * vNz[col]);
                            kvv1[1, 2] = detJWeight * mu * vNz[row] * vNy[col];
                            kvv1[2, 0] = detJWeight * mu * vNx[row] * vNz[col];
                            kvv1[2, 1] = detJWeight * mu * vNy[row] * vNz[col];
                            kvv1[2, 2] = detJWeight * mu * (vNz[row] * vNz[col] +
                                vNx[row] * vNx[col] + vNy[row] * vNy[col] + vNz[row] * vNz[col]);

                            double[,] kvv2 = new double[vDof, vDof];
                            kvv2[0, 0] = detJWeight * rho * vN[row] * (
                                v[0] * vNx[col] + v[1] * vNy[col] + v[2] * vNz[col]);
                            kvv2[0, 1] = 0;
                            kvv2[0, 2] = 0;
                            kvv2[1, 0] = 0;
                            kvv2[1, 1] = detJWeight * rho * vN[row] * (
                                v[0] * vNx[col] + v[1] * vNy[col] + v[2] * vNz[col]);
                            kvv2[1, 2] = 0;
                            kvv2[2, 0] = 0;
                            kvv2[2, 1] = 0;
                            kvv2[2, 2] = detJWeight * rho * vN[row] * (
                                v[0] * vNx[col] + v[1] * vNy[col] + v[2] * vNz[col]);

                            double[,] m = new double[vDof, vDof];
                            m[0, 0] = detJWeight * rho * vN[row] * vN[col];
                            m[1, 1] = detJWeight * rho * vN[row] * vN[col];
                            m[2, 2] = detJWeight * rho * vN[row] * vN[col];

                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                for (int colDof = 0; colDof < vDof; colDof++)
                                {
                                    A[rowNodeId * vDof + rowDof, colNodeId * vDof + colDof] +=
                                        kvv1[rowDof, colDof] + kvv2[rowDof, colDof] +
                                        (gamma / (beta * dt)) * m[rowDof, colDof];

                                    B[rowNodeId * vDof + rowDof] +=
                                        m[rowDof, colDof] * (
                                            (gamma / (beta * dt)) * u[colDof] -
                                            (1.0 - gamma / beta) * vel[colDof] -
                                            dt * (1.0 - gamma / (2.0 * beta)) * acc[colDof]
                                        );
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
                            kvp[2, 0] = -detJWeight * vNz[row] * pN[col];

                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                A[rowNodeId * vDof + rowDof, offset + colNodeId] += kvp[rowDof, 0];
                                A[offset + colNodeId, rowNodeId * vDof + rowDof] += -kvp[rowDof, 0];
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
    }
}
