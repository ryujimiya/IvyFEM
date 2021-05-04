using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class FluidRANS2DTDFEM
    {
        private void CalcABRNGKEpsilon(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint vQuantityId = 0;
            uint pQuantityId = 1;
            uint kpQuantityId = 2;
            uint epQuantityId = 3;
            int vDof = 2;
            int pDof = 1;
            int kpDof = 1;
            int epDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int kpNodeCnt = (int)World.GetNodeCount(kpQuantityId);
            int epNodeCnt = (int)World.GetNodeCount(epQuantityId);
            int pOffset = vNodeCnt * vDof;
            int kpOffset = pOffset + pNodeCnt;
            int epOffset = kpOffset + kpNodeCnt;

            double dt = TimeStep;
            double beta = NewmarkBeta;
            double gamma = NewmarkGamma;
            var FV = World.GetFieldValue(ValueId);
            var kpFV = World.GetFieldValue(KpValueId);
            var epFV = World.GetFieldValue(EpValueId);
            IList<uint> feIds = World.GetTriangleFEIds(vQuantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE vTriFE = World.GetTriangleFE(vQuantityId, feId);
                TriangleFE pTriFE = World.GetTriangleFE(pQuantityId, feId);
                TriangleFE kpTriFE = World.GetTriangleFE(kpQuantityId, feId);
                TriangleFE epTriFE = World.GetTriangleFE(epQuantityId, feId);
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
                int[] kpCoIds = kpTriFE.NodeCoordIds;
                uint kpElemNodeCnt = kpTriFE.NodeCount;
                int[] kpNodes = new int[kpElemNodeCnt];
                for (int iNode = 0; iNode < kpElemNodeCnt; iNode++)
                {
                    int coId = kpCoIds[iNode];
                    int nodeId = World.Coord2Node(kpQuantityId, coId);
                    kpNodes[iNode] = nodeId;
                }
                int[] epCoIds = epTriFE.NodeCoordIds;
                uint epElemNodeCnt = epTriFE.NodeCount;
                int[] epNodes = new int[epElemNodeCnt];
                for (int iNode = 0; iNode < epElemNodeCnt; iNode++)
                {
                    int coId = epCoIds[iNode];
                    int nodeId = World.Coord2Node(epQuantityId, coId);
                    epNodes[iNode] = nodeId;
                }

                Material ma0 = World.GetMaterial(vTriFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is NewtonFluidMaterial);
                var ma = ma0 as NewtonFluidMaterial;
                double rho = ma.MassDensity;
                double mu = ma.Mu;
                double[] g = { ma.GravityX, ma.GravityY };

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
                    double[] kpN = kpTriFE.CalcN(L);
                    double[][] kpNu = kpTriFE.CalcNu(L);
                    double[] kpNx = kpNu[0];
                    double[] kpNy = kpNu[1];
                    double[] epN = epTriFE.CalcN(L);
                    double[][] epNu = epTriFE.CalcNu(L);
                    double[] epNx = epNu[0];
                    double[] epNy = epNu[1];

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
                            double vValue = U[nodeId * vDof + iDof];
                            v[iDof] += vValue * vN[iNode];
                            vx[iDof] += vValue * vNx[iNode];
                            vy[iDof] += vValue * vNy[iNode];
                        }
                    }
                    double[][] vu = { vx, vy };

                    double press = 0.0;
                    double pressX = 0.0;
                    double pressY = 0.0;
                    for (int iNode = 0; iNode < pElemNodeCnt; iNode++)
                    {
                        int nodeId = pNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }

                        double pValue = U[pOffset + nodeId];
                        press += pValue * pN[iNode];
                        pressX += pValue * pNx[iNode];
                        pressY += pValue * pNy[iNode];
                    }

                    double kappa = 0.0;
                    double kappaX = 0.0;
                    double kappaY = 0.0;
                    for (int iNode = 0; iNode < kpElemNodeCnt; iNode++)
                    {
                        int nodeId = kpNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }

                        double kpValue = U[kpOffset + nodeId];
                        kappa += kpValue * kpN[iNode];
                        kappaX += kpValue * kpNx[iNode];
                        kappaY += kpValue * kpNy[iNode];
                    }
                    double[] kappaU = { kappaX, kappaY };

                    double epsilon = 0.0;
                    double epsilonX = 0.0;
                    double epsilonY = 0.0;
                    for (int iNode = 0; iNode < epElemNodeCnt; iNode++)
                    {
                        int nodeId = epNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }

                        double epValue = U[epOffset + nodeId];
                        epsilon += epValue * epN[iNode];
                        epsilonX += epValue * epNx[iNode];
                        epsilonY += epValue * epNy[iNode];
                    }
                    double[] epsilonU = { epsilonX, epsilonY };
                    //---------------------------
                    if (Math.Abs(kappa) < 1.0e-20)
                    {
                        kappa = 1.0e-20;
                    }
                    if (Math.Abs(epsilon) < 1.0e-20)
                    {
                        epsilon = 1.0e-20;
                    }
                    //---------------------------

                    double c1ep = 1.42;
                    double c2ep = 1.68;
                    double alphaK = 1.39;
                    double alphaE = 1.39;
                    double cmu = 0.0485;
                    double mut = rho * cmu * (kappa * kappa / epsilon);
                    double[,] sij = new double[vDof, vDof];
                    for (int iDof = 0; iDof < vDof; iDof++)
                    {
                        for (int jDof = 0; jDof < vDof; jDof++)
                        {
                            sij[iDof, jDof] = (1.0 / 2.0) * (vu[jDof][iDof] + vu[iDof][jDof]);
                        }
                    }
                    double gk = 0.0;
                    for (int iDof = 0; iDof < vDof; iDof++)
                    {
                        for (int jDof = 0; jDof < vDof; jDof++)
                        {
                            gk += 2.0 * mut * sij[iDof, jDof] * sij[iDof, jDof];
                        }
                    }
                    double eta0 = 4.377;
                    double betaRNG = 0.012;
                    double eta = 0.0;
                    for (int iDof = 0; iDof < vDof; iDof++)
                    {
                        for (int jDof = 0; jDof < vDof; jDof++)
                        {
                            eta += 2.0 * sij[iDof, jDof] * sij[iDof, jDof];
                        }
                    }
                    eta = (kappa / epsilon) * Math.Sqrt(eta);

                    //---------------------------------------------------------------
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
                            for (int iDof = 0; iDof < vDof; iDof++)
                            {
                                {
                                    int jDof = iDof;
                                    for (int kDof = 0; kDof < vDof; kDof++)
                                    {
                                        kvv1[iDof, jDof] +=
                                            detJWeight * rho * vN[row] * v[kDof] * vNu[kDof][col];
                                    }
                                }
                            }

                            double[,] kvv3 = new double[vDof, vDof];
                            for (int iDof = 0; iDof < vDof; iDof++)
                            {
                                for (int jDof = 0; jDof < vDof; jDof++)
                                {
                                    kvv3[iDof, jDof] = 
                                        detJWeight * (mu + mut) * vNu[jDof][row] * vNu[iDof][col];
                                }
                            }

                            double[,] kvv4 = new double[vDof, vDof];
                            for (int iDof = 0; iDof < vDof; iDof++)
                            {
                                {
                                    int jDof = iDof;
                                    for (int kDof = 0; kDof < vDof; kDof++)
                                    {
                                        kvv4[iDof, jDof] += 
                                            detJWeight * (mu + mut) * vNu[kDof][row] * vNu[kDof][col];
                                    }
                                }
                            }

                            double[,] m = new double[vDof, vDof];
                            m[0, 0] = detJWeight * rho * vN[row] * vN[col];
                            m[1, 1] = detJWeight * rho * vN[row] * vN[col];

                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                for (int colDof = 0; colDof < vDof; colDof++)
                                {
                                    A[rowNodeId * vDof + rowDof, colNodeId * vDof + colDof] +=
                                        kvv1[rowDof, colDof] + kvv3[rowDof, colDof] + kvv4[rowDof, colDof] +
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

                            double[,] kvp2 = new double[vDof, pDof];
                            for (int iDof = 0; iDof < vDof; iDof++)
                            {
                                kvp2[iDof, 0] = -1.0 * detJWeight * vNu[iDof][row] * pN[col];
                            }

                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                A[rowNodeId * vDof + rowDof, pOffset + colNodeId] +=
                                    kvp2[rowDof, 0];

                            }
                        }
                    }
                    // f
                    for (int row = 0; row < vElemNodeCnt; row++)
                    {
                        int rowNodeId = vNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }

                        double[] f = new double[vDof];
                        for (int iDof = 0; iDof < vDof; iDof++)
                        {
                            f[iDof] = detJWeight * (2.0 / 3.0) * rho * vN[iDof] * kappaU[iDof];
                        }

                        for (int rowDof = 0; rowDof < vDof; rowDof++)
                        {
                            B[rowNodeId * vDof + rowDof] += - f[rowDof];
                        }
                    }

                    //---------------------------------------------------------------
                    for (int row = 0; row < pElemNodeCnt; row++)
                    {
                        int rowNodeId = pNodes[row];
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
                            for (int jDof = 0; jDof < vDof; jDof++)
                            {
                                kpv[0, jDof] = detJWeight * rho * pN[row] * vNu[jDof][col];
                            }

                            for (int colDof = 0; colDof < vDof; colDof++)
                            {
                                A[pOffset + rowNodeId, colNodeId * vDof + colDof] +=
                                    kpv[0, colDof];

                            }
                        }
                    }

                    ///////////////////////////////////////////////////////////
                    //---------------------------------------------------------------
                    for (int row = 0; row < kpElemNodeCnt; row++)
                    {
                        int rowNodeId = kpNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < kpElemNodeCnt; col++)
                        {
                            int colNodeId = kpNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }
                            int colCoId = kpCoIds[col];
                            double[] u = kpFV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                            double[] vel = kpFV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                            double[] acc = kpFV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);

                            double kkk1 = 0.0;
                            for (int kDof = 0; kDof < vDof; kDof++)
                            {
                                kkk1 += detJWeight * rho * kpN[row] * v[kDof] * kpNu[kDof][col];
                            }
                            double kkk2 = 0.0;
                            for (int kDof = 0; kDof < vDof; kDof++)
                            {
                                kkk2 +=
                                    detJWeight * (alphaK * (mu + mut)) * kpNu[kDof][row] * kpNu[kDof][col];
                            }

                            double kkk3 = detJWeight * rho * (epsilon / kappa) * kpN[row] * kpN[col];

                            double m = detJWeight * rho * kpN[row] * kpN[col];

                            A[kpOffset + rowNodeId, kpOffset + colNodeId] +=
                                kkk1 + kkk2 + kkk3 +
                                (gamma / (beta * dt)) * m;

                            B[kpOffset + rowNodeId] +=
                                m * (
                                    (gamma / (beta * dt)) * u[0] -
                                    (1.0 - gamma / beta) * vel[0] -
                                    dt * (1.0 - gamma / (2.0 * beta)) * acc[0]
                                );
                        }
                    }

                    // f
                    for (int row = 0; row < kpElemNodeCnt; row++)
                    {
                        int rowNodeId = kpNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }

                        double f = -1.0 * detJWeight * kpN[row] * gk;
                        B[kpOffset + rowNodeId] += - f;
                    }

                    //---------------------------------------------------------------
                    for (int row = 0; row < epElemNodeCnt; row++)
                    {
                        int rowNodeId = epNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < epElemNodeCnt; col++)
                        {
                            int colNodeId = epNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }
                            int colCoId = epCoIds[col];
                            double[] u = epFV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                            double[] vel = epFV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                            double[] acc = epFV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);

                            double kee1 = 0.0;
                            for (int kDof = 0; kDof < vDof; kDof++)
                            {
                                kee1 += detJWeight * rho * epN[row] * v[kDof] * epNu[kDof][col];
                            }
                            double kee2 = 0.0;
                            for (int kDof = 0; kDof < vDof; kDof++)
                            {
                                kee2 +=
                                    detJWeight * (alphaE * (mu + mut)) * epNu[kDof][row] * epNu[kDof][col];
                            }

                            double kee3 =
                                -1.0 * detJWeight * epN[row] * (epsilon / kappa) *
                                (- c2ep * rho * epN[col]);

                            double m = detJWeight * rho * epN[row] * epN[col];

                            A[epOffset + rowNodeId, epOffset + colNodeId] +=
                                kee1 + kee2 + kee3 +
                                (gamma / (beta * dt)) * m;

                            B[epOffset + rowNodeId] +=
                                m * (
                                    (gamma / (beta * dt)) * u[0] -
                                    (1.0 - gamma / beta) * vel[0] -
                                    dt * (1.0 - gamma / (2.0 * beta)) * acc[0]
                                );
                        }
                    }
                    // f
                    for (int row = 0; row < epElemNodeCnt; row++)
                    {
                        int rowNodeId = epNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        double f1 =
                            -1.0 * detJWeight * epN[row] * (epsilon / kappa) *
                            (c1ep * gk);
                        double f2 =
                            detJWeight * epN[row] *
                            (eta * (1.0 - eta / eta0) / (1.0 + betaRNG * eta * eta * eta)) *
                            (epsilon / kappa) * gk;
                        B[epOffset + rowNodeId] +=
                            - f1 - f2;
                    }
                }
            }
        }
    }
}
