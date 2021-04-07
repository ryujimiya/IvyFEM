using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class FluidFIC2DTDFEM
    {
        private void InitFIC()
        {
            int nDim = 2;
            uint vQuantityId = 0;
            uint pQuantityId = 1;
            uint qQuantityId = 2;
            //uint cQuantityId = 3;
            int vDof = 2;
            int pDof = 1;
            int qDof = 2;
            int cDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int qNodeCnt = (int)World.GetNodeCount(qQuantityId);
            //int cNodeCnt = (int)World.GetNodeCount(cQuantityId);
            int pOffset = vNodeCnt * vDof;
            int qOffset = pOffset + pNodeCnt;
            //int cOffset = qOffset + qNodeCnt * qDof;

            UpdatedCoord = new double[vNodeCnt * nDim];
            for (int nodeId = 0; nodeId < vNodeCnt; nodeId++)
            {
                int coId = World.Node2Coord(vQuantityId, nodeId);
                double[] co = World.GetCoord(vQuantityId, coId);
                for (int iDof = 0; iDof < nDim; iDof++)
                {
                    UpdatedCoord[nodeId * nDim + iDof] = co[iDof];
                }
            }
        }

        private void UpdateFIC()
        {
            double dt = TimeStep;
            int nDim = 2;
            uint vQuantityId = 0;
            uint pQuantityId = 1;
            uint qQuantityId = 2;
            //uint cQuantityId = 3;
            int vDof = 2;
            int pDof = 1;
            int qDof = 2;
            //int cDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int qNodeCnt = (int)World.GetNodeCount(qQuantityId);
            //int cNodeCnt = (int)World.GetNodeCount(cQuantityId);
            int pOffset = vNodeCnt * vDof;
            int qOffset = pOffset + pNodeCnt;
            int cOffset = qOffset + qNodeCnt * qDof;

            for (int nodeId = 0; nodeId < vNodeCnt; nodeId++)
            {
                double[] v = new double[vDof];
                for (int iDof = 0; iDof < vDof; iDof++)
                {
                    v[iDof] = U[nodeId * vDof + iDof];
                }

                for (int iDof = 0; iDof < vDof; iDof++)
                {
                    UpdatedCoord[nodeId * nDim + iDof] += U[nodeId * vDof + iDof] * dt;
                }
            }
        }

        private void UpdateFEDisplacements()
        {
            int quantityCnt = World.GetQuantityCount();
            uint vQuantityId = 0;
            uint pQuantityId = 1;
            uint qQuantityId = 2;
            uint cQuantityId = 3;

            _UpdateFEDisplacements(vQuantityId);
            _UpdateFEDisplacements(pQuantityId);
            _UpdateFEDisplacements(qQuantityId);
            if (cQuantityId < quantityCnt)
            {
                _UpdateFEDisplacements(cQuantityId);
            }
        }

        private void _UpdateFEDisplacements(uint quantity)
        {
            double dt = TimeStep;
            uint vQuantityId = 0;
            uint pQuantityId = 1;
            uint qQuantityId = 2;
            //uint cQuantityId = 3;
            int vDof = 2;
            int pDof = 1;
            int qDof = 2;
            //int cDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int qNodeCnt = (int)World.GetNodeCount(qQuantityId);
            //int cNodeCnt = (int)World.GetNodeCount(cQuantityId);
            int pOffset = vNodeCnt * vDof;
            int qOffset = pOffset + pNodeCnt;
            //int cdOffset = qOffset + qNodeCnt * qDof;

            double[] displacementU = new double[vNodeCnt * vDof];

            for (int nodeId = 0; nodeId < vNodeCnt; nodeId++)
            {
                int coId = World.Node2Coord(vQuantityId, nodeId);
                double[] co1 = World.GetCoord(vQuantityId, coId);
                double[] co2 = new double[vDof];
                for (int iDof = 0; iDof < vDof; iDof++)
                {
                    co2[iDof] = UpdatedCoord[nodeId * vDof + iDof] + U[nodeId * vDof + iDof] * dt;
                }
                for (int iDof = 0; iDof < vDof; iDof++)
                {
                    displacementU[nodeId * vDof + iDof] = co2[iDof] - co1[iDof];
                }
            }
            World.UpdateFEDisplacements(quantity, displacementU);
        }

        private void ClearFEDisplacements()
        {
            int quantityCnt = World.GetQuantityCount();
            uint vQuantityId = 0;
            uint pQuantityId = 1;
            uint qQuantityId = 2;
            uint cQuantityId = 3;

            World.ClearFEDisplacements(vQuantityId);
            World.ClearFEDisplacements(pQuantityId);
            World.ClearFEDisplacements(qQuantityId);
            if (cQuantityId < quantityCnt)
            {
                World.ClearFEDisplacements(cQuantityId);
            }
        }

        private void CalcABFIC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            int nDim = 2;
            uint vQuantityId = 0;
            uint pQuantityId = 1;
            uint qQuantityId = 2;
            //uint cQuantityId = 3;
            int vDof = 2;
            int pDof = 1;
            int qDof = 2;
            //int cDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int qNodeCnt = (int)World.GetNodeCount(qQuantityId);
            //int cNodeCnt = (int)World.GetNodeCount(cQuantityId);
            int pOffset = vNodeCnt * vDof;
            int qOffset = pOffset + pNodeCnt;
            //int cOffset = qOffset + qNodeCnt * qDof;

            double dt = TimeStep;
            double beta = NewmarkBeta;
            double gamma = NewmarkGamma;
            var FV = World.GetFieldValue(ValueId);
            IList<uint> feIds = World.GetTriangleFEIds(vQuantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE vTriFE = World.GetTriangleFE(vQuantityId, feId);
                TriangleFE pTriFE = World.GetTriangleFE(pQuantityId, feId);
                TriangleFE qTriFE = World.GetTriangleFE(qQuantityId, feId);
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
                int[] qCoIds = qTriFE.NodeCoordIds;
                uint qElemNodeCnt = qTriFE.NodeCount;
                int[] qNodes = new int[qElemNodeCnt];
                for (int iNode = 0; iNode < qElemNodeCnt; iNode++)
                {
                    int coId = qCoIds[iNode];
                    int nodeId = World.Coord2Node(qQuantityId, coId);
                    qNodes[iNode] = nodeId;
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
                    double[] qN = qTriFE.CalcN(L);
                    double[][] qNu = qTriFE.CalcNu(L);
                    double[] qNx = qNu[0];
                    double[] qNy = qNu[1];

                    double detJ = vTriFE.GetDetJacobian(L);
                    double weight = ip.Weights[ipPt];
                    double detJWeight = (1.0 / 2.0) * weight * detJ;

                    double[] v = new double[vDof];
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
                        }
                    }

                    double Ae = vTriFE.GetArea();
                    var tau = new double[nDim];
                    for (int iDim = 0; iDim < nDim; iDim++)
                    {
                        double hi = Math.Sqrt(Ae);
                        //tau[iDim] = (3.0 * hi * hi) / (8.0 * mu);
                        double vi = Math.Sqrt(v[0] * v[0] + v[1] * v[1]);
                        double invTau1 = 2.0 * rho * vi / hi;
                        double invTau2 = 8.0 * mu / (3.0 * hi * hi);
                        tau[iDim] = 1.0 / (invTau1 + invTau2);
                    }

                    //---------------------------------------------------------------
                    // (1)
                    // M
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

                            double[,] m = new double[vDof, vDof];
                            m[0, 0] = detJWeight * rho * vN[row] * vN[col];
                            m[1, 1] = detJWeight * rho * vN[row] * vN[col];

                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                for (int colDof = 0; colDof < vDof; colDof++)
                                {
                                    A[rowNodeId * vDof + rowDof, colNodeId * vDof + colDof] +=
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
                    ////////////////////////////////////
                    // K
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

                            var bMati = new IvyFEM.Lapack.DoubleMatrix(3, 2);
                            bMati[0, 0] = vNx[row];
                            bMati[0, 1] = 0.0;
                            bMati[1, 0] = 0.0;
                            bMati[1, 1] = vNy[row];
                            bMati[2, 0] = vNy[row];
                            bMati[2, 1] = vNx[row];
                            bMati.Transpose();

                            var bMatj = new IvyFEM.Lapack.DoubleMatrix(3, 2);
                            bMatj[0, 0] = vNx[col];
                            bMatj[0, 1] = 0.0;
                            bMatj[1, 0] = 0.0;
                            bMatj[1, 1] = vNy[col];
                            bMatj[2, 0] = vNy[col];
                            bMatj[2, 1] = vNx[col];

                            var dMat = new IvyFEM.Lapack.DoubleMatrix(3, 3);
                            dMat[0, 0] = mu * 2.0;
                            dMat[1, 1] = mu * 2.0;
                            dMat[2, 2] = mu * 1.0;

                            var tmp = bMati * dMat;
                            var k = tmp * bMatj;

                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                for (int colDof = 0; colDof < vDof; colDof++)
                                {
                                    A[rowNodeId * vDof + rowDof, colNodeId * vDof + colDof] +=
                                        detJWeight * k[rowDof, colDof];
                               }
                            }
                        }
                    }
                    ////////////////////////////////////
                    // G
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

                            double[] mVec = { 1.0, 1.0, 0.0 }; // for 2D

                            var bMati = new IvyFEM.Lapack.DoubleMatrix(3, 2);
                            bMati[0, 0] = vNx[row];
                            bMati[0, 1] = 0.0;
                            bMati[1, 0] = 0.0;
                            bMati[1, 1] = vNy[row];
                            bMati[2, 0] = vNy[row];
                            bMati[2, 1] = vNx[row];
                            bMati.Transpose();

                            var tmp = new double[2];
                            for (int iDim = 0; iDim < 2; iDim++)
                            {
                                for (int jDim = 0; jDim < 3; jDim++)
                                {
                                    tmp[iDim] += bMati[iDim, jDim] * mVec[jDim];
                                }
                            }

                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                A[rowNodeId * vDof + rowDof, pOffset + colNodeId] +=
                                    -1.0 * detJWeight * tmp[rowDof] * pN[col];

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

                        double[] bodyF = new double[nDim];
                        for (int iDim = 0; iDim < nDim; iDim++)
                        {
                            bodyF[iDim] = detJWeight * vN[iDim] * rho * g[iDim];
                        }

                        for (int rowDof = 0; rowDof < vDof; rowDof++)
                        {
                            B[rowNodeId * vDof + rowDof] += bodyF[rowDof];
                        }
                    }
                    //---------------------------------------------------------------
                    // (2)
                    // G
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

                            double[] mVec = { 1.0, 1.0, 0.0 }; // for 2D

                            var bMati = new IvyFEM.Lapack.DoubleMatrix(3, 2);
                            bMati[0, 0] = vNx[col];
                            bMati[0, 1] = 0.0;
                            bMati[1, 0] = 0.0;
                            bMati[1, 1] = vNy[col];
                            bMati[2, 0] = vNy[col];
                            bMati[2, 1] = vNx[col];
                            bMati.Transpose();

                            var tmp = new double[2];
                            for (int iDim = 0; iDim < 2; iDim++)
                            {
                                for (int jDim = 0; jDim < 3; jDim++)
                                {
                                    tmp[iDim] += bMati[iDim, jDim] * mVec[jDim];
                                }
                            }

                            for (int colDof = 0; colDof < vDof; colDof++)
                            {
                                A[pOffset + rowNodeId, colNodeId * vDof + colDof] +=
                                    detJWeight * tmp[colDof] * pN[row];

                            }
                        }
                    }
                    // L
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

                            for (int k = 0; k < nDim; k++)
                            {
                                A[pOffset + rowNodeId, pOffset + colNodeId] +=
                                    detJWeight * tau[k] * pNu[k][row] * pNu[k][col];
                            }
                        }
                    }
                    // Q
                    for (int row = 0; row < pElemNodeCnt; row++)
                    {
                        int rowNodeId = pNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < qElemNodeCnt; col++)
                        {
                            int colNodeId = qNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double[] qMatij = new double[2];
                            for (int k = 0; k < qDof; k++)
                            {
                                qMatij[k] = detJWeight * tau[k] * pNu[k][row] * qN[col];
                            }

                            for (int colDof = 0; colDof < qDof; colDof++)
                            {
                                A[pOffset + rowNodeId, qOffset + colNodeId * qDof + colDof] += qMatij[colDof];
                            }
                        }
                    }
                    //---------------------------------------------------------------
                    // (3)
                    // Q
                    for (int row = 0; row < qElemNodeCnt; row++)
                    {
                        int rowNodeId = qNodes[row];
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

                            double[] qMatij = new double[2];
                            for (int k = 0; k < qDof; k++)
                            {
                                qMatij[k] = detJWeight * tau[k] * pNu[k][col] * qN[row];
                            }

                            for (int rowDof = 0; rowDof < qDof; rowDof++)
                            {
                                A[qOffset + rowNodeId * qDof + rowDof, pOffset + colNodeId] += qMatij[rowDof];
                            }
                        }
                    }
                    // hatM
                    for (int row = 0; row < qElemNodeCnt; row++)
                    {
                        int rowNodeId = qNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < qElemNodeCnt; col++)
                        {
                            int colNodeId = qNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double[] hatmij = new double[2];
                            for (int k = 0; k < 2; k++)
                            {
                                hatmij[k] = detJWeight * tau[k] * qN[row] * qN[col];
                            }

                            for (int rowDof = 0; rowDof < qDof; rowDof++)
                            {
                                int colDof = rowDof;
                                A[qOffset + rowNodeId * qDof + rowDof, qOffset + colNodeId * qDof + colDof] +=
                                    hatmij[rowDof];
                            }
                        }
                    }
                }
            }
        }
    }
}
