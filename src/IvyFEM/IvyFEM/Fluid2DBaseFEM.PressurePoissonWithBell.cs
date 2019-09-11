using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class Fluid2DBaseFEM
    {
        protected void SetPressurePoissonWithBellSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B,
            double dt, double beta, double gamma, uint valueId = 0)
        {
            // for fixed velocity
            SetPressurePoissonWithBellNoConstraintNeumannBCOfPressure(A, B, dt, beta, gamma, valueId);
            SetPressurePoissonWithBellNormalInflowNeumannBCOfPressure(A, B, dt, beta, gamma, valueId);
            // for outflow
            SetPressurePoissonWithBellOutflowNeumannBCOfPressure(A, B, dt, beta, gamma, valueId);
        }

        private void SetPressurePoissonWithBellNoConstraintNeumannBCOfPressure(
            IvyFEM.Linear.DoubleSparseMatrix A, double[] B, 
            double dt, double beta, double gamma, uint valueId)
        {
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

            if (World.GetPortCount(pQuantityId) == 0)
            {
                return;
            }
            FieldValue FV = null;
            if (valueId != 0)
            {
                // Time Domain
                FV = World.GetFieldValue(valueId);
            }
            uint portCnt = World.GetPortCount(pQuantityId);
            IList<PortCondition> portConditions = World.GetPortConditions(pQuantityId);
            IList<uint> feIds = World.GetLineFEIds(pQuantityId);
            for (int portId = 0; portId < portCnt; portId++)
            {
                PortCondition portCondition = portConditions[portId];
                IList<int> intParam = portCondition.IntAdditionalParameters;
                System.Diagnostics.Debug.Assert(intParam.Count == 1);
                FlowPressureBCType bcType = (FlowPressureBCType)intParam[0];
                if (bcType != FlowPressureBCType.NoConstraint)
                {
                    continue;
                }
                IList<uint> bcEIds = portCondition.EIds;

                foreach (uint feId in feIds)
                {
                    LineFE lineFE = World.GetLineFE(pQuantityId, feId);
                    uint meshId = lineFE.MeshId;
                    //int meshElemId = lineFE.MeshElemId;
                    uint eId;
                    {
                        uint elemCount;
                        MeshType meshType;
                        int loc;
                        uint cadId;
                        World.Mesh.GetMeshInfo(meshId, out elemCount, out meshType, out loc, out cadId);
                        System.Diagnostics.Debug.Assert(meshType == MeshType.Bar);
                        eId = cadId;
                    }
                    if (bcEIds.Contains(eId))
                    {
                        // BCを適用する辺
                    }
                    else
                    {
                        continue;
                    }

                    // BC
                    uint elemNodeCnt = lineFE.NodeCount;
                    int[] pCoIds = lineFE.NodeCoordIds;
                    int[] vNodes = new int[elemNodeCnt];
                    int[] pNodes = new int[elemNodeCnt];
                    int[] pxNodes = new int[elemNodeCnt];
                    int[] pyNodes = new int[elemNodeCnt];
                    int[] pxxNodes = new int[elemNodeCnt];
                    int[] pxyNodes = new int[elemNodeCnt];
                    int[] pyyNodes = new int[elemNodeCnt];
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int coId = pCoIds[iNode];
                        vNodes[iNode] = World.Coord2Node(vQuantityId, coId);
                        pNodes[iNode] = World.Coord2Node(pQuantityId, coId);
                        pxNodes[iNode] = World.Coord2Node(pxQuantityId, coId);
                        pyNodes[iNode] = World.Coord2Node(pyQuantityId, coId);
                        pxxNodes[iNode] = World.Coord2Node(pxxQuantityId, coId);
                        pxyNodes[iNode] = World.Coord2Node(pxyQuantityId, coId);
                        pyyNodes[iNode] = World.Coord2Node(pyyQuantityId, coId);
                    }

                    Material ma0 = World.GetMaterial(lineFE.MaterialId);
                    System.Diagnostics.Debug.Assert(ma0 is NewtonFluidMaterial);
                    var ma = ma0 as NewtonFluidMaterial;
                    double rho = ma.MassDensity;
                    double mu = ma.Mu;
                    double[] g = { ma.GravityX, ma.GravityY };

                    double[] normal = lineFE.GetNormal();
                    //double[,] sNN = lineFE.CalcSNN();

                    {
                        int coId1 = pCoIds[0];
                        int coId2 = pCoIds[1];
                        IList<uint> triFEIds = World.GetTriangleFEIdsFromEdgeCoord(pQuantityId, coId1, coId2);
                        uint triFEId = triFEIds[0];
                        TriangleFE pTriFE = World.GetTriangleFE(pQuantityId, triFEId);
                        uint pTriFENodeCnt = pTriFE.NodeCount;
                        int[] pTriFECoIds = pTriFE.NodeCoordIds;
                        int[] pTriFENodes = new int[pTriFENodeCnt];
                        for (int iNode = 0; iNode < pTriFENodeCnt; iNode++)
                        {
                            int coId = pTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pQuantityId, coId);
                            pTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE pxTriFE = World.GetTriangleFE(pxQuantityId, triFEId);
                        uint pxTriFENodeCnt = pxTriFE.NodeCount;
                        int[] pxTriFECoIds = pxTriFE.NodeCoordIds;
                        int[] pxTriFENodes = new int[pxTriFENodeCnt];
                        for (int iNode = 0; iNode < pxTriFENodeCnt; iNode++)
                        {
                            int coId = pxTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pxQuantityId, coId);
                            pxTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE pyTriFE = World.GetTriangleFE(pyQuantityId, triFEId);
                        uint pyTriFENodeCnt = pyTriFE.NodeCount;
                        int[] pyTriFECoIds = pyTriFE.NodeCoordIds;
                        int[] pyTriFENodes = new int[pyTriFENodeCnt];
                        for (int iNode = 0; iNode < pyTriFENodeCnt; iNode++)
                        {
                            int coId = pyTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pyQuantityId, coId);
                            pyTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE pxxTriFE = World.GetTriangleFE(pxxQuantityId, triFEId);
                        uint pxxTriFENodeCnt = pxxTriFE.NodeCount;
                        int[] pxxTriFECoIds = pxxTriFE.NodeCoordIds;
                        int[] pxxTriFENodes = new int[pxxTriFENodeCnt];
                        for (int iNode = 0; iNode < pxxTriFENodeCnt; iNode++)
                        {
                            int coId = pxxTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pxxQuantityId, coId);
                            pxxTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE pxyTriFE = World.GetTriangleFE(pxyQuantityId, triFEId);
                        uint pxyTriFENodeCnt = pxyTriFE.NodeCount;
                        int[] pxyTriFECoIds = pxyTriFE.NodeCoordIds;
                        int[] pxyTriFENodes = new int[pxyTriFENodeCnt];
                        for (int iNode = 0; iNode < pxyTriFENodeCnt; iNode++)
                        {
                            int coId = pxyTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pxyQuantityId, coId);
                            pxyTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE pyyTriFE = World.GetTriangleFE(pyyQuantityId, triFEId);
                        uint pyyTriFENodeCnt = pyyTriFE.NodeCount;
                        int[] pyyTriFECoIds = pyyTriFE.NodeCoordIds;
                        int[] pyyTriFENodes = new int[pyyTriFENodeCnt];
                        for (int iNode = 0; iNode < pyyTriFENodeCnt; iNode++)
                        {
                            int coId = pyyTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pyyQuantityId, coId);
                            pyyTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE vTriFE = World.GetTriangleFE(vQuantityId, triFEId);
                        uint vTriFENodeCnt = vTriFE.NodeCount;
                        int[] vTriFECoIds = vTriFE.NodeCoordIds;
                        int[] vTriFENodes = new int[vTriFENodeCnt];
                        for (int iNode = 0; iNode < vTriFENodeCnt; iNode++)
                        {
                            int coId = vTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(vQuantityId, coId);
                            vTriFENodes[iNode] = nodeId;
                        }

                        IntegrationPoints ip = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);
                        System.Diagnostics.Debug.Assert(ip.Ls.Length == 5);
                        for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                        {
                            double[] L = ip.Ls[ipPt];
                            double lineLen = lineFE.GetLineLength();
                            double weight = ip.Weights[ipPt];
                            double detJWeight = (lineLen / 2.0) * weight;

                            double[] pt = lineFE.L2Coord(L);

                            double[] v = new double[vDof];
                            double[] vL = vTriFE.Coord2L(pt);
                            double[] vN = vTriFE.CalcN(vL);
                            double[][] vNu = vTriFE.CalcNu(vL);
                            double[] vNx = vNu[0];
                            double[] vNy = vNu[1];
                            double[,][] vNuv = vTriFE.CalcNuv(vL);
                            for (int ivNode = 0; ivNode < vTriFENodeCnt; ivNode++)
                            {
                                int vNodeId = vTriFENodes[ivNode];
                                if (vNodeId == -1)
                                {
                                    continue;
                                }
                                for (int iDof = 0; iDof < vDof; iDof++)
                                {
                                    double vValue = U[vNodeId * vDof + iDof];
                                    v[iDof] += vValue * vN[ivNode];
                                }
                            }
                            double[] pL = pTriFE.Coord2L(pt);
                            double[] pN = pTriFE.CalcN(pL);
                            double[] pxN = pxTriFE.CalcN(pL);
                            double[] pyN = pyTriFE.CalcN(pL);
                            double[] pxxN = pxxTriFE.CalcN(pL);
                            double[] pxyN = pxyTriFE.CalcN(pL);
                            double[] pyyN = pyyTriFE.CalcN(pL);
                            double[][] pNu = pTriFE.CalcNu(pL);
                            double[] pNx = pNu[0];
                            double[] pNy = pNu[1];
                            double[][] pxNu = pxTriFE.CalcNu(pL);
                            double[] pxNx = pxNu[0];
                            double[] pxNy = pxNu[1];
                            double[][] pyNu = pyTriFE.CalcNu(pL);
                            double[] pyNx = pyNu[0];
                            double[] pyNy = pyNu[1];
                            double[][] pxxNu = pxxTriFE.CalcNu(pL);
                            double[] pxxNx = pxxNu[0];
                            double[] pxxNy = pxxNu[1];
                            double[][] pxyNu = pxyTriFE.CalcNu(pL);
                            double[] pxyNx = pxyNu[0];
                            double[] pxyNy = pxyNu[1];
                            double[][] pyyNu = pyyTriFE.CalcNu(pL);
                            double[] pyyNx = pyyNu[0];
                            double[] pyyNy = pyyNu[1];
                            int[] offsetps = {
                                offsetp, offsetpx, offsetpy,
                                offsetpxx, offsetpxy, offsetpyy
                            };
                            uint[] pTriFENodeCnts = {
                                pTriFENodeCnt, pxTriFENodeCnt, pyTriFENodeCnt,
                                pxxTriFENodeCnt, pxyTriFENodeCnt, pyyTriFENodeCnt
                            };
                            int[][] pTriFENodess = {
                                pTriFENodes, pxTriFENodes, pyTriFENodes,
                                pxxTriFENodes, pxyTriFENodes, pyyTriFENodes
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

                            // kpv
                            for (int rowQuantityId = 0; rowQuantityId < 6; rowQuantityId++)
                            {
                                int rowoffsetp = offsetps[rowQuantityId];
                                uint rowpTriFENodeCnt = pTriFENodeCnts[rowQuantityId];
                                int[] rowpTriFENodes = pTriFENodess[rowQuantityId];
                                int[] rowpNodes = pNodess[rowQuantityId];
                                double[] rowpN = pNs[rowQuantityId];
                                double[] rowpNx = pNxs[rowQuantityId];
                                double[] rowpNy = pNys[rowQuantityId];

                                for (int row = 0; row < rowpTriFENodeCnt; row++)
                                {
                                    int rowNodeId = rowpTriFENodes[row];
                                    if (rowNodeId == -1)
                                    {
                                        continue;
                                    }
                                    // 境界上の節点のみ
                                    if (!rowpNodes.Contains(rowNodeId))
                                    {
                                        continue;
                                    }
                                    for (int col = 0; col < vTriFENodeCnt; col++)
                                    {
                                        int colNodeId = vTriFENodes[col];
                                        if (colNodeId == -1)
                                        {
                                            continue;
                                        }

                                        double[,] kpv1 = new double[pDof, vDof];
                                        kpv1[0, 0] = detJWeight * mu * (-vNy[col]) * (
                                            normal[1] * rowpNx[row] - normal[0] * rowpNy[row]);
                                        kpv1[0, 1] = detJWeight * mu * vNx[col] * (
                                            normal[1] * rowpNx[row] - normal[0] * rowpNy[row]);

                                        double[,] kpv2 = new double[pDof, vDof];
                                        if (valueId != 0)
                                        {
                                            // Time Domain
                                            kpv2[0, 0] = detJWeight * rho * rowpN[row] * vN[col] *
                                                normal[0] * (gamma / (beta * dt));
                                            kpv2[0, 1] = detJWeight * rho * pN[row] * vN[col] *
                                                normal[1] * (gamma / (beta * dt));
                                        }

                                        for (int colDof = 0; colDof < vDof; colDof++)
                                        {
                                            A[rowoffsetp + rowNodeId, colNodeId * vDof + colDof] +=
                                                kpv1[0, colDof] + kpv2[0, colDof];
                                        }
                                    }
                                }
                            }
                            // qp
                            for (int rowQuantityId = 0; rowQuantityId < 6; rowQuantityId++)
                            {
                                int rowoffsetp = offsetps[rowQuantityId];
                                uint rowpTriFENodeCnt = pTriFENodeCnts[rowQuantityId];
                                int[] rowpTriFENodes = pTriFENodess[rowQuantityId];
                                int[] rowpNodes = pNodess[rowQuantityId];
                                double[] rowpN = pNs[rowQuantityId];
                                double[] rowpNx = pNxs[rowQuantityId];
                                double[] rowpNy = pNys[rowQuantityId];

                                for (int row = 0; row < rowpTriFENodeCnt; row++)
                                {
                                    int rowNodeId = rowpTriFENodes[row];
                                    if (rowNodeId == -1)
                                    {
                                        continue;
                                    }
                                    // 境界上の節点のみ
                                    if (!rowpNodes.Contains(rowNodeId))
                                    {
                                        continue;
                                    }
                                    double qp1 = 0;
                                    double qp2 = 0;
                                    if (valueId != 0)
                                    {
                                        // Time Domain
                                        for (int col = 0; col < vTriFENodeCnt; col++)
                                        {
                                            int colCoId = vTriFECoIds[col];
                                            double[] u = FV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                                            double[] vel = FV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                                            double[] acc = FV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);
                                            for (int colDof = 0; colDof < vDof; colDof++)
                                            {
                                                qp2 += detJWeight * rho * rowpN[row] * vN[col] *
                                                    normal[colDof] * (
                                                    -(gamma / (beta * dt)) * u[colDof] +
                                                    (1.0 - (gamma / beta)) * vel[colDof] +
                                                    dt * (1.0 - (gamma / (2.0 * beta))) * acc[colDof]);
                                            }
                                        }
                                    }

                                    B[rowoffsetp + rowNodeId] += -qp1 - qp2;
                                }
                            }
                        }
                    }
                }
            }
        }

        private void SetPressurePoissonWithBellNormalInflowNeumannBCOfPressure(
            IvyFEM.Linear.DoubleSparseMatrix A, double[] B,
            double dt, double beta, double gamma, uint valueId)
        {
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

            if (World.GetPortCount(pQuantityId) == 0)
            {
                return;
            }
            FieldValue FV = null;
            if (valueId != 0)
            {
                // Time Domain
                FV = World.GetFieldValue(valueId);
            }
            uint portCnt = World.GetPortCount(pQuantityId);
            IList<PortCondition> portConditions = World.GetPortConditions(pQuantityId);
            IList<uint> feIds = World.GetLineFEIds(pQuantityId);
            for (int portId = 0; portId < portCnt; portId++)
            {
                PortCondition portCondition = portConditions[portId];
                IList<int> intParam = portCondition.IntAdditionalParameters;
                System.Diagnostics.Debug.Assert(intParam.Count == 1);
                FlowPressureBCType bcType = (FlowPressureBCType)intParam[0];
                if (bcType != FlowPressureBCType.NormalInflow)
                {
                    continue;
                }
                IList<uint> bcEIds = portCondition.EIds;

                foreach (uint feId in feIds)
                {
                    LineFE lineFE = World.GetLineFE(pQuantityId, feId);
                    uint meshId = lineFE.MeshId;
                    //int meshElemId = lineFE.MeshElemId;
                    uint eId;
                    {
                        uint elemCount;
                        MeshType meshType;
                        int loc;
                        uint cadId;
                        World.Mesh.GetMeshInfo(meshId, out elemCount, out meshType, out loc, out cadId);
                        System.Diagnostics.Debug.Assert(meshType == MeshType.Bar);
                        eId = cadId;
                    }
                    if (bcEIds.Contains(eId))
                    {
                        // BCを適用する辺
                    }
                    else
                    {
                        continue;
                    }

                    // BC
                    uint elemNodeCnt = lineFE.NodeCount;
                    int[] pCoIds = lineFE.NodeCoordIds;
                    int[] vNodes = new int[elemNodeCnt];
                    int[] pNodes = new int[elemNodeCnt];
                    int[] pxNodes = new int[elemNodeCnt];
                    int[] pyNodes = new int[elemNodeCnt];
                    int[] pxxNodes = new int[elemNodeCnt];
                    int[] pxyNodes = new int[elemNodeCnt];
                    int[] pyyNodes = new int[elemNodeCnt];
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int coId = pCoIds[iNode];
                        vNodes[iNode] = World.Coord2Node(vQuantityId, coId);
                        pNodes[iNode] = World.Coord2Node(pQuantityId, coId);
                        pxNodes[iNode] = World.Coord2Node(pxQuantityId, coId);
                        pyNodes[iNode] = World.Coord2Node(pyQuantityId, coId);
                        pxxNodes[iNode] = World.Coord2Node(pxxQuantityId, coId);
                        pxyNodes[iNode] = World.Coord2Node(pxyQuantityId, coId);
                        pyyNodes[iNode] = World.Coord2Node(pyyQuantityId, coId);
                    }

                    Material ma0 = World.GetMaterial(lineFE.MaterialId);
                    System.Diagnostics.Debug.Assert(ma0 is NewtonFluidMaterial);
                    var ma = ma0 as NewtonFluidMaterial;
                    double rho = ma.MassDensity;
                    double mu = ma.Mu;
                    double[] g = { ma.GravityX, ma.GravityY };

                    double[] normal = lineFE.GetNormal();
                    //double[,] sNN = lineFE.CalcSNN();

                    {
                        int coId1 = pCoIds[0];
                        int coId2 = pCoIds[1];
                        IList<uint> triFEIds = World.GetTriangleFEIdsFromEdgeCoord(pQuantityId, coId1, coId2);
                        uint triFEId = triFEIds[0];
                        TriangleFE pTriFE = World.GetTriangleFE(pQuantityId, triFEId);
                        uint pTriFENodeCnt = pTriFE.NodeCount;
                        int[] pTriFECoIds = pTriFE.NodeCoordIds;
                        int[] pTriFENodes = new int[pTriFENodeCnt];
                        for (int iNode = 0; iNode < pTriFENodeCnt; iNode++)
                        {
                            int coId = pTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pQuantityId, coId);
                            pTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE pxTriFE = World.GetTriangleFE(pxQuantityId, triFEId);
                        uint pxTriFENodeCnt = pxTriFE.NodeCount;
                        int[] pxTriFECoIds = pxTriFE.NodeCoordIds;
                        int[] pxTriFENodes = new int[pxTriFENodeCnt];
                        for (int iNode = 0; iNode < pxTriFENodeCnt; iNode++)
                        {
                            int coId = pxTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pxQuantityId, coId);
                            pxTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE pyTriFE = World.GetTriangleFE(pyQuantityId, triFEId);
                        uint pyTriFENodeCnt = pyTriFE.NodeCount;
                        int[] pyTriFECoIds = pyTriFE.NodeCoordIds;
                        int[] pyTriFENodes = new int[pyTriFENodeCnt];
                        for (int iNode = 0; iNode < pyTriFENodeCnt; iNode++)
                        {
                            int coId = pyTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pyQuantityId, coId);
                            pyTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE pxxTriFE = World.GetTriangleFE(pxxQuantityId, triFEId);
                        uint pxxTriFENodeCnt = pxxTriFE.NodeCount;
                        int[] pxxTriFECoIds = pxxTriFE.NodeCoordIds;
                        int[] pxxTriFENodes = new int[pxxTriFENodeCnt];
                        for (int iNode = 0; iNode < pxxTriFENodeCnt; iNode++)
                        {
                            int coId = pxxTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pxxQuantityId, coId);
                            pxxTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE pxyTriFE = World.GetTriangleFE(pxyQuantityId, triFEId);
                        uint pxyTriFENodeCnt = pxyTriFE.NodeCount;
                        int[] pxyTriFECoIds = pxyTriFE.NodeCoordIds;
                        int[] pxyTriFENodes = new int[pxyTriFENodeCnt];
                        for (int iNode = 0; iNode < pxyTriFENodeCnt; iNode++)
                        {
                            int coId = pxyTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pxyQuantityId, coId);
                            pxyTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE pyyTriFE = World.GetTriangleFE(pyyQuantityId, triFEId);
                        uint pyyTriFENodeCnt = pyyTriFE.NodeCount;
                        int[] pyyTriFECoIds = pyyTriFE.NodeCoordIds;
                        int[] pyyTriFENodes = new int[pyyTriFENodeCnt];
                        for (int iNode = 0; iNode < pyyTriFENodeCnt; iNode++)
                        {
                            int coId = pyyTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pyyQuantityId, coId);
                            pyyTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE vTriFE = World.GetTriangleFE(vQuantityId, triFEId);
                        uint vTriFENodeCnt = vTriFE.NodeCount;
                        int[] vTriFECoIds = vTriFE.NodeCoordIds;
                        int[] vTriFENodes = new int[vTriFENodeCnt];
                        for (int iNode = 0; iNode < vTriFENodeCnt; iNode++)
                        {
                            int coId = vTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(vQuantityId, coId);
                            vTriFENodes[iNode] = nodeId;
                        }

                        IntegrationPoints ip = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);
                        System.Diagnostics.Debug.Assert(ip.Ls.Length == 5);
                        for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                        {
                            double[] L = ip.Ls[ipPt];
                            double lineLen = lineFE.GetLineLength();
                            double weight = ip.Weights[ipPt];
                            double detJWeight = (lineLen / 2.0) * weight;

                            double[] pt = lineFE.L2Coord(L);

                            double[] v = new double[vDof];
                            double[] vL = vTriFE.Coord2L(pt);
                            double[] vN = vTriFE.CalcN(vL);
                            double[][] vNu = vTriFE.CalcNu(vL);
                            double[] vNx = vNu[0];
                            double[] vNy = vNu[1];
                            double[,][] vNuv = vTriFE.CalcNuv(vL);
                            for (int ivNode = 0; ivNode < vTriFENodeCnt; ivNode++)
                            {
                                int vNodeId = vTriFENodes[ivNode];
                                if (vNodeId == -1)
                                {
                                    continue;
                                }
                                for (int iDof = 0; iDof < vDof; iDof++)
                                {
                                    double vValue = U[vNodeId * vDof + iDof];
                                    v[iDof] += vValue * vN[ivNode];
                                }
                            }
                            double[] pL = pTriFE.Coord2L(pt);
                            double[] pN = pTriFE.CalcN(pL);
                            double[] pxN = pxTriFE.CalcN(pL);
                            double[] pyN = pyTriFE.CalcN(pL);
                            double[] pxxN = pxxTriFE.CalcN(pL);
                            double[] pxyN = pxyTriFE.CalcN(pL);
                            double[] pyyN = pyyTriFE.CalcN(pL);
                            double[][] pNu = pTriFE.CalcNu(pL);
                            double[] pNx = pNu[0];
                            double[] pNy = pNu[1];
                            double[][] pxNu = pxTriFE.CalcNu(pL);
                            double[] pxNx = pxNu[0];
                            double[] pxNy = pxNu[1];
                            double[][] pyNu = pyTriFE.CalcNu(pL);
                            double[] pyNx = pyNu[0];
                            double[] pyNy = pyNu[1];
                            double[][] pxxNu = pxxTriFE.CalcNu(pL);
                            double[] pxxNx = pxxNu[0];
                            double[] pxxNy = pxxNu[1];
                            double[][] pxyNu = pxyTriFE.CalcNu(pL);
                            double[] pxyNx = pxyNu[0];
                            double[] pxyNy = pxyNu[1];
                            double[][] pyyNu = pyyTriFE.CalcNu(pL);
                            double[] pyyNx = pyyNu[0];
                            double[] pyyNy = pyyNu[1];
                            int[] offsetps = {
                                offsetp, offsetpx, offsetpy,
                                offsetpxx, offsetpxy, offsetpyy
                            };
                            uint[] pTriFENodeCnts = {
                                pTriFENodeCnt, pxTriFENodeCnt, pyTriFENodeCnt,
                                pxxTriFENodeCnt, pxyTriFENodeCnt, pyyTriFENodeCnt
                            };
                            int[][] pTriFENodess = {
                                pTriFENodes, pxTriFENodes, pyTriFENodes,
                                pxxTriFENodes, pxyTriFENodes, pyyTriFENodes
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

                            // kpv
                            for (int rowQuantityId = 0; rowQuantityId < 6; rowQuantityId++)
                            {
                                int rowoffsetp = offsetps[rowQuantityId];
                                uint rowpTriFENodeCnt = pTriFENodeCnts[rowQuantityId];
                                int[] rowpTriFENodes = pTriFENodess[rowQuantityId];
                                int[] rowpNodes = pNodess[rowQuantityId];
                                double[] rowpN = pNs[rowQuantityId];
                                double[] rowpNx = pNxs[rowQuantityId];
                                double[] rowpNy = pNys[rowQuantityId];

                                for (int row = 0; row < rowpTriFENodeCnt; row++)
                                {
                                    int rowNodeId = rowpTriFENodes[row];
                                    if (rowNodeId == -1)
                                    {
                                        continue;
                                    }
                                    // 境界上の節点のみ
                                    if (!rowpNodes.Contains(rowNodeId))
                                    {
                                        continue;
                                    }
                                    for (int col = 0; col < vTriFENodeCnt; col++)
                                    {
                                        int colNodeId = vTriFENodes[col];
                                        if (colNodeId == -1)
                                        {
                                            continue;
                                        }

                                        double[,] kpv = new double[pDof, vDof];
                                        kpv[0, 0] = -detJWeight * rho * rowpN[row] *
                                            normal[0] * (v[0] * vNx[col] + v[1] * vNy[col]);
                                        kpv[0, 1] = -detJWeight * rho * rowpN[row] *
                                            normal[1] * (v[0] * vNx[col] + v[1] * vNy[col]);

                                        for (int colDof = 0; colDof < vDof; colDof++)
                                        {
                                            A[rowoffsetp + rowNodeId, colNodeId * vDof + colDof] += kpv[0, colDof];
                                        }
                                    }
                                }
                            }
                            // kpp
                            for (int rowQuantityId = 0; rowQuantityId < 6; rowQuantityId++)
                            {
                                int rowoffsetp = offsetps[rowQuantityId];
                                uint rowpTriFENodeCnt = pTriFENodeCnts[rowQuantityId];
                                int[] rowpTriFENodes = pTriFENodess[rowQuantityId];
                                int[] rowpNodes = pNodess[rowQuantityId];
                                double[] rowpN = pNs[rowQuantityId];
                                double[] rowpNx = pNxs[rowQuantityId];
                                double[] rowpNy = pNys[rowQuantityId];

                                for (int row = 0; row < rowpTriFENodeCnt; row++)
                                {
                                    int rowNodeId = rowpTriFENodes[row];
                                    if (rowNodeId == -1)
                                    {
                                        continue;
                                    }
                                    // 境界上の節点のみ
                                    if (!rowpNodes.Contains(rowNodeId))
                                    {
                                        continue;
                                    }

                                    for (int colQuantityId = 0; colQuantityId < 6; colQuantityId++)
                                    {
                                        int coloffsetp = offsetps[colQuantityId];
                                        uint colpTriFENodeCnt = pTriFENodeCnts[colQuantityId];
                                        int[] colpTriFENodes = pTriFENodess[colQuantityId];
                                        int[] colpNodes = pNodess[colQuantityId];
                                        double[] colpN = pNs[colQuantityId];
                                        double[] colpNx = pNxs[colQuantityId];
                                        double[] colpNy = pNys[colQuantityId];

                                        for (int col = 0; col < colpTriFENodeCnt; col++)
                                        {
                                            int colNodeId = colpTriFENodes[col];
                                            if (colNodeId == -1)
                                            {
                                                continue;
                                            }

                                            double[,] kpp = new double[pDof, pDof];
                                            //kpp[0, 0] = -detJWeight * rowpN[row] * (
                                            //    normal[0] * colpNx[col] + normal[1] * colpNy[col]);
                                            // dp/dn = 0

                                            A[rowoffsetp + rowNodeId, coloffsetp + colNodeId] += kpp[0, 0];
                                        }
                                    }
                                }
                            }
                            // qp
                            for (int rowQuantityId = 0; rowQuantityId < 6; rowQuantityId++)
                            {
                                int rowoffsetp = offsetps[rowQuantityId];
                                uint rowpTriFENodeCnt = pTriFENodeCnts[rowQuantityId];
                                int[] rowpTriFENodes = pTriFENodess[rowQuantityId];
                                int[] rowpNodes = pNodess[rowQuantityId];
                                double[] rowpN = pNs[rowQuantityId];
                                double[] rowpNx = pNxs[rowQuantityId];
                                double[] rowpNy = pNys[rowQuantityId];

                                for (int row = 0; row < rowpTriFENodeCnt; row++)
                                {
                                    int rowNodeId = rowpTriFENodes[row];
                                    if (rowNodeId == -1)
                                    {
                                        continue;
                                    }
                                    // 境界上の節点のみ
                                    if (!rowpNodes.Contains(rowNodeId))
                                    {
                                        continue;
                                    }
                                    double qp = detJWeight * rho * rowpN[row] * (
                                        normal[0] * g[0] + normal[1] * g[1]);

                                    B[rowoffsetp + rowNodeId] += -qp;
                                }
                            }
                        }
                    }
                }
            }
        }

        private void SetPressurePoissonWithBellOutflowNeumannBCOfPressure(
            IvyFEM.Linear.DoubleSparseMatrix A, double[] B,
            double dt, double beta, double gamma, uint valueId)
        {
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

            if (World.GetPortCount(pQuantityId) == 0)
            {
                return;
            }
            FieldValue FV = null;
            if (valueId != 0)
            {
                // Time Domain
                FV = World.GetFieldValue(valueId);
            }
            uint portCnt = World.GetPortCount(pQuantityId);
            IList<PortCondition> portConditions = World.GetPortConditions(pQuantityId);
            IList<uint> feIds = World.GetLineFEIds(pQuantityId);
            for (int portId = 0; portId < portCnt; portId++)
            {
                PortCondition portCondition = portConditions[portId];
                IList<int> intParam = portCondition.IntAdditionalParameters;
                System.Diagnostics.Debug.Assert(intParam.Count == 1);
                FlowPressureBCType bcType = (FlowPressureBCType)intParam[0];
                if (bcType != FlowPressureBCType.Outflow)
                {
                    continue;
                }

                IList<uint> bcEIds = portCondition.EIds;

                foreach (uint feId in feIds)
                {
                    LineFE lineFE = World.GetLineFE(pQuantityId, feId);
                    uint meshId = lineFE.MeshId;
                    //int meshElemId = lineFE.MeshElemId;
                    uint eId;
                    {
                        uint elemCount;
                        MeshType meshType;
                        int loc;
                        uint cadId;
                        World.Mesh.GetMeshInfo(meshId, out elemCount, out meshType, out loc, out cadId);
                        System.Diagnostics.Debug.Assert(meshType == MeshType.Bar);
                        eId = cadId;
                    }
                    if (bcEIds.Contains(eId))
                    {
                        // BCを適用する辺
                    }
                    else
                    {
                        continue;
                    }

                    // BC
                    uint elemNodeCnt = lineFE.NodeCount;
                    int[] pCoIds = lineFE.NodeCoordIds;
                    int[] vNodes = new int[elemNodeCnt];
                    int[] pNodes = new int[elemNodeCnt];
                    int[] pxNodes = new int[elemNodeCnt];
                    int[] pyNodes = new int[elemNodeCnt];
                    int[] pxxNodes = new int[elemNodeCnt];
                    int[] pxyNodes = new int[elemNodeCnt];
                    int[] pyyNodes = new int[elemNodeCnt];
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int coId = pCoIds[iNode];
                        vNodes[iNode] = World.Coord2Node(vQuantityId, coId);
                        pNodes[iNode] = World.Coord2Node(pQuantityId, coId);
                        pxNodes[iNode] = World.Coord2Node(pxQuantityId, coId);
                        pyNodes[iNode] = World.Coord2Node(pyQuantityId, coId);
                        pxxNodes[iNode] = World.Coord2Node(pxxQuantityId, coId);
                        pxyNodes[iNode] = World.Coord2Node(pxyQuantityId, coId);
                        pyyNodes[iNode] = World.Coord2Node(pyyQuantityId, coId);
                    }

                    Material ma0 = World.GetMaterial(lineFE.MaterialId);
                    System.Diagnostics.Debug.Assert(ma0 is NewtonFluidMaterial);
                    var ma = ma0 as NewtonFluidMaterial;
                    double rho = ma.MassDensity;
                    double mu = ma.Mu;
                    double[] g = { ma.GravityX, ma.GravityY };

                    double[] normal = lineFE.GetNormal();
                    double[,] sNN = lineFE.CalcSNN();

                    {
                        int coId1 = pCoIds[0];
                        int coId2 = pCoIds[1];
                        IList<uint> triFEIds = World.GetTriangleFEIdsFromEdgeCoord(pQuantityId, coId1, coId2);
                        uint triFEId = triFEIds[0];
                        TriangleFE pTriFE = World.GetTriangleFE(pQuantityId, triFEId);
                        uint pTriFENodeCnt = pTriFE.NodeCount;
                        int[] pTriFECoIds = pTriFE.NodeCoordIds;
                        int[] pTriFENodes = new int[pTriFENodeCnt];
                        for (int iNode = 0; iNode < pTriFENodeCnt; iNode++)
                        {
                            int coId = pTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pQuantityId, coId);
                            pTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE pxTriFE = World.GetTriangleFE(pxQuantityId, triFEId);
                        uint pxTriFENodeCnt = pxTriFE.NodeCount;
                        int[] pxTriFECoIds = pxTriFE.NodeCoordIds;
                        int[] pxTriFENodes = new int[pxTriFENodeCnt];
                        for (int iNode = 0; iNode < pxTriFENodeCnt; iNode++)
                        {
                            int coId = pxTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pxQuantityId, coId);
                            pxTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE pyTriFE = World.GetTriangleFE(pyQuantityId, triFEId);
                        uint pyTriFENodeCnt = pyTriFE.NodeCount;
                        int[] pyTriFECoIds = pyTriFE.NodeCoordIds;
                        int[] pyTriFENodes = new int[pyTriFENodeCnt];
                        for (int iNode = 0; iNode < pyTriFENodeCnt; iNode++)
                        {
                            int coId = pyTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pyQuantityId, coId);
                            pyTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE pxxTriFE = World.GetTriangleFE(pxxQuantityId, triFEId);
                        uint pxxTriFENodeCnt = pxxTriFE.NodeCount;
                        int[] pxxTriFECoIds = pxxTriFE.NodeCoordIds;
                        int[] pxxTriFENodes = new int[pxxTriFENodeCnt];
                        for (int iNode = 0; iNode < pxxTriFENodeCnt; iNode++)
                        {
                            int coId = pxxTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pxxQuantityId, coId);
                            pxxTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE pxyTriFE = World.GetTriangleFE(pxyQuantityId, triFEId);
                        uint pxyTriFENodeCnt = pxyTriFE.NodeCount;
                        int[] pxyTriFECoIds = pxyTriFE.NodeCoordIds;
                        int[] pxyTriFENodes = new int[pxyTriFENodeCnt];
                        for (int iNode = 0; iNode < pxyTriFENodeCnt; iNode++)
                        {
                            int coId = pxyTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pxyQuantityId, coId);
                            pxyTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE pyyTriFE = World.GetTriangleFE(pyyQuantityId, triFEId);
                        uint pyyTriFENodeCnt = pyyTriFE.NodeCount;
                        int[] pyyTriFECoIds = pyyTriFE.NodeCoordIds;
                        int[] pyyTriFENodes = new int[pyyTriFENodeCnt];
                        for (int iNode = 0; iNode < pyyTriFENodeCnt; iNode++)
                        {
                            int coId = pyyTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(pyyQuantityId, coId);
                            pyyTriFENodes[iNode] = nodeId;
                        }
                        TriangleFE vTriFE = World.GetTriangleFE(vQuantityId, triFEId);
                        uint vTriFENodeCnt = vTriFE.NodeCount;
                        int[] vTriFECoIds = vTriFE.NodeCoordIds;
                        int[] vTriFENodes = new int[vTriFENodeCnt];
                        for (int iNode = 0; iNode < vTriFENodeCnt; iNode++)
                        {
                            int coId = vTriFECoIds[iNode];
                            int nodeId = World.Coord2Node(vQuantityId, coId);
                            vTriFENodes[iNode] = nodeId;
                        }

                        IntegrationPoints ip = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);
                        System.Diagnostics.Debug.Assert(ip.Ls.Length == 5);
                        for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                        {
                            double[] L = ip.Ls[ipPt];
                            double lineLen = lineFE.GetLineLength();
                            double weight = ip.Weights[ipPt];
                            double detJWeight = (lineLen / 2.0) * weight;

                            double[] pt = lineFE.L2Coord(L);

                            double[] v = new double[vDof];
                            double[] vx = new double[vDof];
                            double[] vy = new double[vDof];
                            double[] vL = vTriFE.Coord2L(pt);
                            double[] vN = vTriFE.CalcN(vL);
                            double[][] vNu = vTriFE.CalcNu(vL);
                            double[] vNx = vNu[0];
                            double[] vNy = vNu[1];
                            for (int ivNode = 0; ivNode < vTriFENodeCnt; ivNode++)
                            {
                                int vNodeId = vTriFENodes[ivNode];
                                if (vNodeId == -1)
                                {
                                    continue;
                                }
                                for (int iDof = 0; iDof < vDof; iDof++)
                                {
                                    double vValue = U[vNodeId * vDof + iDof];
                                    v[iDof] += vValue * vN[ivNode];
                                    vx[iDof] += vValue * vNx[ivNode];
                                    vy[iDof] += vValue * vNy[ivNode];
                                }
                            }
                            double[] pL = pTriFE.Coord2L(pt);
                            double[] pN = pTriFE.CalcN(pL);
                            double[] pxN = pxTriFE.CalcN(pL);
                            double[] pyN = pyTriFE.CalcN(pL);
                            double[] pxxN = pxxTriFE.CalcN(pL);
                            double[] pxyN = pxyTriFE.CalcN(pL);
                            double[] pyyN = pyyTriFE.CalcN(pL);
                            double[][] pNu = pTriFE.CalcNu(pL);
                            double[] pNx = pNu[0];
                            double[] pNy = pNu[1];
                            double[][] pxNu = pxTriFE.CalcNu(pL);
                            double[] pxNx = pxNu[0];
                            double[] pxNy = pxNu[1];
                            double[][] pyNu = pyTriFE.CalcNu(pL);
                            double[] pyNx = pyNu[0];
                            double[] pyNy = pyNu[1];
                            double[][] pxxNu = pxxTriFE.CalcNu(pL);
                            double[] pxxNx = pxxNu[0];
                            double[] pxxNy = pxxNu[1];
                            double[][] pxyNu = pxyTriFE.CalcNu(pL);
                            double[] pxyNx = pxyNu[0];
                            double[] pxyNy = pxyNu[1];
                            double[][] pyyNu = pyyTriFE.CalcNu(pL);
                            double[] pyyNx = pyyNu[0];
                            double[] pyyNy = pyyNu[1];
                            int[] offsetps = {
                                offsetp, offsetpx, offsetpy,
                                offsetpxx, offsetpxy, offsetpyy
                            };
                            uint[] pTriFENodeCnts = {
                                pTriFENodeCnt, pxTriFENodeCnt, pyTriFENodeCnt,
                                pxxTriFENodeCnt, pxyTriFENodeCnt, pyyTriFENodeCnt
                            };
                            int[][] pTriFENodess = {
                                pTriFENodes, pxTriFENodes, pyTriFENodes,
                                pxxTriFENodes, pxyTriFENodes, pyyTriFENodes
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

                            // kpv
                            for (int rowQuantityId = 0; rowQuantityId < 6; rowQuantityId++)
                            {
                                int rowoffsetp = offsetps[rowQuantityId];
                                uint rowpTriFENodeCnt = pTriFENodeCnts[rowQuantityId];
                                int[] rowpTriFENodes = pTriFENodess[rowQuantityId];
                                int[] rowpNodes = pNodess[rowQuantityId];
                                double[] rowpN = pNs[rowQuantityId];
                                double[] rowpNx = pNxs[rowQuantityId];
                                double[] rowpNy = pNys[rowQuantityId];

                                for (int row = 0; row < rowpTriFENodeCnt; row++)
                                {
                                    int rowNodeId = rowpTriFENodes[row];
                                    if (rowNodeId == -1)
                                    {
                                        continue;
                                    }
                                    // 境界上の節点のみ
                                    if (!rowpNodes.Contains(rowNodeId))
                                    {
                                        continue;
                                    }
                                    for (int col = 0; col < vTriFENodeCnt; col++)
                                    {
                                        int colNodeId = vTriFENodes[col];
                                        if (colNodeId == -1)
                                        {
                                            continue;
                                        }

                                        double[,] kpv = new double[pDof, vDof];
                                        if (valueId != 0)
                                        {
                                            // Time Domain
                                            kpv[0, 0] = detJWeight * rho * rowpN[row] * vN[col] *
                                                normal[0] * (gamma / (beta * dt));
                                            kpv[0, 1] = detJWeight * rho * rowpN[row] * vN[col] *
                                                normal[1] * (gamma / (beta * dt));
                                        }

                                        for (int colDof = 0; colDof < vDof; colDof++)
                                        {
                                            A[rowoffsetp + rowNodeId, colNodeId * vDof + colDof] += kpv[0, colDof];
                                        }
                                    }
                                }
                            }
                            // kpp
                            for (int rowQuantityId = 0; rowQuantityId < 6; rowQuantityId++)
                            {
                                int rowoffsetp = offsetps[rowQuantityId];
                                uint rowpTriFENodeCnt = pTriFENodeCnts[rowQuantityId];
                                int[] rowpTriFENodes = pTriFENodess[rowQuantityId];
                                int[] rowpNodes = pNodess[rowQuantityId];
                                double[] rowpN = pNs[rowQuantityId];
                                double[] rowpNx = pNxs[rowQuantityId];
                                double[] rowpNy = pNys[rowQuantityId];

                                for (int row = 0; row < rowpTriFENodeCnt; row++)
                                {
                                    int rowNodeId = rowpTriFENodes[row];
                                    if (rowNodeId == -1)
                                    {
                                        continue;
                                    }
                                    // 境界上の節点のみ
                                    if (!rowpNodes.Contains(rowNodeId))
                                    {
                                        continue;
                                    }

                                    for (int colQuantityId = 0; colQuantityId < 6; colQuantityId++)
                                    {
                                        int coloffsetp = offsetps[colQuantityId];
                                        uint colpTriFENodeCnt = pTriFENodeCnts[colQuantityId];
                                        int[] colpTriFENodes = pTriFENodess[colQuantityId];
                                        int[] colpNodes = pNodess[colQuantityId];
                                        double[] colpN = pNs[colQuantityId];
                                        double[] colpNx = pNxs[colQuantityId];
                                        double[] colpNy = pNys[colQuantityId];

                                        for (int col = 0; col < colpTriFENodeCnt; col++)
                                        {
                                            int colNodeId = colpTriFENodes[col];
                                            if (colNodeId == -1)
                                            {
                                                continue;
                                            }

                                            double[,] kpp = new double[pDof, pDof];
                                            kpp[0, 0] = -detJWeight * rowpN[row] * (
                                                normal[0] * colpNx[col] + normal[1] * colpNy[col]);

                                            A[rowoffsetp + rowNodeId, coloffsetp + colNodeId] += kpp[0, 0];
                                        }
                                    }
                                }
                            }
                            // qp
                            for (int rowQuantityId = 0; rowQuantityId < 6; rowQuantityId++)
                            {
                                int rowoffsetp = offsetps[rowQuantityId];
                                uint rowpTriFENodeCnt = pTriFENodeCnts[rowQuantityId];
                                int[] rowpTriFENodes = pTriFENodess[rowQuantityId];
                                int[] rowpNodes = pNodess[rowQuantityId];
                                double[] rowpN = pNs[rowQuantityId];
                                double[] rowpNx = pNxs[rowQuantityId];
                                double[] rowpNy = pNys[rowQuantityId];

                                for (int row = 0; row < rowpTriFENodeCnt; row++)
                                {
                                    int rowNodeId = rowpTriFENodes[row];
                                    if (rowNodeId == -1)
                                    {
                                        continue;
                                    }
                                    // 境界上の節点のみ
                                    if (!rowpNodes.Contains(rowNodeId))
                                    {
                                        continue;
                                    }
                                    double qp1 = 0;
                                    double qp2 = 0;
                                    if (valueId != 0)
                                    {
                                        // Time Domain
                                        for (int col = 0; col < vTriFENodeCnt; col++)
                                        {
                                            int colCoId = vTriFECoIds[col];
                                            double[] u = FV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                                            double[] vel = FV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                                            double[] acc = FV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);
                                            for (int colDof = 0; colDof < vDof; colDof++)
                                            {
                                                qp2 += detJWeight * rho * rowpN[row] * vN[col] *
                                                    normal[colDof] * (
                                                    -(gamma / (beta * dt)) * u[colDof] +
                                                    (1.0 - (gamma / beta)) * vel[colDof] +
                                                    dt * (1.0 - (gamma / (2.0 * beta))) * acc[colDof]);
                                            }
                                        }
                                    }

                                    B[rowoffsetp + rowNodeId] += -qp1 - qp2;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
