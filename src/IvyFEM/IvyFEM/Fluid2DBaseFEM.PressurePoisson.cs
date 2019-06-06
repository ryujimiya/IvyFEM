using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class Fluid2DBaseFEM
    {
        protected void SetPressurePoissonSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B,
            double dt, double beta, double gamma, uint valueId = 0)
        {
            // for fixed velocity
            SetPressurePoissonNoConstraintNeumannBCOfPressure(A, B, dt, beta, gamma, valueId);
            SetPressurePoissonNormalInflowNeumannBCOfPressure(A, B, dt, beta, gamma, valueId);
            // for outflow
            SetPressurePoissonOutflowNeumannBCOfPressure(A, B, dt, beta, gamma, valueId);
        }

        private void SetPressurePoissonNoConstraintNeumannBCOfPressure(
            IvyFEM.Linear.DoubleSparseMatrix A, double[] B, 
            double dt, double beta, double gamma, uint valueId)
        {
            uint vQuantityId = 0;
            uint pQuantityId = 1;
            int vDof = 2;
            int pDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int offset = vNodeCnt * vDof;

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
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int coId = pCoIds[iNode];
                        vNodes[iNode] = World.Coord2Node(vQuantityId, coId);
                        pNodes[iNode] = World.Coord2Node(pQuantityId, coId);
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
                            double[][] pNu = pTriFE.CalcNu(pL);
                            double[] pNx = pNu[0];
                            double[] pNy = pNu[1];

                            // kpv
                            for (int row = 0; row < pTriFENodeCnt; row++)
                            {
                                int rowNodeId = pTriFENodes[row];
                                if (rowNodeId == -1)
                                {
                                    continue;
                                }
                                // 境界上の節点のみ
                                if (!pNodes.Contains(rowNodeId))
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
                                        normal[1] * pNx[row] - normal[0] * pNy[row]);
                                    kpv1[0, 1] = detJWeight * mu * vNx[col] * (
                                        normal[1] * pNx[row] - normal[0] * pNy[row]);

                                    double[,] kpv2 = new double[pDof, vDof];
                                    if (valueId != 0)
                                    {
                                        // Time Domain
                                        kpv2[0, 0] = detJWeight * rho * pN[row] * vN[col] *
                                            normal[0] * (gamma / (beta * dt));
                                        kpv2[0, 1] = detJWeight * rho * pN[row] * vN[col] *
                                            normal[1] * (gamma / (beta * dt));
                                    }

                                    for (int colDof = 0; colDof < vDof; colDof++)
                                    {
                                        A[offset + rowNodeId, colNodeId * vDof + colDof] +=
                                            kpv1[0, colDof] + kpv2[0, colDof];
                                    }
                                }
                            }
                            // qp
                            for (int row = 0; row < pTriFENodeCnt; row++)
                            {
                                int rowNodeId = pTriFENodes[row];
                                if (rowNodeId == -1)
                                {
                                    continue;
                                }
                                // 境界上の節点のみ
                                if (!pNodes.Contains(rowNodeId))
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
                                            qp2 += detJWeight * rho * pN[row] * vN[col] *
                                                normal[colDof] * (
                                                -(gamma / (beta * dt)) * u[colDof] +
                                                (1.0 - (gamma / beta)) * vel[colDof] +
                                                dt * (1.0 - (gamma / (2.0 * beta))) * acc[colDof]);
                                        }
                                    }
                                }

                                B[offset + rowNodeId] += -qp1 -qp2;
                            }
                        }
                    }
                }
            }
        }

        private void SetPressurePoissonNormalInflowNeumannBCOfPressure(
            IvyFEM.Linear.DoubleSparseMatrix A, double[] B,
            double dt, double beta, double gamma, uint valueId)
        {
            uint vQuantityId = 0;
            uint pQuantityId = 1;
            int vDof = 2;
            int pDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int offset = vNodeCnt * vDof;

            if (World.GetPortCount(pQuantityId) == 0)
            {
                return;
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
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int coId = pCoIds[iNode];
                        vNodes[iNode] = World.Coord2Node(vQuantityId, coId);
                        pNodes[iNode] = World.Coord2Node(pQuantityId, coId);
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
                                    vx[iDof] += vValue * vNx[ivNode];
                                    vy[iDof] += vValue * vNy[ivNode];
                                }
                            }
                            double[] pL = pTriFE.Coord2L(pt);
                            double[] pN = pTriFE.CalcN(pL);
                            double[][] pNu = pTriFE.CalcNu(pL);
                            double[] pNx = pNu[0];
                            double[] pNy = pNu[1];

                            // kpv
                            for (int row = 0; row < pTriFENodeCnt; row++)
                            {
                                int rowNodeId = pTriFENodes[row];
                                if (rowNodeId == -1)
                                {
                                    continue;
                                }
                                // 境界上の節点のみ
                                if (!pNodes.Contains(rowNodeId))
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
                                    kpv[0, 0] = -detJWeight * rho * pN[row] *
                                        normal[0] * (v[0] * vNx[col] + v[1] * vNy[col]);
                                    kpv[0, 1] = -detJWeight * rho * pN[row] *
                                        normal[1] * (v[0] * vNx[col] + v[1] * vNy[col]);

                                    for (int colDof = 0; colDof < vDof; colDof++)
                                    {
                                        A[offset + rowNodeId, colNodeId * vDof + colDof] += kpv[0, colDof];
                                    }
                                }
                            }
                            // kpp
                            for (int row = 0; row < pTriFENodeCnt; row++)
                            {
                                int rowNodeId = pTriFENodes[row];
                                if (rowNodeId == -1)
                                {
                                    continue;
                                }
                                // 境界上の節点のみ
                                if (!pNodes.Contains(rowNodeId))
                                {
                                    continue;
                                }
                                for (int col = 0; col < pTriFENodeCnt; col++)
                                {
                                    int colNodeId = pTriFENodes[col];
                                    if (colNodeId == -1)
                                    {
                                        continue;
                                    }

                                    double[,] kpp = new double[pDof, pDof];
                                    //kpp[0, 0] = -detJWeight * pN[row] * (normal[0] * pNx[col] + normal[1] * pNy[col]);
                                    // dp/dn = 0

                                    A[offset + rowNodeId, offset + colNodeId] += kpp[0, 0];
                                }
                            }
                            // qp
                            for (int row = 0; row < pTriFENodeCnt; row++)
                            {
                                int rowNodeId = pTriFENodes[row];
                                if (rowNodeId == -1)
                                {
                                    continue;
                                }
                                // 境界上の節点のみ
                                if (!pNodes.Contains(rowNodeId))
                                {
                                    continue;
                                }
                                double qp = detJWeight * rho * pN[row] * (normal[0] * g[0] + normal[1] * g[1]);
                                B[offset + rowNodeId] += -qp;
                            }
                        }
                    }
                }
            }
        }

        private void SetPressurePoissonOutflowNeumannBCOfPressure(
            IvyFEM.Linear.DoubleSparseMatrix A, double[] B,
            double dt, double beta, double gamma, uint valueId)
        {
            uint vQuantityId = 0;
            uint pQuantityId = 1;
            int vDof = 2;
            int pDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int offset = vNodeCnt * vDof;

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
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int coId = pCoIds[iNode];
                        vNodes[iNode] = World.Coord2Node(vQuantityId, coId);
                        pNodes[iNode] = World.Coord2Node(pQuantityId, coId);
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
                            double[][] pNu = pTriFE.CalcNu(pL);
                            double[] pNx = pNu[0];
                            double[] pNy = pNu[1];
                            // kpv
                            for (int row = 0; row < pTriFENodeCnt; row++)
                            {
                                int rowNodeId = pTriFENodes[row];
                                if (rowNodeId == -1)
                                {
                                    continue;
                                }
                                // 境界上の節点のみ
                                if (!pNodes.Contains(rowNodeId))
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
                                        kpv[0, 0] = detJWeight * rho * pN[row] * vN[col] *
                                            normal[0] * (gamma / (beta * dt));
                                        kpv[0, 1] = detJWeight * rho * pN[row] * vN[col] *
                                            normal[1] * (gamma / (beta * dt));
                                    }

                                    for (int colDof = 0; colDof < vDof; colDof++)
                                    {
                                        A[offset + rowNodeId, colNodeId * vDof + colDof] += kpv[0, 0];
                                    }
                                }
                            }
                            // kpp
                            for (int row = 0; row < pTriFENodeCnt; row++)
                            {
                                int rowNodeId = pTriFENodes[row];
                                if (rowNodeId == -1)
                                {
                                    continue;
                                }
                                // 境界上の節点のみ
                                if (!pNodes.Contains(rowNodeId))
                                {
                                    continue;
                                }
                                for (int col = 0; col < pTriFENodeCnt; col++)
                                {
                                    int colNodeId = pTriFENodes[col];
                                    if (colNodeId == -1)
                                    {
                                        continue;
                                    }

                                    double[,] kpp = new double[pDof, pDof];
                                    kpp[0, 0] = -detJWeight * pN[row] * (normal[0] * pNx[col] + normal[1] * pNy[col]);

                                    A[offset + rowNodeId, offset + colNodeId] += kpp[0, 0];
                                }
                            }
                            // qp
                            for (int row = 0; row < pTriFENodeCnt; row++)
                            {
                                int rowNodeId = pTriFENodes[row];
                                if (rowNodeId == -1)
                                {
                                    continue;
                                }
                                // 境界上の節点のみ
                                if (!pNodes.Contains(rowNodeId))
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
                                            qp2 += detJWeight * rho * pN[row] * vN[col] * 
                                                normal[colDof] * (
                                                -(gamma / (beta * dt)) * u[colDof] +
                                                (1.0 - (gamma / beta)) * vel[colDof] +
                                                dt * (1.0 - (gamma / (2.0 * beta))) * acc[colDof]);
                                        }
                                    }
                                }

                                B[offset + rowNodeId] += -qp1 - qp2;
                            }
                        }
                    }
                }
            }
        }
    }
}
