using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Fluid2DRKTDFEM
    {
        private void SetPressurePoisson1stEquationSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            // Nothing
        }

        private void SetPressurePoisson2ndEquationSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            // for fixed velocity
            SetPressurePoisson2ndEquationNoConstraintNeumannBCOfPressure(A, B);
            SetPressurePoisson2ndEquationNormalInflowNeumannBCOfPressure(A, B);
            // for outflow
            SetPressurePoisson2ndEquationOutflowNeumannBCOfPressure(A, B);
        }

        private void SetPressurePoisson2ndEquationNoConstraintNeumannBCOfPressure(
            IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
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

                                    for (int colDof = 0; colDof < vDof; colDof++)
                                    {
                                        double vValue = U[colNodeId * vDof + colDof];
                                        B[rowNodeId] +=
                                            -kpv1[0, colDof] * vValue;
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
                                double qp = 0;
                                B[rowNodeId] += -qp;
                            }
                        }
                    }
                }
            }
        }

        private void SetPressurePoisson2ndEquationNormalInflowNeumannBCOfPressure(
            IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
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
                                        double vValue = U[colNodeId * vDof + colDof];
                                        B[rowNodeId] += -kpv[0, colDof] * vValue;
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

                                    A[rowNodeId, colNodeId] += kpp[0, 0];
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
                                B[rowNodeId] += -qp;
                            }
                        }
                    }
                }
            }
        }

        private void SetPressurePoisson2ndEquationOutflowNeumannBCOfPressure(
            IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
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
            uint valueId = ValueId;
            double dt = TimeStep;
            FieldValue FV = null;
            FV = World.GetFieldValue(valueId);

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

                                    // Time Domain
                                    int colCoId = vTriFECoIds[col];
                                    double[] u = FV.GetDoubleValue(colCoId, FieldDerivativeType.Value);

                                    double[,] kpv = new double[pDof, vDof];
                                    kpv[0, 0] = detJWeight * rho * pN[row] * vN[col] *
                                        normal[0];
                                    kpv[0, 1] = detJWeight * rho * pN[row] * vN[col] *
                                        normal[1];

                                    for (int colDof = 0; colDof < vDof; colDof++)
                                    {
                                        double vValue = U[colNodeId * vDof + colDof];
                                        B[rowNodeId] += 
                                            -kpv[0, 0] * (vValue - u[colDof]) / dt;
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

                                    A[rowNodeId, colNodeId] += kpp[0, 0];
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
                                double qp = 0;
                                B[rowNodeId] += -qp;
                            }
                        }
                    }
                }
            }
        }
    }
}
