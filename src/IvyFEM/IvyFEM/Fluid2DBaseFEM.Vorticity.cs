using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class Fluid2DBaseFEM
    {
        protected void SetVorticitySpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            SetVorticityDirichletBCOfVorticityForTangentialFlow(A, B);
            SetVorticityDirichletBCOfVorticityForOutflow(A, B);
            SetVorticityDirichletBCOfStreamForOutflow(A, B);
        }

        protected void VorticityPostSolve()
        {
            GetVelocityFromVorticityStream();
        }

        // 境界が流線方向と同じとき
        // ψの法線成分（速度の接線成分に比例）がある
        // ψ = const or 0
        // Taylor級数展開 2次までの項を考慮
        private void SetVorticityDirichletBCOfVorticityForTangentialFlow(
            IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint wQuantityId = 0;
            uint pQuantityId = 1;
            int wNodeCnt = (int)World.GetNodeCount(wQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int offset = wNodeCnt;

            if (World.GetPortCount(wQuantityId) == 0)
            {
                return;
            }
            uint portCnt = World.GetPortCount(wQuantityId);
            IList<PortCondition> portConditions = World.GetPortConditions(wQuantityId);
            IList<uint> feIds = World.GetLineFEIds(wQuantityId);
            for (int portId = 0; portId < portCnt; portId++)
            {
                PortCondition portCondition = portConditions[portId];
                IList<int> intParam = portCondition.IntAdditionalParameters;
                System.Diagnostics.Debug.Assert(intParam.Count == 1);
                FlowVorticityBCType bcType = (FlowVorticityBCType)intParam[0];
                if (bcType != FlowVorticityBCType.TangentialFlow)
                {
                    continue;
                }
                IList<uint> bcEIds = portCondition.EIds;

                foreach (uint feId in feIds)
                {
                    LineFE lineFE = World.GetLineFE(wQuantityId, feId);
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
                    int[] wCoIds = lineFE.NodeCoordIds;
                    int[] wNodes = new int[elemNodeCnt];
                    int[] pNodes = new int[elemNodeCnt];
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int coId = wCoIds[iNode];
                        wNodes[iNode] = World.Coord2Node(wQuantityId, coId);
                        pNodes[iNode] = World.Coord2Node(pQuantityId, coId);
                    }

                    Material ma0 = World.GetMaterial(lineFE.MaterialId);
                    System.Diagnostics.Debug.Assert(ma0 is NewtonFluidMaterial);
                    var ma = ma0 as NewtonFluidMaterial;
                    double rho = ma.MassDensity;
                    double mu = ma.Mu;
                    double[] g = { ma.GravityX, ma.GravityY };

                    double[] normal = lineFE.GetNormal();

                    int coId1 = wCoIds[0];
                    int coId2 = wCoIds[1];
                    int adjCoId = -1; // adjacent
                    {
                        IList<uint> triFEIds = World.GetTriangleFEIdsFromEdgeCoord(wQuantityId, coId1, coId2);
                        uint triFEId = triFEIds[0];
                        TriangleFE triFE = World.GetTriangleFE(wQuantityId, triFEId);
                        int[] triFECoIds = triFE.NodeCoordIds;
                        foreach (int coId in triFECoIds)
                        {
                            if (coId != coId1 && coId != coId2)
                            {
                                adjCoId = coId;
                                break;
                            }
                        }
                        System.Diagnostics.Debug.Assert(adjCoId != -1);
                    }
                    //int wAdjNodeId = -1;
                    //wAdjNodeId = World.Coord2Node(wQuantityId, adjCoId); // 特殊な場合-1はありえる
                    int pAdjNodeId = -1;
                    pAdjNodeId = World.Coord2Node(pQuantityId, adjCoId); // 特殊な場合-1はありえる
                    double[] pt1 = World.GetCoord(wQuantityId, coId1);
                    double[] pt2 = World.GetCoord(wQuantityId, coId2);
                    double[] adjPt = World.GetCoord(wQuantityId, adjCoId);
                    double hn = Math.Abs(IvyFEM.CadUtils2D.TriHeight(
                        new OpenTK.Vector2d(adjPt[0], adjPt[1]), // 点
                        new OpenTK.Vector2d(pt1[0], pt1[1]), // 辺の点1
                        new OpenTK.Vector2d(pt2[0], pt2[1]) // 辺の点2
                        ));

                    for (int row = 0; row < elemNodeCnt; row++)
                    {
                        int rowCoId = wCoIds[row];
                        int rowNodeId = wNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        IList<double> param = portCondition.GetDoubleAdditionalParameters(rowCoId);
                        System.Diagnostics.Debug.Assert(param.Count == 2); // 0: dψ/dx 1: dψ/dy
                        double pxValue = param[0];
                        double pyValue = param[1];

                        for (int colNodeId = 0; colNodeId < wNodeCnt; colNodeId++)
                        {
                            int colCoId = World.Node2Coord(wQuantityId, colNodeId);
                            if (colCoId == rowCoId)
                            {
                                A[rowNodeId, colNodeId] = 1.0;
                            }
                            else
                            {
                                A[rowNodeId, colNodeId] = 0;
                            }
                        }
                        for (int colNodeId = 0; colNodeId < pNodeCnt; colNodeId++)
                        {
                            int colCoId = World.Node2Coord(pQuantityId, colNodeId);
                            if (colCoId == rowCoId)
                            {
                                A[rowNodeId, offset + colNodeId] = -2.0 / (hn * hn);
                            }
                            else if (pAdjNodeId != -1 && colNodeId == pAdjNodeId)
                            {
                                A[rowNodeId, offset + colNodeId] = 2.0 / (hn * hn);
                            }
                            else
                            {
                                A[rowNodeId, offset + colNodeId] = 0;
                            }
                        }
                        B[rowNodeId] = -(normal[0] * pxValue + normal[1] * pyValue) * (2.0 / hn);
                    }
                }
            }
        }

        private void SetVorticityDirichletBCOfVorticityForOutflow(
            IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint wQuantityId = 0;
            uint pQuantityId = 1;
            int wNodeCnt = (int)World.GetNodeCount(wQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int offset = wNodeCnt;
            int vDof = 2; // 速度

            if (World.GetPortCount(wQuantityId) == 0)
            {
                return;
            }
            uint portCnt = World.GetPortCount(wQuantityId);
            IList<PortCondition> portConditions = World.GetPortConditions(wQuantityId);
            IList<uint> feIds = World.GetLineFEIds(wQuantityId);
            for (int portId = 0; portId < portCnt; portId++)
            {
                PortCondition portCondition = portConditions[portId];
                IList<int> intParam = portCondition.IntAdditionalParameters;
                System.Diagnostics.Debug.Assert(intParam.Count == 1);
                FlowVorticityBCType bcType = (FlowVorticityBCType)intParam[0];
                if (bcType != FlowVorticityBCType.Outflow)
                {
                    continue;
                }
                IList<uint> bcEIds = portCondition.EIds;
                //IList<double> param = portCondition.DoubleAdditionalParameters;
                //System.Diagnostics.Debug.Assert(param.Count == 0);

                foreach (uint feId in feIds)
                {
                    LineFE lineFE = World.GetLineFE(wQuantityId, feId);
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
                    int[] wCoIds = lineFE.NodeCoordIds;
                    int[] wNodes = new int[elemNodeCnt];
                    int[] pNodes = new int[elemNodeCnt];
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int coId = wCoIds[iNode];
                        wNodes[iNode] = World.Coord2Node(wQuantityId, coId);
                        pNodes[iNode] = World.Coord2Node(pQuantityId, coId);
                    }

                    Material ma0 = World.GetMaterial(lineFE.MaterialId);
                    System.Diagnostics.Debug.Assert(ma0 is NewtonFluidMaterial);
                    var ma = ma0 as NewtonFluidMaterial;
                    double rho = ma.MassDensity;
                    double mu = ma.Mu;
                    double[] g = { ma.GravityX, ma.GravityY };

                    double[] normal = lineFE.GetNormal();
                    double[] tan = { normal[1], -normal[0] };

                    int coId1 = wCoIds[0];
                    int coId2 = wCoIds[1];
                    double[] pt1 = World.GetCoord(wQuantityId, coId1);
                    double[] pt2 = World.GetCoord(wQuantityId, coId2);
                    double[] dir = IvyFEM.CadUtils2D.GetDirection2D(pt1, pt2);

                    for (int row = 0; row < elemNodeCnt; row++)
                    {
                        int rowCoId = wCoIds[row];
                        int rowNodeId = wNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        int rowCoId2 = wCoIds[(row + 1) % elemNodeCnt];
                        double[] rowPt = World.GetCoord(wQuantityId, rowCoId);
                        double[] rowPt2 = World.GetCoord(wQuantityId, rowCoId2);
                        double[] rowDir = IvyFEM.CadUtils2D.GetDirection2D(rowPt, rowPt2);
                        double h = dir[0] * (rowPt2[0] - rowPt[0]) + dir[1] * (rowPt2[1] - rowPt[0]); // マイナスもある

                        // ψの補間なのでωの節点にはそのままでは使えない
                        //double[] v = { CoordV[rowCoId * vDof + 0], CoordV[rowCoId * vDof + 1] };
                        double[] v = GetVelocityAtVorticityCoord(rowCoId, lineFE);
                        double vn = normal[0] * v[0] + normal[1] * v[1];

                        for (int colNodeId = 0; colNodeId < wNodeCnt; colNodeId++)
                        {
                            int colCoId = World.Node2Coord(wQuantityId, colNodeId);
                            if (colCoId == rowCoId)
                            {
                                A[rowNodeId, colNodeId] = 1.0;
                            }
                            else
                            {
                                A[rowNodeId, colNodeId] = 0;
                            }
                        }
                        for (int colNodeId = 0; colNodeId < pNodeCnt; colNodeId++)
                        {
                            int colCoId = World.Node2Coord(pQuantityId, colNodeId);
                            if (colCoId == rowCoId)
                            {
                                A[rowNodeId, offset + colNodeId] = -2.0 / (h * h);
                            }
                            else if (colCoId == rowCoId2)
                            {
                                A[rowNodeId, offset + colNodeId] = 2.0 / (h * h);
                            }
                            else
                            {
                                A[rowNodeId, offset + colNodeId] = 0;
                            }
                        }
                        B[rowNodeId] = vn * (2.0 / h);
                    }
                }
            }
        }

        private double[] GetVelocityAtVorticityCoord(int wCoId, LineFE wLineFE)
        {
            uint wQuantityId = 0;
            uint pQuantityId = 1;
            int wCoCnt = (int)World.GetCoordCount(wQuantityId);
            int pCoCnt = (int)World.GetCoordCount(pQuantityId);
            uint vDof = 2;
            double[] v = new double[vDof];
            if (wCoId < pCoCnt)
            {
                v[0] = CoordVelocity[wCoId * vDof];
                v[1] = CoordVelocity[wCoId * vDof + 1];
            }
            else
            {
                // ωの補間次数 > ψの補間次数で頂点以外の節点
                System.Diagnostics.Debug.Assert(wLineFE.Order == 2);
                uint elemNodeCnt = wLineFE.NodeCount;
                uint elemVertexCnt = wLineFE.VertexCount;
                int[] coIds = wLineFE.NodeCoordIds;
                int nodeIndex = -1;
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    if (coIds[iNode] == wCoId)
                    {
                        nodeIndex = iNode;
                        break;
                    }
                }
                System.Diagnostics.Debug.Assert(nodeIndex != -1);
                System.Diagnostics.Debug.Assert(nodeIndex == 2);
                int coId1 = coIds[0];
                int coId2 = coIds[1];
                double[] v1 = { CoordVelocity[coId1 * vDof], CoordVelocity[coId1 * vDof + 1] };
                double[] v2 = { CoordVelocity[coId2 * vDof], CoordVelocity[coId2 * vDof + 1] };
                v[0] = (v1[0] + v2[0]) / 2.0;
                v[1] = (v1[1] + v2[1]) / 2.0;
            }
            return v;
        }

        private void SetVorticityDirichletBCOfStreamForOutflow(
            IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint wQuantityId = 0;
            uint pQuantityId = 1;
            int wNodeCnt = (int)World.GetNodeCount(wQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int offset = wNodeCnt;

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
                FlowVorticityBCType bcType = (FlowVorticityBCType)intParam[0];
                if (bcType != FlowVorticityBCType.Outflow)
                {
                    continue;
                }
                IList<uint> bcEIds = portCondition.EIds;
                //IList<double> param = portCondition.DoubleAdditionalParameters;
                //System.Diagnostics.Debug.Assert(param.Count == 0);

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
                    int[] wNodes = new int[elemNodeCnt];
                    int[] pNodes = new int[elemNodeCnt];
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int coId = pCoIds[iNode];
                        wNodes[iNode] = World.Coord2Node(wQuantityId, coId);
                        pNodes[iNode] = World.Coord2Node(pQuantityId, coId);
                    }

                    Material ma0 = World.GetMaterial(lineFE.MaterialId);
                    System.Diagnostics.Debug.Assert(ma0 is NewtonFluidMaterial);
                    var ma = ma0 as NewtonFluidMaterial;
                    double rho = ma.MassDensity;
                    double mu = ma.Mu;
                    double[] g = { ma.GravityX, ma.GravityY };

                    double[] normal = lineFE.GetNormal();
                    double[] tan = { normal[1], -normal[0] };

                    int coId1 = pCoIds[0];
                    int coId2 = pCoIds[1];
                    int adjCoId = -1; // adjacent
                    uint triFEId1;
                    {
                        IList<uint> triFEIds = World.GetTriangleFEIdsFromEdgeCoord(pQuantityId, coId1, coId2);
                        uint triFEId = triFEIds[0];
                        triFEId1 = triFEId;
                        TriangleFE triFE = World.GetTriangleFE(pQuantityId, triFEId);
                        int[] triFECoIds = triFE.NodeCoordIds;
                        foreach (int coId in triFECoIds)
                        {
                            if (coId != coId1 && coId != coId2)
                            {
                                adjCoId = coId;
                                break;
                            }
                        }
                        System.Diagnostics.Debug.Assert(adjCoId != -1);
                    }
                    int pAdjNodeId = -1;
                    pAdjNodeId = World.Coord2Node(pQuantityId, adjCoId); // 特殊な場合-1はありえる
                    double[] pt1 = World.GetCoord(pQuantityId, coId1);
                    double[] pt2 = World.GetCoord(pQuantityId, coId2);
                    double[] adjPt = World.GetCoord(pQuantityId, adjCoId);
                    double[] dir = IvyFEM.CadUtils2D.GetDirection2D(pt1, pt2);
                    int adjCoId2 = -1;
                    {
                        IList<uint> triFEIds = World.GetTriangleFEIdsFromCoord(pQuantityId, adjCoId);
                        double proj = 0.0;
                        foreach (uint triFEId in triFEIds)
                        {
                            if (triFEId == triFEId1)
                            {
                                continue;
                            }
                            TriangleFE triFE = World.GetTriangleFE(pQuantityId, triFEId);
                            int[] triFECoIds = triFE.NodeCoordIds;
                            foreach (int coId in triFECoIds)
                            {
                                if (coId == adjCoId)
                                {
                                    continue;
                                }
                                double[] tmpPt = World.GetCoord(pQuantityId, coId);
                                double[] tmpDir = IvyFEM.CadUtils2D.GetDirection2D(adjPt, tmpPt);
                                double tmpProj = dir[0] * tmpDir[0] + dir[1] * tmpDir[1];
                                //if (tmpProj > proj) // 同じ方向
                                if (Math.Abs(tmpProj) > Math.Abs(proj)) // projectionが大きければ方向は逆も許容
                                {
                                    proj = tmpProj;
                                    adjCoId2 = coId;
                                }
                            }
                        }
                        System.Diagnostics.Debug.Assert(adjCoId2 != -1);
                    }
                    int pAdjNodeId2 = -1;
                    pAdjNodeId2 = World.Coord2Node(pQuantityId, adjCoId2); // 特殊な場合-1はありえる
                    double[] adjPt2 = World.GetCoord(pQuantityId, adjCoId2);
                    double h2 = dir[0] * (adjPt2[0] - adjPt[0]) + dir[1] * (adjPt2[1] - adjPt[1]); // マイナスもある

                    for (int row = 0; row < elemNodeCnt; row++)
                    {
                        int rowCoId = pCoIds[row];
                        int rowNodeId = pNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        int rowCoId2 = pCoIds[(row + 1) % elemNodeCnt];
                        int rowNodeId2 = pNodes[(row + 1) % elemNodeCnt];

                        double[] rowPt1 = World.GetCoord(pQuantityId, rowCoId);
                        double[] rowPt2 = World.GetCoord(pQuantityId, rowCoId2);
                        double h1 = dir[0] * (rowPt2[0] - rowPt1[0]) + dir[1] * (rowPt2[1] - rowPt1[1]); // マイナスもある

                        for (int colNodeId = 0; colNodeId < wNodeCnt; colNodeId++)
                        {
                            A[offset + rowNodeId, colNodeId] = 0;
                        }
                        for (int colNodeId = 0; colNodeId < pNodeCnt; colNodeId++)
                        {
                            int colCoId = World.Node2Coord(pQuantityId, colNodeId);
                            if (colCoId == rowCoId)
                            {
                                A[offset + rowNodeId, offset + colNodeId] = 1.0 / h1;
                            }
                            else if (colCoId == rowCoId2)
                            {
                                A[offset + rowNodeId, offset + colNodeId] = -1.0 / h1;
                            }
                            else if (pAdjNodeId != -1 && colNodeId == pAdjNodeId)
                            {
                                A[offset + rowNodeId, offset + colNodeId] = -1.0 / h2;
                            }
                            else if (pAdjNodeId2 != -1 && colNodeId == pAdjNodeId2)
                            {
                                A[offset + rowNodeId, offset + colNodeId] = 1.0 / h2;
                            }
                            else
                            {
                                A[offset + rowNodeId, offset + colNodeId] = 0;
                            }
                        }
                        B[offset + rowNodeId] = 0;
                    }
                }
            }
        }

        private void GetVelocityFromVorticityStream()
        {
            uint wQuantityId = 0;
            uint pQuantityId = 1;
            int wNodeCnt = (int)World.GetNodeCount(wQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int offset = wNodeCnt;
            int wCoCnt = (int)World.GetCoordCount(wQuantityId);
            int pCoCnt = (int)World.GetCoordCount(pQuantityId);

            int vDof = 2; // 速度
            CoordVelocity = new double[pCoCnt * vDof]; // nodeでなくcoord
            int[] cnt = new int[pCoCnt];

            IList<uint> feIds = World.GetTriangleFEIds(pQuantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE pTriFE = World.GetTriangleFE(pQuantityId, feId);

                int[] pCoIds = pTriFE.NodeCoordIds;
                uint pElemNodeCnt = pTriFE.NodeCount;
                int[] pNodes = new int[pElemNodeCnt];
                for (int iNode = 0; iNode < pElemNodeCnt; iNode++)
                {
                    int coId = pCoIds[iNode];
                    int nodeId = World.Coord2Node(pQuantityId, coId);
                    pNodes[iNode] = nodeId;
                }

                Material ma0 = World.GetMaterial(pTriFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is NewtonFluidMaterial);
                var ma = ma0 as NewtonFluidMaterial;
                double rho = ma.MassDensity;
                double mu = ma.Mu;
                double[] g = { ma.GravityX, ma.GravityY };

                for (int row = 0; row < pElemNodeCnt; row++)
                {
                    int rowCoId = pCoIds[row];
                    double[] L = pTriFE.GetNodeL(row);
                    double[] pN = pTriFE.CalcN(L);
                    double[][] pNu = pTriFE.CalcNu(L);
                    double[] pNx = pNu[0];
                    double[] pNy = pNu[1];
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

                    for (int rowDof = 0; rowDof < vDof; rowDof++)
                    {
                        CoordVelocity[rowCoId * vDof + rowDof] += v[rowDof];
                    }
                    cnt[rowCoId]++;
                }
            }
            for (int coId = 0; coId < pCoCnt; coId++)
            {
                if (cnt[coId] > 0)
                {
                    for (int iDof = 0; iDof < vDof; iDof++)
                    {
                        CoordVelocity[coId * vDof + iDof] /= cnt[coId];
                    }
                }
            }
        }
    }
}
