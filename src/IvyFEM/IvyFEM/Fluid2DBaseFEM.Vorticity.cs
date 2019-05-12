using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class Fluid2DBaseFEM
    {
        // 境界が流線方向と同じとき
        // ψの法線成分（速度の接線成分に比例）がある
        // ψ = const or 0
        // Taylor級数展開 2次までの項を考慮
        protected void SetVorticityDirichletBCOfVorticityForTangentFlow(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
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
                //IList<int> intParam = portCondition.IntAdditionalParameters;
                //System.Diagnostics.Debug.Assert(intParam.Count == 0);
                IList<uint> bcEIds = portCondition.EIds;
                IList<double> param = portCondition.DoubleAdditionalParameters;
                System.Diagnostics.Debug.Assert(param.Count == 2); // 0: dψ/dx 1: dψ/dy
                double pxValue = param[0];
                double pyValue = param[1];

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
                    double hn = IvyFEM.CadUtils.TriHeight(
                        new OpenTK.Vector2d(adjPt[0], adjPt[1]), // 点
                        new OpenTK.Vector2d(pt1[0], pt1[1]), // 辺の点1
                        new OpenTK.Vector2d(pt2[0], pt2[1]) // 辺の点2
                        );

                    for (int row = 0; row < elemNodeCnt; row++)
                    {
                        int rowCoId = wCoIds[row];
                        int rowNodeId = wNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
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
        /*
        // Taylor級数展開 3次までの項を考慮
        protected void SetVorticityDirichletBCOfVorticityForTangentFlow(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
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
                //IList<int> intParam = portCondition.IntAdditionalParameters;
                //System.Diagnostics.Debug.Assert(intParam.Count == 0);
                IList<uint> bcEIds = portCondition.EIds;
                IList<double> param = portCondition.DoubleAdditionalParameters;
                System.Diagnostics.Debug.Assert(param.Count == 2); // 0: dψ/dx 1: dψ/dy
                double pxValue = param[0];
                double pyValue = param[1];

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
                    int wAdjNodeId = -1;
                    wAdjNodeId = World.Coord2Node(wQuantityId, adjCoId); // 特殊な場合-1はありえる
                    int pAdjNodeId = -1;
                    pAdjNodeId = World.Coord2Node(pQuantityId, adjCoId); // 特殊な場合-1はありえる
                    double[] pt1 = World.GetCoord(wQuantityId, coId1);
                    double[] pt2 = World.GetCoord(wQuantityId, coId2);
                    double[] adjPt = World.GetCoord(wQuantityId, adjCoId);
                    double hn = IvyFEM.CadUtils.TriHeight(
                        new OpenTK.Vector2d(adjPt[0], adjPt[1]), // 点
                        new OpenTK.Vector2d(pt1[0], pt1[1]), // 辺の点1
                        new OpenTK.Vector2d(pt2[0], pt2[1]) // 辺の点2
                        );

                    for (int row = 0; row < elemNodeCnt; row++)
                    {
                        int rowCoId = wCoIds[row];
                        int rowNodeId = wNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int colNodeId = 0; colNodeId < wNodeCnt; colNodeId++)
                        {
                            int colCoId = World.Node2Coord(wQuantityId, colNodeId);
                            if (colCoId == rowCoId)
                            {
                                A[rowNodeId, colNodeId] = 1.0;
                            }
                            else if (wAdjNodeId != -1 && colNodeId == wAdjNodeId)
                            {
                                A[rowNodeId, colNodeId] = -1.0 / 2.0;
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
                                A[rowNodeId, offset + colNodeId] = -3.0 / (hn * hn);
                            }
                            else if (pAdjNodeId != -1 && colNodeId == pAdjNodeId)
                            {
                                A[rowNodeId, offset + colNodeId] = 3.0 / (hn * hn);
                            }
                            else
                            {
                                A[rowNodeId, offset + colNodeId] = 0;
                            }
                        }
                        B[rowNodeId] = -(normal[0] * pxValue + normal[1] * pyValue) * (3.0 / hn);
                    }
                }
            }
        }
        */
        /*
        // Taylor級数展開 3次までの項を考慮 境界積分評価
        protected void SetVorticityDirichletBCOfVorticityForTangentFlow(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
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
                //IList<int> intParam = portCondition.IntAdditionalParameters;
                //System.Diagnostics.Debug.Assert(intParam.Count == 0);
                IList<uint> bcEIds = portCondition.EIds;
                IList<double> param = portCondition.DoubleAdditionalParameters;
                System.Diagnostics.Debug.Assert(param.Count == 2); // 0: dψ/dx 1: dψ/dy
                double pxValue = param[0];
                double pyValue = param[1];

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
                    double[] sN = lineFE.CalcSN();

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
                    int wAdjNodeId = -1;
                    wAdjNodeId = World.Coord2Node(wQuantityId, adjCoId); // 特殊な場合-1はありえる
                    int pAdjNodeId = -1;
                    pAdjNodeId = World.Coord2Node(pQuantityId, adjCoId); // 特殊な場合-1はありえる
                    double[] pt1 = World.GetCoord(wQuantityId, coId1);
                    double[] pt2 = World.GetCoord(wQuantityId, coId2);
                    double[] adjPt = World.GetCoord(wQuantityId, adjCoId);
                    double hn = IvyFEM.CadUtils.TriHeight(
                        new OpenTK.Vector2d(adjPt[0], adjPt[1]), // 点
                        new OpenTK.Vector2d(pt1[0], pt1[1]), // 辺の点1
                        new OpenTK.Vector2d(pt2[0], pt2[1]) // 辺の点2
                        );

                    for (int row = 0; row < elemNodeCnt; row++)
                    {
                        int rowCoId = wCoIds[row];
                        int rowNodeId = wNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int colNodeId = 0; colNodeId < wNodeCnt; colNodeId++)
                        {
                            int colCoId = World.Node2Coord(wQuantityId, colNodeId);
                            if (colCoId == rowCoId)
                            {
                                A[rowNodeId, colNodeId] += -mu * (3.0 / hn) * sN[row];
                            }
                        }
                        for (int colNodeId = 0; colNodeId < pNodeCnt; colNodeId++)
                        {
                            int colCoId = World.Node2Coord(pQuantityId, colNodeId);
                            if (colCoId == rowCoId)
                            {
                                A[rowNodeId, offset + colNodeId] += mu * (6.0 / (hn * hn * hn)) * sN[row];
                            }
                            else if (pAdjNodeId != -1 && colNodeId == pAdjNodeId)
                            {
                                A[rowNodeId, offset + colNodeId] += -mu * (6.0 / (hn * hn * hn)) * sN[row];
                            }
                        }
                        B[rowNodeId] += mu * (normal[0] * pxValue + normal[1] * pyValue) * (6.0 / (hn * hn)) * sN[row];
                    }
                }
            }
        }
        */

        protected void GetVelocityFromVorticityStream()
        {
            uint wQuantityId = 0;
            uint pQuantityId = 1;
            int wNodeCnt = (int)World.GetNodeCount(wQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int offset = wNodeCnt;
            int wCoCnt = (int)World.GetCoordCount(wQuantityId);
            int pCoCnt = (int)World.GetCoordCount(pQuantityId);

            int vDof = 2;
            CoordV = new double[pCoCnt * vDof]; // nodeでなくcoord

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
                        CoordV[rowCoId * vDof + rowDof] = v[rowDof];
                    }
                }
            }
        }
    }
}
