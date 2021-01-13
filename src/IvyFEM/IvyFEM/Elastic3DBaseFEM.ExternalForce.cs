using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class Elastic3DBaseFEM
    {
        protected void SetExternalForceSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            int quantityCnt = World.GetQuantityCount();
            for (uint quantityId = 0; quantityId < quantityCnt; quantityId++)
            {
                IList<PortCondition> portConditions = World.GetPortConditions(quantityId);
                if (World.GetPortCount(quantityId) == 0)
                {
                    continue;
                }
                CadElementType cadElemType = portConditions[0].CadElemType;
                if (cadElemType == CadElementType.Loop)
                {
                    SetExternalForceQuantitySpecialBC(quantityId, A, B);
                }
                else if (cadElemType == CadElementType.Edge)
                {
                    SetExternalForceQuantitySpecialBCForPlate(quantityId, A, B);
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
            }
        }

        private void SetExternalForceQuantitySpecialBC(
            uint quantityId, IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            if (World.GetPortCount(quantityId) == 0)
            {
                return;
            }
            int nodeCnt = (int)World.GetNodeCount(quantityId);
            int dof = (int)World.GetDof(quantityId);
            int offset = World.GetOffset(quantityId);

            uint portCnt = World.GetPortCount(quantityId);
            IList<PortCondition> portConditions = World.GetPortConditions(quantityId);
            IList<uint> feIds = World.GetTriangleFEIds(quantityId);
            for (int portId = 0; portId < portCnt; portId++)
            {
                PortCondition portCondition = portConditions[portId];
                IList<int> intParam = portCondition.IntAdditionalParameters;
                System.Diagnostics.Debug.Assert(intParam.Count == 1);
                ElasticBCType bcType = (ElasticBCType)intParam[0];
                if (bcType != ElasticBCType.ExternalForce)
                {
                    continue;
                }
                IList<uint> bcLIds = portCondition.LIds;

                foreach (uint feId in feIds)
                {
                    TriangleFE triFE = World.GetTriangleFE(quantityId, feId);
                    uint meshId = triFE.MeshId;
                    //int meshElemId = triFE.MeshElemId;
                    uint lId;
                    {
                        uint elemCount;
                        MeshType meshType;
                        int loc;
                        uint cadId;
                        World.Mesh.GetMeshInfo(meshId, out elemCount, out meshType, out loc, out cadId);
                        System.Diagnostics.Debug.Assert(meshType == MeshType.Tri);
                        lId = cadId;
                    }
                    if (bcLIds.Contains(lId))
                    {
                        // BCを適用する面
                    }
                    else
                    {
                        continue;
                    }

                    // BC
                    uint elemNodeCnt = triFE.NodeCount;
                    int[] coIds = triFE.NodeCoordIds;
                    int[] nodes = new int[elemNodeCnt];
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int coId = coIds[iNode];
                        nodes[iNode] = World.Coord2Node(quantityId, coId);
                    }

                    double[] sN = triFE.CalcSN();
                    for (int row = 0; row < elemNodeCnt; row++)
                    {
                        int rowCoId = coIds[row];
                        int rowNodeId = nodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        IList<double> param = portCondition.GetDoubleAdditionalParameters(rowCoId);
                        System.Diagnostics.Debug.Assert(param.Count == dof);
                        double[] f = param.ToArray();
                        for (int rowDof = 0; rowDof < dof; rowDof++)
                        {
                            B[offset + rowNodeId * dof + rowDof] += sN[row] * f[rowDof];
                        }
                    }
                }
            }
        }

        private void SetExternalForceQuantitySpecialBCForPlate(
            uint quantityId, IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            if (World.GetPortCount(quantityId) == 0)
            {
                return;
            }
            int nodeCnt = (int)World.GetNodeCount(quantityId);
            int dof = (int)World.GetDof(quantityId);
            int offset = World.GetOffset(quantityId);

            uint portCnt = World.GetPortCount(quantityId);
            IList<PortCondition> portConditions = World.GetPortConditions(quantityId);
            IList<uint> feIds = World.GetLineFEIds(quantityId);
            for (int portId = 0; portId < portCnt; portId++)
            {
                PortCondition portCondition = portConditions[portId];
                IList<int> intParam = portCondition.IntAdditionalParameters;
                System.Diagnostics.Debug.Assert(intParam.Count == 1);
                ElasticBCType bcType = (ElasticBCType)intParam[0];
                if (bcType != ElasticBCType.ExternalForce)
                {
                    continue;
                }
                IList<uint> bcEIds = portCondition.EIds;

                foreach (uint feId in feIds)
                {
                    LineFE lineFE = World.GetLineFE(quantityId, feId);
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
                    int[] coIds = lineFE.NodeCoordIds;
                    int[] nodes = new int[elemNodeCnt];
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int coId = coIds[iNode];
                        nodes[iNode] = World.Coord2Node(quantityId, coId);
                    }

                    double[] sN = lineFE.CalcSN();
                    for (int row = 0; row < elemNodeCnt; row++)
                    {
                        int rowCoId = coIds[row];
                        int rowNodeId = nodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        IList<double> param = portCondition.GetDoubleAdditionalParameters(rowCoId);
                        System.Diagnostics.Debug.Assert(param.Count == dof);
                        // param: DKTなら d1Quantityのとき 0: fx 1: fy, d2Quantityのとき 0: fz 
                        double[] f = param.ToArray();
                        for (int rowDof = 0; rowDof < dof; rowDof++)
                        {
                            B[offset + rowNodeId * dof + rowDof] += sN[row] * f[rowDof];
                        }
                    }
                }
            }
        }
    }
}
