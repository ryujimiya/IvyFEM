using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class Elastic2DBaseFEM
    {
        protected void SetExternalForceSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint quantityId = 0;
            SetExternalForceQuantitySpecialBC(quantityId, A, B);
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
            int offset = GetOffset(quantityId);

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
                        System.Diagnostics.Debug.Assert(param.Count == 2); // 0: fx 1: fy
                        System.Diagnostics.Debug.Assert(dof == 2);
                        double fxValue = param[0];
                        double fyValue = param[1];
                        double[] f = { fxValue, fyValue };
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
