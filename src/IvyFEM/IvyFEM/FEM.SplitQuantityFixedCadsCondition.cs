using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class FEM
    {
        protected void DoubleSetSplitQuantityFixedCadsCondtion(
            uint sQuantityId, uint eQuantityId, IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            int nodeCnt = 0;
            int dof = (int)World.GetDof(sQuantityId);
            int[] offsets = new int[eQuantityId - sQuantityId + 1];
            for (uint quantityId = sQuantityId; quantityId <= eQuantityId; quantityId++)
            {
                offsets[quantityId - sQuantityId] = nodeCnt;
                nodeCnt += (int)World.GetNodeCount(quantityId) * dof;
                System.Diagnostics.Debug.Assert(World.GetDof(quantityId) == dof);
            }

            // A11
            for (uint quantityId = sQuantityId; quantityId <= eQuantityId; quantityId++)
            {
                int quantityNodeCnt = (int)World.GetNodeCount(quantityId);
                int rowNodeCnt = quantityNodeCnt;
                int rowDofCnt = dof;
                for (int rowNodeId = 0; rowNodeId < rowNodeCnt; rowNodeId++)
                {
                    int rowCoId = World.Node2Coord(quantityId, rowNodeId);
                    var rowFixedCads = World.GetFixedCadsFromCoord(quantityId, rowCoId);
                    if (rowFixedCads.Count == 0)
                    {
                        continue;
                    }
                    foreach (var rowFixedCad in rowFixedCads)
                    {
                        foreach (uint rowFixedDof in rowFixedCad.FixedDofIndexs)
                        {
                            double[] values = rowFixedCad.GetDoubleValues(rowCoId);
                            double value = values[rowFixedDof];
                            {
                                int colNodeCnt = quantityNodeCnt;
                                int colDofCnt = dof;
                                int colOffset = offsets[quantityId];
                                for (int colNodeId = 0; colNodeId < colNodeCnt; colNodeId++)
                                {
                                    for (int colDof = 0; colDof < colDofCnt; colDof++)
                                    {
                                        double a = (colNodeId == rowNodeId &&
                                            colDof == rowFixedDof) ? 1 : 0;
                                        A[rowNodeId * rowDofCnt + (int)rowFixedDof,
                                            colOffset + colNodeId * colDofCnt + colDof] = a;
                                    }
                                }
                                B[rowNodeId * rowDofCnt + rowFixedDof] = value;
                            }
                        }
                    }
                }
            }
        }
    }
}
