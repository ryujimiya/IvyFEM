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
            uint quantityId, IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            int nodeCnt = (int)World.GetNodeCount(quantityId);
            int dof = (int)World.GetDof(quantityId);

            // A11
            {
                int rowNodeCnt = nodeCnt;
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
                                int colNodeCnt = nodeCnt;
                                int colDofCnt = dof;
                                int colOffset = 0;
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
