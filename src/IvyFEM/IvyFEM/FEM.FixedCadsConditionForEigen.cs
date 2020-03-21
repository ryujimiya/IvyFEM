using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class FEM
    {
        protected void DoubleSetFixedCadsCondtionForEigen(IvyFEM.Lapack.DoubleMatrix A, IvyFEM.Lapack.DoubleMatrix B)
        {
            int quantityCnt = World.GetQuantityCount();
            int[] dofs = new int[quantityCnt];
            int[] nodeCnts = new int[quantityCnt];
            for (uint quantityId = 0; quantityId < quantityCnt; quantityId++)
            {
                dofs[quantityId] = (int)World.GetDof(quantityId);
                nodeCnts[quantityId] = (int)World.GetNodeCount(quantityId);
            }

            // A11, A12
            for (uint rowQuantityId = 0; rowQuantityId < dofs.Length; rowQuantityId++)
            {
                int rowNodeCnt = nodeCnts[rowQuantityId];
                int rowDofCnt = dofs[rowQuantityId];
                int rowOffset = 0;
                for (int i = 0; i < rowQuantityId; i++)
                {
                    rowOffset += nodeCnts[i] * dofs[i];
                }
                for (int rowNodeId = 0; rowNodeId < rowNodeCnt; rowNodeId++)
                {
                    int rowCoId = World.Node2Coord(rowQuantityId, rowNodeId);
                    var rowFixedCads = World.GetFixedCadsFromCoord(rowQuantityId, rowCoId);
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
                            // 強制境界に限る
                            System.Diagnostics.Debug.Assert(Math.Abs(value) < IvyFEM.Constants.PrecisionLowerLimit);
                            for (int colQuantityId = 0; colQuantityId < dofs.Length; colQuantityId++)
                            {
                                int colNodeCnt = nodeCnts[colQuantityId];
                                int colDofCnt = dofs[colQuantityId];
                                int colOffset = 0;
                                for (int i = 0; i < colQuantityId; i++)
                                {
                                    colOffset += nodeCnts[i] * dofs[i];
                                }
                                for (int colNodeId = 0; colNodeId < colNodeCnt; colNodeId++)
                                {
                                    for (int colDof = 0; colDof < colDofCnt; colDof++)
                                    {
                                        double a = (colQuantityId == rowQuantityId &&
                                            colNodeId == rowNodeId &&
                                            colDof == rowFixedDof) ? 1.0 : 0.0;
                                        A[rowOffset + rowNodeId * rowDofCnt + (int)rowFixedDof,
                                            colOffset + colNodeId * colDofCnt + colDof] = a;
                                        B[rowOffset + rowNodeId * rowDofCnt + (int)rowFixedDof,
                                            colOffset + colNodeId * colDofCnt + colDof] = 0.0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
