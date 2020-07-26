using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class FEM
    {
        protected void DoubleSetFixedCadsCondtionForEigen(
            IvyFEM.Lapack.DoubleMatrix A, IvyFEM.Lapack.DoubleMatrix B,
            uint maxQuantityId, bool isDoubleSize, int portId)
        {
            // ポートの固有値問題?
            bool isPort = portId != -1;
 
            int quantityCnt = World.GetQuantityCount();
            int[] dofs = new int[quantityCnt];
            int[] nodeCnts = new int[quantityCnt];
            for (uint quantityId = 0; quantityId < quantityCnt; quantityId++)
            {
                dofs[quantityId] = (int)World.GetDof(quantityId);
                if (isPort)
                {
                    nodeCnts[quantityId] = (int)World.GetPortNodeCount(quantityId, (uint)portId);
                }
                else
                {
                    nodeCnts[quantityId] = (int)World.GetNodeCount(quantityId);
                }
            }
            // for DoubleSize
            int halfOffset = -1;
            if (isDoubleSize)
            {
                int halfRowLen = A.RowLength / 2;
                halfOffset = halfRowLen;
            }

            // A11, A12
            for (uint rowQuantityId = 0; rowQuantityId < dofs.Length; rowQuantityId++)
            {
                if (rowQuantityId > maxQuantityId)
                {
                    continue;
                }
                int rowNodeCnt = nodeCnts[rowQuantityId];
                int rowDofCnt = dofs[rowQuantityId];
                int rowOffset = 0;
                for (int i = 0; i < rowQuantityId; i++)
                {
                    rowOffset += nodeCnts[i] * dofs[i];
                }
                for (int rowNodeId = 0; rowNodeId < rowNodeCnt; rowNodeId++)
                {
                    int rowCoId;
                    if (isPort)
                    {
                        rowCoId = World.PortNode2Coord(rowQuantityId, (uint)portId, rowNodeId);
                    }
                    else
                    {
                        rowCoId = World.Node2Coord(rowQuantityId, rowNodeId);
                    }
                    var rowFixedCads = World.GetFixedCadsFromCoord(rowQuantityId, rowCoId);
                    if (rowFixedCads.Count == 0)
                    {
                        continue;
                    }
                    foreach (var rowFixedCad in rowFixedCads)
                    {
                        foreach (uint rowFixedDof in rowFixedCad.FixedDofIndexs)
                        {
                            // checkしない
                            //double[] values = rowFixedCad.GetDoubleValues(rowCoId);
                            //double value = values[rowFixedDof];
                            //// 強制境界に限る
                            //System.Diagnostics.Debug.Assert(Math.Abs(value) < IvyFEM.Constants.PrecisionLowerLimit);
                            for (int colQuantityId = 0; colQuantityId < dofs.Length; colQuantityId++)
                            {
                                if (colQuantityId > maxQuantityId)
                                {
                                    continue;
                                }
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
                                        if (isDoubleSize)
                                        {
                                            // 非線形の線形化でサイズが2倍のもの
                                            A[rowOffset + rowNodeId * rowDofCnt + (int)rowFixedDof,
                                                colOffset + colNodeId * colDofCnt + colDof] = a;
                                            A[rowOffset + rowNodeId * rowDofCnt + (int)rowFixedDof,
                                                colOffset + colNodeId * colDofCnt + colDof + halfOffset] = 0.0;
                                            A[rowOffset + rowNodeId * rowDofCnt + (int)rowFixedDof + halfOffset,
                                                colOffset + colNodeId * colDofCnt + colDof] = 0.0;
                                            A[rowOffset + rowNodeId * rowDofCnt + (int)rowFixedDof + halfOffset,
                                                colOffset + colNodeId * colDofCnt + colDof + halfOffset] = a;
                                            if (B != null)
                                            {
                                                B[rowOffset + rowNodeId * rowDofCnt + (int)rowFixedDof,
                                                    colOffset + colNodeId * colDofCnt + colDof] = 0.0;
                                                B[rowOffset + rowNodeId * rowDofCnt + (int)rowFixedDof,
                                                    colOffset + colNodeId * colDofCnt + colDof + halfOffset] = 0.0;
                                                B[rowOffset + rowNodeId * rowDofCnt + (int)rowFixedDof + halfOffset,
                                                    colOffset + colNodeId * colDofCnt + colDof] = 0.0;
                                                B[rowOffset + rowNodeId * rowDofCnt + (int)rowFixedDof + halfOffset,
                                                    colOffset + colNodeId * colDofCnt + colDof + halfOffset] = 0.0;
                                            }
                                        }
                                        else
                                        {
                                            // 通常
                                            A[rowOffset + rowNodeId * rowDofCnt + (int)rowFixedDof,
                                                colOffset + colNodeId * colDofCnt + colDof] = a;
                                            if (B != null)
                                            {
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
    }
}
