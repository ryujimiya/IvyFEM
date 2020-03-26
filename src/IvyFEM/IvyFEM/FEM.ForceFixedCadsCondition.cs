using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class FEM
    {
        protected void DoubleSetForceFixedCadsCondtion(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            int quantityCnt = World.GetQuantityCount();
            int[] dofs = new int[quantityCnt];
            int[] nodeCnts = new int[quantityCnt];
            for (uint quantityId = 0; quantityId < quantityCnt; quantityId++)
            {
                dofs[quantityId] = (int)World.GetDof(quantityId);
                nodeCnts[quantityId] = (int)World.GetNodeCount(quantityId);
            }

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
                    var rowFixedCads = World.GetForceFixedCadsFromCoord(rowQuantityId, rowCoId);
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
                            B[rowOffset + rowNodeId * rowDofCnt + rowFixedDof] += value;
                        }
                    }
                }
            }
        }

        protected void ComplexSetForceFixedCadsCondtion(
            IvyFEM.Linear.ComplexSparseMatrix A, System.Numerics.Complex[] B,
            int[] nodeCnts, int[] dofs)
        {
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
                    var rowFixedCads = World.GetForceFixedCadsFromCoord(rowQuantityId, rowCoId);
                    if (rowFixedCads.Count == 0)
                    {
                        continue;
                    }
                    foreach (var rowFixedCad in rowFixedCads)
                    {
                        foreach (uint rowFixedDof in rowFixedCad.FixedDofIndexs)
                        {
                            System.Numerics.Complex[] values = rowFixedCad.GetComplexValues(rowCoId);
                            System.Numerics.Complex value = values[rowFixedDof];
                            B[rowOffset + rowNodeId * rowDofCnt + rowFixedDof] += value;
                        }
                    }
                }
            }
        }
    }
}
