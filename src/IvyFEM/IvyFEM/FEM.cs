using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    abstract public class FEM
    {
        public FEWorld World { get; set; } = null;
        public IvyFEM.Linear.IEquationSolver Solver { get; set; } = null;

        public abstract void Solve();

        protected static void SetFixedCadsCondtion(FEWorld world, IvyFEM.Linear.DoubleSparseMatrix A, double[]B,
            int[] nodeCnts, int[] dofs)
        {
            // A21の右辺移行
            // Note:速度改善のためcolを先にしている
            for (uint colQuantityId = 0; colQuantityId < dofs.Length; colQuantityId++)
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
                    // fixed節点、自由度
                    int colCoId = world.Node2Coord(colQuantityId, colNodeId);
                    var colFixedCads = world.GetFixedCadsFromCoord(colQuantityId, colCoId);
                    if (colFixedCads.Count == 0)
                    {
                        continue;
                    }
                    foreach (var colFixedCad in colFixedCads)
                    {
                        foreach (uint colFixedDofIndex in colFixedCad.FixedDofIndexs)
                        {
                            double value = colFixedCad.DoubleValues[colFixedDofIndex];
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
                                    int rowCoId = world.Node2Coord(rowQuantityId, rowNodeId);
                                    IList<FieldFixedCad> rowFixedCads = world.GetFixedCadsFromCoord(rowQuantityId, rowCoId);
                                    // fixedでない節点、自由度
                                    IList<uint> rowDofs = new List<uint>();
                                    for (uint rowDof = 0; rowDof < rowDofCnt; rowDof++)
                                    {
                                        rowDofs.Add(rowDof);
                                    }
                                    foreach (var rowfixedCad in rowFixedCads)
                                    {
                                        foreach (uint rowFixedDof in rowfixedCad.FixedDofIndexs)
                                        {
                                            rowDofs.Remove(rowFixedDof);
                                        }
                                    }
                                    if (rowDofs.Count == 0)
                                    {
                                        continue;
                                    }

                                    foreach (int rowDof in rowDofs)
                                    {
                                        double a = A[rowOffset + rowNodeId * rowDofCnt + rowDof,
                                            colOffset + colNodeId * colDofCnt + (int)colFixedDofIndex];
                                        B[rowOffset + rowNodeId * rowDofCnt + rowDof] -= a * value;
                                        A[rowOffset + rowNodeId * rowDofCnt + rowDof,
                                            colOffset + colNodeId * colDofCnt + (int)colFixedDofIndex] = 0;
                                    }
                                }
                            }
                        }
                    }
                }
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
                    int rowCoId = world.Node2Coord(rowQuantityId, rowNodeId);
                    var rowFixedCads = world.GetFixedCadsFromCoord(rowQuantityId, rowCoId);
                    if (rowFixedCads.Count == 0)
                    {
                        continue;
                    }
                    foreach (var rowFixedCad in rowFixedCads)
                    {
                        foreach (uint rowFixedDof in rowFixedCad.FixedDofIndexs)
                        {
                            double value = rowFixedCad.DoubleValues[rowFixedDof];
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
                                            colDof == rowFixedDof) ? 1 : 0;
                                        A[rowOffset + rowNodeId * rowDofCnt + (int)rowFixedDof,
                                            colOffset + colNodeId * colDofCnt + colDof] = a;
                                    }
                                }
                                B[rowOffset + rowNodeId * rowDofCnt + rowFixedDof] = value;
                            }
                        }
                    }
                }
            }
        }
    }
}
