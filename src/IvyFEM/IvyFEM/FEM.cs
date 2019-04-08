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

        protected void DoubleSetFixedCadsCondtion(IvyFEM.Linear.DoubleSparseMatrix A, double[] B,
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
                    int colCoId = World.Node2Coord(colQuantityId, colNodeId);
                    var colFixedCads = World.GetFixedCadsFromCoord(colQuantityId, colCoId);
                    if (colFixedCads.Count == 0)
                    {
                        continue;
                    }
                    foreach (var colFixedCad in colFixedCads)
                    {
                        foreach (uint colFixedDofIndex in colFixedCad.FixedDofIndexs)
                        {
                            double[] values = colFixedCad.GetDoubleValues(colCoId);
                            double value = values[colFixedDofIndex];
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
                                    IList<FieldFixedCad> rowFixedCads =
                                        World.GetFixedCadsFromCoord(rowQuantityId, rowCoId);
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

        protected void ComplexSetFixedCadsCondtion(
            IvyFEM.Linear.ComplexSparseMatrix A, System.Numerics.Complex[] B,
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
                    int colCoId = World.Node2Coord(colQuantityId, colNodeId);
                    var colFixedCads = World.GetFixedCadsFromCoord(colQuantityId, colCoId);
                    if (colFixedCads.Count == 0)
                    {
                        continue;
                    }
                    foreach (var colFixedCad in colFixedCads)
                    {
                        foreach (uint colFixedDofIndex in colFixedCad.FixedDofIndexs)
                        {
                            System.Numerics.Complex[] values = colFixedCad.GetComplexValues(colCoId);
                            System.Numerics.Complex value = values[colFixedDofIndex];
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
                                    IList<FieldFixedCad> rowFixedCads =
                                        World.GetFixedCadsFromCoord(rowQuantityId, rowCoId);
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
                                        System.Numerics.Complex a = A[rowOffset + rowNodeId * rowDofCnt + rowDof,
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
                            System.Numerics.Complex[] values = rowFixedCad.GetComplexValues(rowCoId);
                            System.Numerics.Complex value = values[rowFixedDof];
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
                                        System.Numerics.Complex a = (colQuantityId == rowQuantityId &&
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

        protected void UpdateFieldValuesTimeDomain(
            double[] U,
            uint valueId, uint prevValueId,
            double timeStep,
            double newmarkBeta, double newmarkGamma)
        {
            double dt = timeStep;
            double beta = newmarkBeta;
            double gamma = newmarkGamma;

            var FV = World.GetFieldValue(valueId);
            var prevFV = World.GetFieldValue(prevValueId);
            prevFV.Copy(FV);

            World.UpdateFieldValueValuesFromNodeValues(valueId, FieldDerivativeType.Value, U);

            double[] u = FV.GetDoubleValues(FieldDerivativeType.Value);
            double[] vel = FV.GetDoubleValues(FieldDerivativeType.Velocity);
            double[] acc = FV.GetDoubleValues(FieldDerivativeType.Acceleration);
            double[] prevU = prevFV.GetDoubleValues(FieldDerivativeType.Value);
            double[] prevVel = prevFV.GetDoubleValues(FieldDerivativeType.Velocity);
            double[] prevAcc = prevFV.GetDoubleValues(FieldDerivativeType.Acceleration);

            uint coCnt = FV.GetPointCount();
            uint quantityId = FV.QuantityId;
            int dof = (int)FV.Dof;
            System.Diagnostics.Debug.Assert(coCnt == World.GetCoordCount(quantityId));
            System.Diagnostics.Debug.Assert(u.Length == coCnt * dof);
            for (int iPt = 0; iPt < coCnt; iPt++)
            {
                for (int iDof = 0; iDof < dof; iDof++)
                {
                    int index = iPt * dof + iDof;
                    vel[index] =
                        (gamma / (beta * dt)) * (u[index] - prevU[index]) +
                        (1.0 - gamma / beta) * prevVel[index] +
                        dt * (1.0 - gamma / (2.0 * beta)) * prevAcc[index];
                    acc[index] =
                        (1.0 / (beta * dt * dt)) * (u[index] - prevU[index]) -
                        (1.0 / (beta * dt)) * prevVel[index] -
                        (1.0 / (2.0 * beta) - 1.0) * prevAcc[index];
                }
            }
        }
    }
}
