using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Fluid2DRKTDFEM : FEM
    {
        public double ConvRatioToleranceForNonlinearIter { get; set; }
            = 1.0e+2 * IvyFEM.Linear.Constants.ConvRatioTolerance; // 収束しないので収束条件を緩めている

        public double TimeStep { get; private set; } = 0;

        public uint ValueId { get; private set; } = 0;
        public uint PrevValueId { get; private set; } = 0;

        /// <summary>
        /// 方程式のタイプ
        /// </summary>
        public FluidEquationType EquationType { get; set; } = FluidEquationType.StdGPressurePoisson;

        // output
        public double[] U = null;
        public double[] CoordVelocity { get; protected set; } = null; // Vorticityの場合の速度を格納


        public Fluid2DRKTDFEM(FEWorld world, double timeStep, uint valueId, uint prevValueId)
        {
            World = world;
            TimeStep = timeStep;
            ValueId = valueId;
            PrevValueId = prevValueId;
        }

        private void Calc1stEquationKMF(
            IvyFEM.Linear.DoubleSparseMatrix K, IvyFEM.Linear.DoubleSparseMatrix M, double[] F)
        {
            switch (EquationType)
            {
                case FluidEquationType.StdGVorticity:
                    CalcStdGVorticity1stEquationKMF(K, M, F);
                    break;
                case FluidEquationType.StdGPressurePoisson:
                    CalcStdGPressurePoisson1stEquationKMF(K, M, F);
                    break;
                case FluidEquationType.StdGPressurePoissonWithBell:
                    CalcStdGPressurePoissonWithBell1stEquationKMF(K, M, F);
                    break;
                case FluidEquationType.SUPGPressurePoissonWithBell:
                    CalcSUPGPressurePoissonWithBell1stEquationKMF(K, M, F);
                    break;
                default:
                    System.Diagnostics.Debug.Assert(false);
                    break;
            }
        }

        private void Set1stEquationSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            if (EquationType == FluidEquationType.StdGVorticity)
            {
                SetVorticity1stEquationSpecialBC(A, B);
            }
            else if (EquationType == FluidEquationType.StdGPressurePoisson)
            {
                SetPressurePoisson1stEquationSpecialBC(A, B);
            }
            else if (EquationType == FluidEquationType.StdGPressurePoissonWithBell ||
                EquationType == FluidEquationType.SUPGPressurePoissonWithBell)
            {
                SetPressurePoissonWithBell1stEquationSpecialBC(A, B);
            }
        }

        private void Calc2ndEquationAB(
            IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            switch (EquationType)
            {
                case FluidEquationType.StdGVorticity:
                    CalcStdGVorticity2ndEquationAB(A, B);
                    break;
                case FluidEquationType.StdGPressurePoisson:
                    CalcStdGPressurePoisson2ndEquationAB(A, B);
                    break;
                case FluidEquationType.StdGPressurePoissonWithBell:
                    CalcStdGPressurePoissonWithBell2ndEquationAB(A, B);
                    break;
                case FluidEquationType.SUPGPressurePoissonWithBell:
                    CalcSUPGPressurePoissonWithBell2ndEquationAB(A, B);
                    break;
                default:
                    System.Diagnostics.Debug.Assert(false);
                    break;
            }
        }

        private void Set2ndEquationSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            if (EquationType == FluidEquationType.StdGVorticity)
            {
                SetVorticity2ndEquationSpecialBC(A, B);
            }
            else if (EquationType == FluidEquationType.StdGPressurePoisson)
            {
                SetPressurePoisson2ndEquationSpecialBC(A, B);
            }
            else if (EquationType == FluidEquationType.StdGPressurePoissonWithBell ||
                EquationType == FluidEquationType.SUPGPressurePoissonWithBell)
            {
                SetPressurePoissonWithBell2ndEquationSpecialBC(A, B);
            }
        }

        private void PostSolve()
        {
            if (EquationType == FluidEquationType.StdGVorticity)
            {
                VorticityPostSolve();
            }
        }

        private int GetOffset(uint quantityId)
        {
            int offset = 0;
            int quantityCnt = World.GetQuantityCount();
            for (uint tmpQuantityId = 0; tmpQuantityId < quantityId; tmpQuantityId++)
            {
                int nodeCnt = (int)World.GetNodeCount(tmpQuantityId);
                int dof = (int)World.GetDof(tmpQuantityId);
                offset += nodeCnt * dof;
            }
            return offset;
        }

        private void Solve1stEquation()
        {
            double dt = TimeStep;
            var FV = World.GetFieldValue(ValueId);

            uint quantityId1 = 0;
            int quantityDof1 = (int)World.GetDof(quantityId1);
            int quantityNodeCnt1 = (int)World.GetNodeCount(quantityId1);
            int nodeCnt = quantityNodeCnt1 * quantityDof1;
            int offset = GetOffset(quantityId1);

            // 前の時刻の値をセット
            double[] prevU1 = new double[nodeCnt];
            for (int nodeId = 0; nodeId < quantityNodeCnt1; nodeId++)
            {
                int coId = World.Node2Coord(quantityId1, nodeId);
                double[] value = FV.GetDoubleValue(coId, FieldDerivativeType.Value);
                System.Diagnostics.Debug.Assert(value.Length == quantityDof1);
                for (int iDof = 0; iDof < quantityDof1; iDof++)
                {
                    prevU1[nodeId * quantityDof1 + iDof] = value[iDof];
                }
            }

            var K = new IvyFEM.Linear.DoubleSparseMatrix(nodeCnt, nodeCnt);
            var M = new IvyFEM.Linear.DoubleSparseMatrix(nodeCnt, nodeCnt);
            var F = new double[nodeCnt];

            Calc1stEquationKMF(K, M, F);

            // RK4
            int rkStepCnt = 4;
            double[][] rkkn = new double[rkStepCnt][];
            for (int iStep = 0; iStep < rkStepCnt; iStep++)
            {
                var A = new IvyFEM.Linear.DoubleSparseMatrix(M);
                double[] B = new double[nodeCnt];

                double[] tmpU = new double[nodeCnt];
                prevU1.CopyTo(tmpU, 0);
                if (iStep == 0)
                {
                    // nothing
                }
                else if (iStep == 1)
                {
                    for (int i = 0; i < nodeCnt; i++)
                    {
                        tmpU[i] += (1.0 / 2.0) * rkkn[0][i];
                    }
                }
                else if (iStep == 2)
                {
                    for (int i = 0; i < nodeCnt; i++)
                    {
                        tmpU[i] += (1.0 / 2.0) * rkkn[1][i];
                    }
                }
                else if (iStep == 3)
                {
                    for (int i = 0; i < nodeCnt; i++)
                    {
                        tmpU[i] += rkkn[2][i];
                    }
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }

                double[] KU = K * tmpU;
                double[] prevMU = M * prevU1;
                for (int i = 0; i < nodeCnt; i++)
                {
                    B[i] = dt * (-KU[i] + F[i]) + prevMU[i];
                }

                for (int i = 0; i < nodeCnt; i++)
                {
                    for (int j = 0; j < nodeCnt; j++)
                    {
                        A[i, j] *= (1.0 / dt);
                    }
                    B[i] *= (1.0 / dt);
                }
                Set1stEquationSpecialBC(A, B);
                DoubleSetSplitQuantityFixedCadsCondtion(quantityId1, quantityId1, A, B);

                double[] X;
                Solver.DoubleSolve(out X, A, B);
                rkkn[iStep] = new double[nodeCnt];
                for (int i = 0; i < nodeCnt; i++)
                {
                    rkkn[iStep][i] = X[i] - prevU1[i];
                }
            }
            double[] U1 = new double[nodeCnt];
            for (int i = 0; i < nodeCnt; i++)
            {
                U1[i] = prevU1[i] + (1.0 / 6.0) * (rkkn[0][i] + 2.0 * rkkn[1][i] + 2.0 * rkkn[2][i] + rkkn[3][i]);
            }

            U1.CopyTo(U, offset);
        }

        private void Solve2ndEquation(
            int iter, double tolerance, out bool isConverged, ref double sqInvNorm0, out double convRatio)
        {
            int quantityCnt = World.GetQuantityCount();

            isConverged = false;
            convRatio = 0;

            uint vQuantityId = 0;
            uint pQuantityId = 1;
            int nodeCnt = 0;
            for (uint quantityId2 = pQuantityId; quantityId2 < quantityCnt; quantityId2++)
            {
                int quantityDof2 = (int)World.GetDof(quantityId2);
                int quantityNodeCnt2 = (int)World.GetNodeCount(quantityId2);
                nodeCnt += quantityNodeCnt2 * quantityDof2;
            }
            int pOffset = GetOffset(pQuantityId);

            double[] U2 = new double[nodeCnt];

            for (uint quantityId2 = pQuantityId; quantityId2 < quantityCnt; quantityId2++)
            {
                int quantityDof2 = (int)World.GetDof(quantityId2);
                int quantityNodeCnt2 = (int)World.GetNodeCount(quantityId2);
                int offset = GetOffset(quantityId2);
                int offsetp = offset - pOffset;

                // 前の反復のときの値をセット(方程式を解くときに必要はないが、解く前の収束判定で使用する）
                for (int nodeId = 0; nodeId < quantityNodeCnt2; nodeId++)
                {
                    for (int iDof = 0; iDof < quantityDof2; iDof++)
                    {
                        U2[offsetp + nodeId * quantityDof2 + iDof] = U[offset + nodeId * quantityDof2 + iDof];
                    }
                }
            }

            var A = new IvyFEM.Linear.DoubleSparseMatrix(nodeCnt, nodeCnt);
            var B = new double[nodeCnt];

            Calc2ndEquationAB(A, B);
            Set2ndEquationSpecialBC(A, B);
            DoubleSetSplitQuantityFixedCadsCondtion(pQuantityId, (uint)(quantityCnt - 1), A, B);

            double[] AU = A * U2;
            double[] R = IvyFEM.Lapack.Functions.daxpy(-1.0, B, AU);
            double sqNorm = IvyFEM.Lapack.Functions.ddot(R, R);
            if (iter == 0)
            {
                if (sqNorm < IvyFEM.Constants.PrecisionLowerLimit)
                {
                    convRatio = 0;
                    isConverged = true;
                    return;
                }
                sqInvNorm0 = 1.0 / sqNorm;
            }
            else
            {
                convRatio = Math.Sqrt(sqNorm * sqInvNorm0);
                System.Diagnostics.Debug.WriteLine("cur convRatio =" + convRatio);
                if (sqNorm * sqInvNorm0 < tolerance * tolerance)
                {
                    isConverged = true;
                    return;
                }
            }

            double[] X;
            Solver.DoubleSolve(out X, A, B);
            U2 = X;

            U2.CopyTo(U, pOffset);
        }

        public override void Solve()
        {
            int quantityCnt = World.GetQuantityCount();
            //System.Diagnostics.Debug.Assert(quantityCnt == 2);
            // Bell elementの場合2ではない

            // Nonlinear Iter
            bool isConverged = false;
            double sqInvNorm0 = 0;
            double convRatio = ConvRatioToleranceForNonlinearIter;
            double tolerance = convRatio;
            const int maxIter = IvyFEM.Linear.Constants.MaxIter;
            int iter = 0;

            int nodeCnt = 0;
            for (uint quantityId = 0; quantityId < quantityCnt; quantityId++)
            {
                int quantityDof = (int)World.GetDof(quantityId);
                int quantityNodeCnt = (int)World.GetNodeCount(quantityId);
                nodeCnt += quantityNodeCnt * quantityDof;
            }
            U = new double[nodeCnt];

            PostSolve();
            for (iter = 0; iter < maxIter; iter++)
            {
                // 1st Equation
                Solve1stEquation();

                // 2nd Equation
                Solve2ndEquation(iter, tolerance, out isConverged, ref sqInvNorm0, out convRatio);

                PostSolve();
                if (isConverged)
                {
                    break;
                }
            }
            System.Diagnostics.Debug.WriteLine("Nonlinear iter = " + iter + " norm = " + convRatio);
            System.Diagnostics.Debug.Assert(iter < maxIter);
        }

        public void UpdateFieldValuesTimeDomain()
        {
            var FV = World.GetFieldValue(ValueId);
            var prevFV = World.GetFieldValue(PrevValueId);
            prevFV.Copy(FV);

            World.UpdateFieldValueValuesFromNodeValues(ValueId, FieldDerivativeType.Value, U);
        }
    }
}
