using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Fluid2DTDFEM : FEM
    {
        public double ConvRatioToleranceForNewtonRaphson { get; set; }
            = 1.0e+2 * IvyFEM.Linear.Constants.ConvRatioTolerance; // 少し収束条件を緩くしている
        public double TimeStep { get; private set; } = 0;
        public double NewmarkBeta { get; private set; } = 1.0 / 4.0;
        public double NewmarkGamma { get; private set; } = 1.0 / 2.0;
        public uint ValueId { get; private set; } = 0;
        public uint PrevValueId { get; private set; } = 0;

        /// <summary>
        /// 方程式のタイプ
        /// </summary>
        public FluidEquationType EquationType { get; set; } = FluidEquationType.StdGNavierStokes;

        // output
        public double[] U { get; private set; } = null;

        public Fluid2DTDFEM(FEWorld world,
            double timeStep,
            double newmarkBeta, double newmarkGamma,
            uint valueId, uint prevValueId)
        {
            World = world;
            TimeStep = timeStep;
            NewmarkBeta = newmarkBeta;
            NewmarkGamma = newmarkGamma;
            ValueId = valueId;
            PrevValueId = prevValueId;
        }

        public void UpdateFieldValuesTimeDomain()
        {
            UpdateFieldValuesTimeDomain(
                U, ValueId, PrevValueId,
                TimeStep,
                NewmarkBeta, NewmarkGamma);
        }

        private void CalcAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            switch (EquationType)
            {
                case FluidEquationType.Stokes:
                    CalcStokesAB(A, B);
                    break;
                case FluidEquationType.StdGNavierStokes:
                    CalcStdGNavierStokesAB(A, B);
                    break;
                case FluidEquationType.SUPGNavierStokes:
                    //CalcSUPGNavierStokesAB(A, B);
                    CalcSUPGNavierStokesByPicardAB(A, B);
                    break;
                default:
                    System.Diagnostics.Debug.Assert(false);
                    break;
            }
        }

        private bool MustUseNewtonRaphson()
        {
            if (EquationType == FluidEquationType.Stokes)
            {
                return false;
            }
            return true;
        }

        public override void Solve()
        {
            uint vQuantityId = 0;
            int vDof = 2;
            uint pQuantityId = 1;
            int pDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int nodeCnt = vNodeCnt * vDof + pNodeCnt * pDof;

            if (MustUseNewtonRaphson())
            {
                // Newton Raphson
                double sqNorm = 0;
                double sqInvNorm0 = 0;
                double convRatio = ConvRatioToleranceForNewtonRaphson;
                double tolerance = convRatio;
                const int maxIter = IvyFEM.Linear.Constants.MaxIter;
                int iter = 0;

                U = new double[nodeCnt];

                for (iter = 0; iter < maxIter; iter++)
                {
                    var A = new IvyFEM.Linear.DoubleSparseMatrix(nodeCnt, nodeCnt);
                    var B = new double[nodeCnt];

                    CalcAB(A, B);

                    DoubleSetFixedCadsCondtion(A, B, new int[] { vNodeCnt, pNodeCnt }, new int[] { vDof, pDof });

                    double[] AU = A * U;
                    double[] R = IvyFEM.Lapack.Functions.daxpy(-1.0, B, AU);
                    sqNorm = IvyFEM.Lapack.Functions.ddot(R, R);
                    if (iter == 0)
                    {
                        if (sqNorm < IvyFEM.Constants.PrecisionLowerLimit)
                        {
                            convRatio = 0;
                            break;
                        }
                        sqInvNorm0 = 1.0 / sqNorm;
                    }
                    else
                    {
                        convRatio = Math.Sqrt(sqNorm * sqInvNorm0);
                        System.Diagnostics.Debug.WriteLine("cur convRatio =" + convRatio);
                        if (sqNorm * sqInvNorm0 < tolerance * tolerance)
                        {
                            break;
                        }
                    }

                    //---------------------------------------------------
                    double[] X;
                    Solver.DoubleSolve(out X, A, B);
                    U = X;
                    //---------------------------------------------------
                }
                System.Diagnostics.Debug.WriteLine("Newton Raphson iter = " + iter + " norm = " + convRatio);
                System.Diagnostics.Debug.Assert(iter < maxIter);
            }
            else
            {
                var A = new IvyFEM.Linear.DoubleSparseMatrix(nodeCnt, nodeCnt);
                var B = new double[nodeCnt];

                CalcAB(A, B);

                DoubleSetFixedCadsCondtion(A, B, new int[] { vNodeCnt, pNodeCnt }, new int[] { vDof, pDof });
                double[] X;
                Solver.DoubleSolve(out X, A, B);
                U = X;
            }
        }
    }
}
