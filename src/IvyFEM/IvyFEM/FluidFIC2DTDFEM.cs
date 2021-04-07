using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class FluidFIC2DTDFEM : FEM
    {
        public double TimeStep { get; private set; } = 0;
        public double NewmarkBeta { get; private set; } = 1.0 / 4.0;
        public double NewmarkGamma { get; private set; } = 1.0 / 2.0;
        public uint ValueId { get; private set; } = 0;
        public uint PrevValueId { get; private set; } = 0;

        public double ConvRatioToleranceForNonlinearIter { get; set; }
            = 1.0e+2 * IvyFEM.Linear.Constants.ConvRatioTolerance; // 収束しないので収束条件を緩めている

        public int Counter { get; set; } = 0;

        // output
        public double[] U { get; protected set; } = null;
        public double[] UpdatedCoord { get; set; } = null;

        public FluidFIC2DTDFEM(FEWorld world,
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
            UpdateFieldValuesNewmarkBetaTimeDomain(
                U, ValueId, PrevValueId,
                TimeStep,
                NewmarkBeta, NewmarkGamma);
        }

        public override void Solve()
        {
            int nodeCnt = 0;
            int quantityCnt = World.GetQuantityCount();
            for (uint quantityId = 0; quantityId < quantityCnt; quantityId++)
            {
                int quantityDof = (int)World.GetDof(quantityId);
                int quantityNodeCnt = (int)World.GetNodeCount(quantityId);
                nodeCnt += quantityNodeCnt * quantityDof;
            }

            ////////////////////////////////////////////////////
            // Nonlinear Iter
            double sqNorm = 0;
            double sqInvNorm0 = 0;
            double convRatio = ConvRatioToleranceForNonlinearIter;
            double tolerance = convRatio;
            const int maxIter = IvyFEM.Linear.Constants.MaxIter;
            int iter = 0;

            U = new double[nodeCnt];
            if (Counter == 0)
            {
                InitFIC();
            }

            for (iter = 0; iter < maxIter; iter++)
            {
                var A = new IvyFEM.Linear.DoubleSparseMatrix(nodeCnt, nodeCnt);
                var B = new double[nodeCnt];

                UpdateFEDisplacements();
                CalcAB(A, B);
                ClearFEDisplacements();

                SetSpecialBC(A, B);

                DoubleSetFixedCadsCondtion(A, B);

                DoubleSetForceFixedCadsCondtion(A, B);

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

            UpdateFIC();

            System.Diagnostics.Debug.WriteLine("Nonlinear iter = " + iter + " norm = " + convRatio);
            System.Diagnostics.Debug.Assert(iter < maxIter);
            ////////////////////////////////////////////////////
        }

        protected void CalcAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            CalcABFIC(A, B);
        }

        protected void SetSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            SetMultipointConstraintSpecialBC(A, B);
        }
    }
}
