using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class FluidRANS2DTDFEM : FEM
    {
        public double TimeStep { get; private set; } = 0;
        public double NewmarkBeta { get; private set; } = 1.0 / 4.0;
        public double NewmarkGamma { get; private set; } = 1.0 / 2.0;
        public uint ValueId { get; private set; } = 0;
        public uint PrevValueId { get; private set; } = 0;
        public uint KpValueId { get; private set; } = 0;
        public uint PrevKpValueId { get; private set; } = 0;
        public uint EpValueId { get; private set; } = 0;
        public uint PrevEpValueId { get; private set; } = 0;

        public double ConvRatioToleranceForNonlinearIter { get; set; }
            = 1.0e+2 * IvyFEM.Linear.Constants.ConvRatioTolerance; // 収束しないので収束条件を緩めている

        // output
        public double[] U { get; protected set; } = null;

        public FluidRANS2DTDFEM(FEWorld world,
            double timeStep,
            double newmarkBeta, double newmarkGamma,
            uint valueId, uint prevValueId,
            uint kpValueId, uint prevKpValueId,
            uint epValueId, uint prevEpValueId)
        {
            World = world;
            TimeStep = timeStep;
            NewmarkBeta = newmarkBeta;
            NewmarkGamma = newmarkGamma;
            ValueId = valueId;
            PrevValueId = prevValueId;
            KpValueId = kpValueId;
            PrevKpValueId = prevKpValueId;
            EpValueId = epValueId;
            PrevEpValueId = prevEpValueId;
        }

        public void UpdateFieldValuesTimeDomain()
        {
            UpdateFieldValuesNewmarkBetaTimeDomain(
                U, ValueId, PrevValueId,
                TimeStep,
                NewmarkBeta, NewmarkGamma);
            UpdateFieldValuesNewmarkBetaTimeDomain(
                U, KpValueId, PrevKpValueId,
                TimeStep,
                NewmarkBeta, NewmarkGamma);
            UpdateFieldValuesNewmarkBetaTimeDomain(
                U, EpValueId, PrevEpValueId,
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

            U = new double[nodeCnt];
            SolveRANS();
        }

        private void SolveRANS()
        {
            uint vQuantityId = 0;
            uint pQuantityId = 1;
            uint kpQuantityId = 2;
            uint epQuantityId = 3;
            int vDof = 2;
            int pDof = 1;
            int kpDof = 1;
            int epDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int kpNodeCnt = (int)World.GetNodeCount(kpQuantityId);
            int epNodeCnt = (int)World.GetNodeCount(epQuantityId);
            int pOffset = vNodeCnt * vDof;
            int kpOffset = pOffset + pNodeCnt;
            int epOffset = kpOffset + kpNodeCnt;

            int nodeCnt = vNodeCnt * vDof + pNodeCnt +kpNodeCnt + epNodeCnt;
            System.Diagnostics.Debug.Assert(U.Length == nodeCnt);

            ////////////////////////////////////////////////////
            // Nonlinear Iter
            double sqNorm = 0;
            double sqInvNorm0 = 0;
            double convRatio = ConvRatioToleranceForNonlinearIter;
            double tolerance = convRatio;
            const int maxIter = IvyFEM.Linear.Constants.MaxIter;
            int iter = 0;
            for (iter = 0; iter < maxIter; iter++)
            {
                var A = new IvyFEM.Linear.DoubleSparseMatrix(nodeCnt, nodeCnt);
                var B = new double[nodeCnt];

                // RANS
                //CalcABStandardKEpsilon(A, B);
                CalcABRNGKEpsilon(A, B);

                DoubleSetFixedCadsCondtion(A, B);

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
                X.CopyTo(U, 0);
                //---------------------------------------------------
            }
            System.Diagnostics.Debug.WriteLine("(RANS) Nonlinear iter = " + iter + " norm = " + convRatio);
            System.Diagnostics.Debug.Assert(iter < maxIter);
        }
    }
}
