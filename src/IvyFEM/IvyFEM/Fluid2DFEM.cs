using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Fluid2DFEM : FEM
    {
        public double ConvRatioToleranceForNewtonRaphson { get; set; }
            = 1.0e+2 * IvyFEM.Linear.Constants.ConvRatioTolerance; // 収束しないので収束条件を緩めている

        // output
        public double[] U { get; private set; } = null;

        public Fluid2DFEM(FEWorld world)
        {
            World = world;
        }

        private void CalcAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            CalcStandardGalerkinAB(A, B);
            //CalcSUPGAB(A, B);
        }

        public override void Solve()
        {
            // Newton Raphson
            double sqNorm = 0;
            double sqInvNorm0 = 0;
            double convRatio = ConvRatioToleranceForNewtonRaphson;
            double tolerance = convRatio;
            const int maxIter = IvyFEM.Linear.Constants.MaxIter;
            const int minIter = 2;
            int iter = 0;

            uint vQuantityId = 0;
            int vDof = 2;
            uint pQuantityId = 1;
            int pDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int nodeCnt = vNodeCnt * vDof + pNodeCnt * pDof;
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
                    // 1回目で収束してしまうのを防ぐ
                    if (iter >= minIter && sqNorm * sqInvNorm0 < tolerance * tolerance)
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
    }
}
