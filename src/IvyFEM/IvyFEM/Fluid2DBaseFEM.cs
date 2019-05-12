using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class Fluid2DBaseFEM : FEM
    {
        public double ConvRatioToleranceForNewtonRaphson { get; set; }
            = 1.0e+2 * IvyFEM.Linear.Constants.ConvRatioTolerance; // 収束しないので収束条件を緩めている

        /// <summary>
        /// 方程式のタイプ
        /// </summary>
        public FluidEquationType EquationType { get; set; } = FluidEquationType.StdGNavierStokes;

        // output
        public double[] U { get; protected set; } = null;
        public double[] CoordV { get; protected set; } = null; // Vorticityの場合の速度を格納

        protected abstract void CalcAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B);
        protected abstract void SetSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B);
        protected abstract void PostSolve();
        protected abstract bool MustUseNewtonRaphson();

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

                    PostSolve();

                    CalcAB(A, B);

                    SetSpecialBC(A, B);

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

                SetSpecialBC(A, B);

                DoubleSetFixedCadsCondtion(A, B);
                double[] X;
                Solver.DoubleSolve(out X, A, B);
                U = X;

                PostSolve();
            }
        }
    }
}
