using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public delegate void CalcElementDoubleAB(
        uint feId, IvyFEM.Linear.DoubleSparseMatrix A, double[] B);

    abstract public partial class Elastic2DBaseFEM : FEM
    {
        protected int[] NodeCounts { get; set; } = null;
        protected int[] Dofs { get; set; } = null;
        public double ConvRatioToleranceForNewtonRaphson { get; set; }
            = 1.0e+2 * IvyFEM.Linear.Constants.ConvRatioTolerance; // 収束しないので収束条件を緩めている
        // Calc Matrix
        protected IList<CalcElementDoubleAB> CalcElementABs { get; set; } = new List<CalcElementDoubleAB>();

        //Solve
        // Output
        public double[] U { get; protected set; }

        protected void SetupNodeCount()
        {
            int quantityCnt = World.GetQuantityCount();
            Dofs = new int[quantityCnt];
            NodeCounts = new int[quantityCnt];
            for (uint quantityId = 0; quantityId < quantityCnt; quantityId++)
            {
                Dofs[quantityId] = (int)World.GetDof(quantityId);
                NodeCounts[quantityId] = (int)World.GetNodeCount(quantityId);
            }
        }

        protected int GetNodeCount()
        {
            if (Dofs == null)
            {
                return 0;
            }
            if (NodeCounts == null)
            {
                return 0;
            }
            int cnt = 0;
            for (int i = 0; i < Dofs.Length; i++)
            {
                cnt += NodeCounts[i] * Dofs[i];
            }
            return cnt;
        }

        protected int GetOffset(uint quantityIdIndex)
        {
            if (Dofs == null)
            {
                return 0;
            }
            if (NodeCounts == null)
            {
                return 0;
            }
            int cnt = 0;
            for (int i = 0; i < quantityIdIndex; i++)
            {
                cnt += NodeCounts[i] * Dofs[i];
            }
            return cnt;
        }

        protected void CalcAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint quantityId = 0; // Note: 複数変数のときでも要素Idは同じはずなので0指定
            IList<uint> feIds = World.GetTriangleFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                foreach (var calcElementAB in CalcElementABs)
                {
                    calcElementAB(feId, A, B);
                }
            }
            CalcMultipointConstraintAB(A, B);
            CalcTwoBodyContactAB(A, B);
        }

        public override void Solve()
        {
            SetupNodeCount();

            int nodeCnt = GetNodeCount();

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
                    int t;
                    var A = new IvyFEM.Linear.DoubleSparseMatrix(nodeCnt, nodeCnt);
                    var B = new double[nodeCnt];

                    t = System.Environment.TickCount;
                    CalcAB(A, B);
                    System.Diagnostics.Debug.WriteLine("CalcAB: t = " + (System.Environment.TickCount - t));

                    t = System.Environment.TickCount;
                    DoubleSetFixedCadsCondtion(A, B, NodeCounts, Dofs);
                    System.Diagnostics.Debug.WriteLine("Condition: t = " + (System.Environment.TickCount - t));

                    t = System.Environment.TickCount;
                    double[] AU = A * U;
                    double[] R = IvyFEM.Lapack.Functions.daxpy(-1.0, B, AU);
                    sqNorm = IvyFEM.Lapack.Functions.ddot(R, R);
                    System.Diagnostics.Debug.WriteLine("Calc Norm: t = " + (System.Environment.TickCount - t));
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
                        if (sqNorm * sqInvNorm0 < tolerance * tolerance)
                        {
                            break;
                        }
                    }

                    t = System.Environment.TickCount;
                    //---------------------------------------------------
                    double[] X;
                    Solver.DoubleSolve(out X, A, B);
                    U = X;
                    //---------------------------------------------------
                    System.Diagnostics.Debug.WriteLine("Solve: t = " + (System.Environment.TickCount - t));
                }
                System.Diagnostics.Debug.WriteLine("Newton Raphson iter = " + iter + " norm = " + convRatio);
                System.Diagnostics.Debug.Assert(iter < maxIter);
            }
            else
            {
                int t;
                var A = new IvyFEM.Linear.DoubleSparseMatrix(nodeCnt, nodeCnt);
                var B = new double[nodeCnt];

                t = System.Environment.TickCount;
                CalcAB(A, B);
                System.Diagnostics.Debug.WriteLine("CalcAB: t = " + (System.Environment.TickCount - t));

                t = System.Environment.TickCount;
                DoubleSetFixedCadsCondtion(A, B, NodeCounts, Dofs);
                System.Diagnostics.Debug.WriteLine("Condtion: t = " + (System.Environment.TickCount - t));

                t = System.Environment.TickCount;
                //-------------------------------
                double[] X;
                Solver.DoubleSolve(out X, A, B);
                U = X;
                //-------------------------------
                System.Diagnostics.Debug.WriteLine("Solve: t = " + (System.Environment.TickCount - t));
            }

        }

        protected bool MustUseNewtonRaphson()
        {
            if (HasMultipointConstraints())
            {
                return true;
            }
            if (HasTwoBodyContact())
            {
                return true;
            }
            if (HasNonLinearElasticMaterial())
            {
                return true;
            }
            return false;
        }

        private bool HasNonLinearElasticMaterial()
        {
            bool hasNonlinear = false;
            uint quantityId = 0; // Note: 複数変数のときでも要素Idは同じはずなので0指定
            IList<uint> feIds = World.GetTriangleFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = World.GetTriangleFE(quantityId, feId);
                Material ma = World.GetMaterial(triFE.MaterialId);
                if (ma is LinearElasticMaterial)
                {
                    // linear
                }
                else
                {
                    hasNonlinear = true;
                    break;
                }
            }
            return hasNonlinear;
        }

        private bool HasMultipointConstraints()
        {
            bool hasMPC = false;
            for (uint quantityId = 0; quantityId < World.GetQuantityCount(); quantityId++)
            {
                int cnt = World.GetMultipointConstraintCount(quantityId);
                if (cnt > 0)
                {
                    hasMPC = true;
                    break;
                }
            }
            return hasMPC;
        }

        private bool HasTwoBodyContact()
        {
            bool hasContact = false;
            for (uint quantityId = 0; quantityId < World.GetQuantityCount(); quantityId++)
            {
                int slaveCnt = World.GetContactSlaveEIds(quantityId).Count;
                int masterCnt = World.GetContactMasterEIds(quantityId).Count;
                if (slaveCnt > 0 && masterCnt > 0)
                {
                    hasContact = true;
                    break;
                }
            }
            return hasContact;
        }
    }
}
