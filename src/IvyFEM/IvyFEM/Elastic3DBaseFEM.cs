using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class Elastic3DBaseFEM : FEM
    {
        public delegate void CalcElementDoubleAB(
           uint feId, IvyFEM.Linear.DoubleSparseMatrix A, double[] B);

        public double ConvRatioToleranceForNonlinearIter { get; set; }
            = 1.0e+2 * IvyFEM.Linear.Constants.ConvRatioTolerance; // 収束しないので収束条件を緩めている

        // Calc Matrix
        protected IList<CalcElementDoubleAB> CalcElementABs { get; set; } = new List<CalcElementDoubleAB>();

        public IList<uint> DisplacementQuantityIds { get; set; } = new List<uint> { 0 };

        //Solve
        // Output
        public double[] U { get; protected set; }

        protected int GetOffset(uint quantityId)
        {
            int cnt = 0;
            for (uint tmpId = 0; tmpId < quantityId; tmpId++)
            {
                int quantityDof = (int)World.GetDof(tmpId);
                int quantityNodeCnt = (int)World.GetNodeCount(tmpId);
                cnt += quantityDof * quantityNodeCnt;
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
        }

        protected void SetSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            SetExternalForceSpecialBC(A, B);
        }

        public override void Solve()
        {
            int quantityCnt = World.GetQuantityCount();
            int nodeCnt = 0;
            for (uint quantityId = 0; quantityId < quantityCnt; quantityId++)
            {
                int quantityDof = (int)World.GetDof(quantityId);
                int quantityNodeCnt = (int)World.GetNodeCount(quantityId);
                nodeCnt += quantityDof * quantityNodeCnt;
            }

            if (MustUseNonlinearIter())
            {
                // Nonlinear Iter
                double sqNorm = 0;
                double sqInvNorm0 = 0;
                double convRatio = ConvRatioToleranceForNonlinearIter;
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
                    SetSpecialBC(A, B);
                    System.Diagnostics.Debug.WriteLine("SetSpecialBC: t = " + (System.Environment.TickCount - t));

                    t = System.Environment.TickCount;
                    DoubleSetFixedCadsCondtion(A, B);
                    System.Diagnostics.Debug.WriteLine("Condition: t = " + (System.Environment.TickCount - t));

                    t = System.Environment.TickCount;
                    DoubleSetForceFixedCadsCondtion(A, B);
                    System.Diagnostics.Debug.WriteLine("Force Condition: t = " + (System.Environment.TickCount - t));

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
                        System.Diagnostics.Debug.WriteLine("cur convRatio =" + convRatio);
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
                System.Diagnostics.Debug.WriteLine("Nonlinear iter = " + iter + " norm = " + convRatio);
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
                SetSpecialBC(A, B);
                System.Diagnostics.Debug.WriteLine("SetSpecialBC: t = " + (System.Environment.TickCount - t));

                t = System.Environment.TickCount;
                DoubleSetFixedCadsCondtion(A, B);
                System.Diagnostics.Debug.WriteLine("Condtion: t = " + (System.Environment.TickCount - t));

                t = System.Environment.TickCount;
                DoubleSetForceFixedCadsCondtion(A, B);
                System.Diagnostics.Debug.WriteLine("Force Condtion: t = " + (System.Environment.TickCount - t));

                t = System.Environment.TickCount;
                //-------------------------------
                double[] X;
                Solver.DoubleSolve(out X, A, B);
                U = X;
                //-------------------------------
                System.Diagnostics.Debug.WriteLine("Solve: t = " + (System.Environment.TickCount - t));
            }

        }

        protected bool MustUseNonlinearIter()
        {
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
                if (ma is NullMaterial)
                {
                    // null
                }
                else if (ma is DKTPlateMaterial)
                {
                    // linear
                }
                else if (ma is MindlinPlateMaterial)
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
    }
}
