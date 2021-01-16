using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class Elastic3DBaseFEM : FEM
    {
        // for incremental solution
        public delegate void SetNodeValues(uint quantityId, int coId);
        public delegate void SetElementValues(uint feId);

        public delegate void CalcElementDoubleAB(
           uint feId, IvyFEM.Linear.DoubleSparseMatrix A, double[] B);

        public double ConvRatioToleranceForNonlinearIter { get; set; }
            = 1.0e+2 * IvyFEM.Linear.Constants.ConvRatioTolerance; // 収束しないので収束条件を緩めている


        ////////////////////////////////////////////////////////////
        // 3D
        // Calc Matrix
        protected IList<CalcElementDoubleAB> CalcElementABs { get; set; } = new List<CalcElementDoubleAB>();

        ////////////////////////////////////////////////////////////
        // 2D
        // Init/Update
        public bool IsUseInit { get; set; } = false;
        public bool IsUseUpdate { get; set; } = false;
        public int TimeIndexForInit { get; set; } = -1;
        // Init/Update
        public IList<uint> AdditionalValueIds { get; set; } = new List<uint>();

        // Init
        protected IList<SetNodeValues> InitNodeValuessForPlate { get; set; } = new List<SetNodeValues>();
        protected IList<SetElementValues> InitElementValuessForPlate { get; set; } = new List<SetElementValues>();

        // Update
        protected IList<SetNodeValues> UpdateNodeValuessForPlate { get; set; } = new List<SetNodeValues>();
        protected IList<SetElementValues> UpdateElementValuessForPlate { get; set; } = new List<SetElementValues>();

        // Calc Matrix
        protected IList<CalcElementDoubleAB> CalcElementABsForPlate { get; set; } = new List<CalcElementDoubleAB>();

        protected int ConstraintCount { get => HasMultipointConstraints() ? 1 : 0; }

        //Solve
        // Output
        public double[] U { get; protected set; }

        protected void InitValues()
        {
            if (!IsUseInit)
            {
                return;
            }

            ////////////////////////////////////////////////////
            // 3D
            // なし

            ////////////////////////////////////////////////////
            // 2D
            // node
            if (InitNodeValuessForPlate.Count > 0)
            {
                int quantityCnt = World.GetQuantityCount();
                for (uint quantityId = 0; quantityId < quantityCnt; quantityId++)
                {
                    int coCnt = (int)World.GetCoordCount(quantityId);
                    for (int coId = 0; coId < coCnt; coId++)
                    {
                        foreach (var initNodeValuesForPlate in InitNodeValuessForPlate)
                        {
                            initNodeValuesForPlate(quantityId, coId);
                        }
                    }
                }
            }
            // element
            if (InitElementValuessForPlate.Count > 0)
            {
                uint quantityId = 0; // Note: 複数変数のときでも要素Idは同じはずなので0指定
                IList<uint> feIds = World.GetTriangleFEIds(quantityId);
                foreach (uint feId in feIds)
                {
                    foreach (var initElementValuesForPlate in InitElementValuessForPlate)
                    {
                        initElementValuesForPlate(feId);
                    }
                }
            }
        }

        protected void UpdateValues()
        {
            if (!IsUseUpdate)
            {
                return;
            }

            ////////////////////////////////////////////////////
            // 3D
            // なし

            ////////////////////////////////////////////////////
            // 2D
            // node
            if (UpdateNodeValuessForPlate.Count > 0)
            {
                int quantityCnt = World.GetQuantityCount();
                for (uint quantityId = 0; quantityId < quantityCnt; quantityId++)
                {
                    int coCnt = (int)World.GetCoordCount(quantityId);
                    for (int coId = 0; coId < coCnt; coId++)
                    {
                        foreach (var updateNodeValuesForPlate in UpdateNodeValuessForPlate)
                        {
                            updateNodeValuesForPlate(quantityId, coId);
                        }
                    }
                }
            }
            // element
            if (UpdateElementValuessForPlate.Count > 0)
            {
                uint quantityId = 0; // Note: 複数変数のときでも要素Idは同じはずなので0指定
                IList<uint> feIds = World.GetTriangleFEIds(quantityId);
                foreach (uint feId in feIds)
                {
                    foreach (var updateElementValuesForPlate in UpdateElementValuessForPlate)
                    {
                        updateElementValuesForPlate(feId);
                    }
                }
            }
        }

        protected void CalcAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint quantityId = 0; // Note: 複数変数のときでも要素Idは同じはずなので0指定
            if (CalcElementABs.Count > 0)
            {
                IList<uint> feIds = World.GetTetrahedronFEIds(quantityId);
                foreach (uint feId in feIds)
                {
                    foreach (var calcElementAB in CalcElementABs)
                    {
                        calcElementAB(feId, A, B);
                    }
                }
            }
            if (CalcElementABsForPlate.Count > 0)
            {
                IList<uint> triFEIds = World.GetTriangleFEIds(quantityId);
                foreach (uint feId in triFEIds)
                {
                    foreach (var calcElementABForPlate in CalcElementABsForPlate)
                    {
                        calcElementABForPlate(feId, A, B);
                    }
                }
            }
        }

        protected void SetSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            SetExternalForceSpecialBC(A, B);
            SetMultipointConstraintSpecialBC(A, B);
        }

        public override void Solve()
        {
            // Init
            InitValues();

            // Solve
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

            // Update
            UpdateValues();
        }

        protected bool MustUseNonlinearIter()
        {
            if (HasMultipointConstraints())
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
            IList<uint> feIds = World.GetTetrahedronFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = World.GetTetrahedronFE(quantityId, feId);
                Material ma = World.GetMaterial(tetFE.MaterialId);
                if (ma is NullMaterial)
                {
                    // null
                }
                else if (ma is LinearElasticMaterial)
                {
                    // linear
                }
                else
                {
                    hasNonlinear = true;
                    break;
                }
            }

            if (!hasNonlinear)
            {
                IList<uint> triFEIds = World.GetTriangleFEIds(quantityId);
                foreach (uint feId in triFEIds)
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
                    else if (ma is MITCLinearPlateMaterial)
                    {
                        // linear
                    }
                    else if (ma is MITCStVenantPlateMaterial)
                    {
                        // 非線形問題だが、方程式は線形化できる
                        // linear
                    }
                    else
                    {
                        hasNonlinear = true;
                        break;
                    }
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
    }
}
