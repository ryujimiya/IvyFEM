using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public delegate void CalcElementDoubleKM(
        uint feId, IvyFEM.Lapack.DoubleMatrix K, IvyFEM.Lapack.DoubleMatrix M);

    public partial class Elastic2DEigenFEM : FEM
    {
        // Calc Matrix
        protected IList<CalcElementDoubleKM> CalcElementKMs { get; set; } = new List<CalcElementDoubleKM>();
        // Calc Matrix for truss
        protected IList<CalcElementDoubleKM> CalcElementKMsForLine { get; set; } = new List<CalcElementDoubleKM>();

        protected int ConstraintCount => 0;
        public IList<uint> DisplacementQuantityIds { get; set; } = new List<uint> { 0 };

        //Solve
        // Output
        public System.Numerics.Complex[] FrequencyZs { get; protected set; }
        public System.Numerics.Complex[][] EVecZs { get; protected set; }

        public Elastic2DEigenFEM(FEWorld world)
        {
            World = world;
            SetupCalcKMs();
        }

        protected void SetupCalcKMs()
        {
            CalcElementKMs.Clear();
            CalcElementKMsForLine.Clear();

            // Linear/Staint Venant

            // Hyperelastic

            // Truss/Beam
            CalcElementKMsForLine.Add(CalcBeamElementKMForLine);
        }

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

        protected void CalcKM(IvyFEM.Lapack.DoubleMatrix K, IvyFEM.Lapack.DoubleMatrix M)
        {
            uint quantityId = 0; // Note: 複数変数のときでも要素Idは同じはずなので0指定
            IList<uint> feIds = World.GetTriangleFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                foreach (var calcElementKM in CalcElementKMs)
                {
                    calcElementKM(feId, K, M);
                }
            }
            IList<uint> lineFEIds = World.GetLineFEIds(quantityId);
            foreach (uint feId in lineFEIds)
            {
                foreach (var calcElementKMForLine in CalcElementKMsForLine)
                {
                    calcElementKMForLine(feId, K, M);
                }
            }
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

            var A = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            var B = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            CalcKM(A, B);

            DoubleSetFixedCadsCondtionForEigen(A, B);

            System.Numerics.Complex[] eVals = null;
            System.Numerics.Complex[][] eVecs = null;
            int ret = -1;
            try
            {
                ret = IvyFEM.Lapack.Functions.dggev_dirty(A.Buffer, A.RowLength, A.ColumnLength,
                    B.Buffer, B.RowLength, B.ColumnLength,
                    out eVals, out eVecs);
                System.Diagnostics.Debug.Assert(ret == 0);
            }
            catch (InvalidOperationException exception)
            {
                //System.Diagnostics.Debug.Assert(false);
                System.Diagnostics.Debug.WriteLine("!!!!!!!ERROR!!!!!!!!!");
                System.Diagnostics.Debug.WriteLine(exception.Message);
                System.Diagnostics.Debug.WriteLine(exception.StackTrace);
                ret = -1;
            }
            if (ret != 0)
            {
                // fail safe
                int n = A.RowLength;
                eVals = new System.Numerics.Complex[n];
                eVecs = new System.Numerics.Complex[n][];
                for (int iMode = 0; iMode < n; iMode++)
                {
                    eVecs[iMode] = new System.Numerics.Complex[n];
                }
            }

            SortEVals(eVals, eVecs);
            AdjustPhaseEVecs(eVecs);

            int modeCnt = eVals.Length;
            var freqZs = new List<System.Numerics.Complex>();
            var eVecZs = new List<System.Numerics.Complex[]>();
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                System.Numerics.Complex omegaZ = System.Numerics.Complex.Sqrt(eVals[iMode]);
                System.Numerics.Complex freqZ = omegaZ / (2.0 * Math.PI);
                System.Numerics.Complex[] eVecZ = eVecs[iMode];

                freqZs.Add(freqZ);
                eVecZs.Add(eVecZ);
            }
            FrequencyZs = freqZs.ToArray();
            EVecZs = eVecZs.ToArray();
        }

        private void SortEVals(System.Numerics.Complex[] eVals, System.Numerics.Complex[][] eVecs)
        {
            int modeCnt = eVals.Length;
            var eValEVecs = new List<KeyValuePair<System.Numerics.Complex, System.Numerics.Complex[]>>();
            for (int i = 0; i < modeCnt; i++)
            {
                eValEVecs.Add(new KeyValuePair<System.Numerics.Complex, System.Numerics.Complex[]>(eVals[i], eVecs[i]));
            }
            eValEVecs.Sort((a, b) =>
            {
                // eVal(ω^2) の実部を比較
                double diff = a.Key.Real - b.Key.Real;
                // 昇順(ωの小さい方から)
                if (diff > 0)
                {
                    return 1;
                }
                else if (diff < 0)
                {
                    return -1;
                }
                return 0;
            });

            for (int i = 0; i < modeCnt; i++)
            {
                eVals[i] = eValEVecs[i].Key;
                eVecs[i] = eValEVecs[i].Value;
            }
        }

        private void AdjustPhaseEVecs(System.Numerics.Complex[][] eVecs)
        {
            int modeCnt = eVecs.Length;

            int quantityCnt = World.GetQuantityCount();
            int nodeCnt = 0;
            for (uint quantityId = 0; quantityId < quantityCnt; quantityId++)
            {
                int quantityDof = (int)World.GetDof(quantityId);
                int quantityNodeCnt = (int)World.GetNodeCount(quantityId);
                nodeCnt += quantityDof * quantityNodeCnt;
            }

            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var eVec = eVecs[iMode];
                // uの最大値を求める
                System.Numerics.Complex maxValue = new System.Numerics.Complex(0, 0);
                double maxAbs = 0;
                foreach (uint uQuantityId in DisplacementQuantityIds)
                {
                    int uDof = (int)World.GetDof(uQuantityId);
                    int uNodeCnt = (int)World.GetNodeCount(uQuantityId);
                    int offset = GetOffset(uQuantityId);

                    for (int iNode = 0; iNode < uNodeCnt * uDof; iNode++)
                    {
                        System.Numerics.Complex value = eVec[offset + iNode];
                        double abs = value.Magnitude;
                        if (abs > maxAbs)
                        {
                            maxAbs = abs;
                            maxValue = value;
                        }
                    }
                }
                System.Numerics.Complex phase = maxValue / maxAbs;

                for (int i = 0; i < eVec.Length; i++)
                {
                    eVec[i] /= phase;
                }
            }
        }
    }
}
