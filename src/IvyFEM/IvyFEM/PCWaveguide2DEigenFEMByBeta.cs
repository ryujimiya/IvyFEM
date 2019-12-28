using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    // βを与えてk0を解く
    public class PCWaveguide2DEigenFEMByBeta : FEM
    {
        public uint QuantityId { get; private set; } = 0;
        public uint PortId { get; private set; } = 0;
        public bool IsTMMode { get; set; } = false;
        public PCWaveguidePortInfo WgPortInfo { get; private set; }

        private IvyFEM.Lapack.DoubleMatrix KMat0 = null;
        private IvyFEM.Lapack.DoubleMatrix MMat0 = null;

        // Solve
        // Input
        public double BetaX { get; set; } = 0.0;
        public double BetaY { get; set; } = 0.0; // for photonic band
        // Output
        public double[] Frequencys { get; private set; }
        public System.Numerics.Complex[][] EVecs { get; private set; }

        public PCWaveguide2DEigenFEMByBeta(
            FEWorld world, uint quantityId, uint portId, PCWaveguidePortInfo wgPortInfo)
        {
            World = world;
            QuantityId = quantityId;
            PortId = portId;
            WgPortInfo = wgPortInfo;
        }

        private void CalcMatrixs()
        {
            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            IList<uint> feIds = World.GetPeriodicPortTriangleFEIds(QuantityId, PortId);

            KMat0 = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            MMat0 = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);

            foreach (uint feId in feIds)
            {
                TriangleFE triFE = World.GetTriangleFE(QuantityId, feId);
                uint elemNodeCnt = triFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = triFE.NodeCoordIds[iNode];
                    int nodeId = World.PortCoord2Node(QuantityId, PortId, coId);
                    nodes[iNode] = nodeId;
                }

                Material ma0 = World.GetMaterial(triFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is DielectricMaterial);
                var ma = ma0 as DielectricMaterial;
                double maPxx = 0;
                double maPyy = 0;
                double maQzz = 0;
                if (IsTMMode)
                {
                    // TMモード
                    maPxx = 1.0 / ma.Epxx;
                    maPyy = 1.0 / ma.Epyy;
                    maQzz = ma.Muzz;
                }
                else
                {
                    // TEモード
                    maPxx = 1.0 / ma.Muxx;
                    maPyy = 1.0 / ma.Muyy;
                    maQzz = ma.Epzz;
                }

                double[,] sNN = triFE.CalcSNN();
                double[,][,] sNuNv = triFE.CalcSNuNv();
                double[,] sNxNx = sNuNv[0, 0];
                double[,] sNyNy = sNuNv[1, 1];
                double[][,] sNuN = triFE.CalcSNuN();
                double[,] sNxN = sNuN[0];
                double[,] sNyN = sNuN[1];

                // 要素剛性行列、要素質量行列を作る
                //  { [K]e - jβ[C]e - β^2[M]e }{Φ}= {0}
                for (int row = 0; row < elemNodeCnt; row++)
                {
                    int rowNodeId = nodes[row];
                    if (rowNodeId == -1)
                    {
                        continue;
                    }
                    for (int col = 0; col < elemNodeCnt; col++)
                    {
                        int colNodeId = nodes[col];
                        if (colNodeId == -1)
                        {
                            continue;
                        }
                        double kVal = maPxx * sNyNy[row, col] + maPyy * sNxNx[row, col];
                        double mVal = maQzz * sNN[row, col];

                        KMat0[rowNodeId, colNodeId] += kVal;
                        MMat0[rowNodeId, colNodeId] += mVal;
                    }
                }
            }
        }

        public override void Solve()
        {
            Frequencys = null;
            EVecs = null;

            System.Diagnostics.Debug.Assert(!WgPortInfo.IsSVEA);


            if (WgPortInfo.BcEdgeIds3.Count == 0)
            {
                // for waveguide
                SolveWaveguideEigen();
            }
            else
            {
                // for photonic band
                SolvePhotonicBand();
            }
        }

        // 導波路固有値問題
        private void SolveWaveguideEigen()
        {
            // input
            double beta = BetaX;

            double minFreq = WgPortInfo.MinFrequency;
            double maxFreq = WgPortInfo.MaxFrequency;
            double minK0 = 0;
            double maxK0 = double.MaxValue;
            minK0 = (2.0 * Math.PI / Constants.C0) * minFreq;
            if (maxFreq >= double.MaxValue)
            {
                maxK0 = double.MaxValue;
            }
            else
            {
                maxK0 = (2.0 * Math.PI / Constants.C0) * maxFreq;
            }

            // 領域
            CalcMatrixs();

            double[] eVals;
            System.Numerics.Complex[][] eVecs;
            // 導波路 (Φを直接解く)
            SolveWaveguideNonSVEAModeAsQuadraticGeneralizedEigen(beta, out eVals, out eVecs);
            SortEVals(eVals, eVecs);
            AdjustPhaseEVecs(eVecs);

            double[] defectK0s;
            System.Numerics.Complex[][] defectEVecs;
            GetDefectModes(
                minK0, maxK0,
                eVals, eVecs,
                out defectK0s, out defectEVecs);
            // モード追跡
            TraceMode(defectK0s, defectEVecs);

            Frequencys = new double[defectK0s.Length];
            for (int i = 0; i < defectK0s.Length; i++)
            {
                double k0 = defectK0s[i];
                double omega = k0 * Constants.C0;
                double freq = omega / (2.0 * Math.PI);
                Frequencys[i] = freq;
            }
            EVecs = defectEVecs;
        }

        // Φを直接解く
        private void SolveWaveguideNonSVEAModeAsQuadraticGeneralizedEigen(
            double beta,
            out double[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            eVals = null;
            eVecs = null;

            int nodeCnt = KMat0.RowLength;
            int bcNodeCnt = WgPortInfo.BcNodess[0].Count;
            System.Diagnostics.Debug.Assert(WgPortInfo.BcNodess[0].Count == WgPortInfo.BcNodess[1].Count);
            int nodeCnt1 = nodeCnt - bcNodeCnt; // 境界2を除く
            bool isPortBc2Reverse = WgPortInfo.IsPortBc2Reverse;
            double periodicDistance = WgPortInfo.PeriodicDistanceX;

            int innerNodeCnt = nodeCnt - bcNodeCnt * 2;
            var KMat = new IvyFEM.Lapack.ComplexMatrix(nodeCnt1, nodeCnt1);
            var MMat = new IvyFEM.Lapack.ComplexMatrix(nodeCnt1, nodeCnt1);

            System.Numerics.Complex expA = 
                System.Numerics.Complex.Exp(-1.0 * System.Numerics.Complex.ImaginaryOne * beta * periodicDistance);

            for (int i = 0; i < bcNodeCnt; i++)
            {
                int iB2 = isPortBc2Reverse ? (bcNodeCnt * 2 - 1 - i) : (bcNodeCnt + i);
                for (int j = 0; j < bcNodeCnt; j++)
                {
                    int jB2 = isPortBc2Reverse ? (bcNodeCnt * 2 - 1 - j) : (bcNodeCnt + j);
                    // [K11]
                    KMat[i, j] += KMat0[i, j];
                    MMat[i, j] += MMat0[i, j];
                    // [K12] ->[K11]
                    KMat[i, j] += expA * KMat0[i, jB2];
                    MMat[i, j] += expA * MMat0[i, jB2];
                    // [K21] ->[K11]
                    KMat[i, j] += (1.0 / expA) * KMat0[iB2, j];
                    MMat[i, j] += (1.0 / expA) * MMat0[iB2, j];
                    // [K22] ->[K11]
                    // Note: (1.0 /expA) * expA
                    KMat[i, j] += KMat0[iB2, jB2];
                    MMat[i, j] += MMat0[iB2, jB2];
                }
                for (int j = 0; j < innerNodeCnt; j++)
                {
                    // [K10]
                    KMat[i, j + bcNodeCnt] += KMat0[i, j + bcNodeCnt * 2];
                    MMat[i, j + bcNodeCnt] += MMat0[i, j + bcNodeCnt * 2];
                    // [K20] -> [K10]
                    KMat[i, j + bcNodeCnt] += (1.0 / expA) * KMat0[iB2, j + bcNodeCnt * 2];
                    MMat[i, j + bcNodeCnt] += (1.0 / expA) * MMat0[iB2, j + bcNodeCnt * 2];
                }
            }
            for (int i = 0; i < innerNodeCnt; i++)
            {
                for (int j = 0; j < bcNodeCnt; j++)
                {
                    int jB2 = isPortBc2Reverse ? (bcNodeCnt * 2 - 1 - j) : (bcNodeCnt + j);
                    // [K01]
                    KMat[i + bcNodeCnt, j] += KMat0[i + bcNodeCnt * 2, j];
                    MMat[i + bcNodeCnt, j] += MMat0[i + bcNodeCnt * 2, j];
                    // [K02] -> [K01]
                    KMat[i + bcNodeCnt, j] += expA * KMat0[i + bcNodeCnt * 2, jB2];
                    MMat[i + bcNodeCnt, j] += expA * MMat0[i + bcNodeCnt * 2, jB2];
                }
                for (int j = 0; j < innerNodeCnt; j++)
                {
                    // [K00]
                    KMat[i + bcNodeCnt, j + bcNodeCnt] += KMat0[i + bcNodeCnt * 2, j + bcNodeCnt * 2];
                    MMat[i + bcNodeCnt, j + bcNodeCnt] += MMat0[i + bcNodeCnt * 2, j + bcNodeCnt * 2];
                }
            }

            System.Diagnostics.Debug.WriteLine("solve eigen");
            double[] eVals0;
            System.Numerics.Complex[][] eVecs0;
            SolveComplexHermiteBandGeneralizedEigen(KMat, MMat, out eVals0, out eVecs0);

            int modeCnt = eVals0.Length;
            // B2を復元
            eVecs = new System.Numerics.Complex[modeCnt][];
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                System.Numerics.Complex[] eVec = new System.Numerics.Complex[nodeCnt];
                eVecs[iMode] = eVec;

                System.Numerics.Complex[] eVec0 = eVecs0[iMode];
                System.Diagnostics.Debug.Assert(eVec0.Length == nodeCnt1);
                for (int i = 0; i < bcNodeCnt; i++)
                {
                    int iB2 = isPortBc2Reverse ? (bcNodeCnt * 2 - 1 - i) : (bcNodeCnt + i);

                    System.Numerics.Complex value = eVec0[i];
                    eVec[i] = value;
                    eVec[iB2] = expA * value;
                }
                for (int i = 0; i < innerNodeCnt; i++)
                {
                    eVec[i + bcNodeCnt * 2] = eVec0[i + bcNodeCnt];
                }

            }
            // k0に変換
            {
                eVals = new double[modeCnt];
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    // λ = k0^2と置いた場合( k0 = sqrt(λ) )
                    eVals[iMode] = Math.Sqrt(eVals0[iMode]);
                }
            }
        }

        private void SolveComplexHermiteBandGeneralizedEigen(
            IvyFEM.Lapack.ComplexMatrix A, IvyFEM.Lapack.ComplexMatrix B,
            out double[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            eVals = null;
            eVecs = null;

            int matLen = A.RowLength;
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            System.Diagnostics.Debug.Assert(B.RowLength == B.ColumnLength);
            System.Diagnostics.Debug.Assert(A.RowLength == B.RowLength);

            // 行列要素check
            {
                double th = 1.0e-12;
                for (int i = 0; i < matLen; i++)
                {
                    // [B]の正定値行列チェック
                    System.Diagnostics.Debug.Assert(B[i, i].Real > 0);
                    for (int j = i; j < matLen; j++)
                    {
                        // [A]はエルミート行列
                        System.Diagnostics.Debug.Assert(
                            (A[i, j] - System.Numerics.Complex.Conjugate(A[j, i])).Magnitude < th);
                        // [B]はエルミート行列
                        System.Diagnostics.Debug.Assert(
                            (B[i, j] - System.Numerics.Complex.Conjugate(B[j, i])).Magnitude < th);
                    }
                }
            }

            // Bのバンド幅を縮小
            System.Numerics.Complex[] dummyVec1 = new System.Numerics.Complex[matLen];
            IvyFEM.Linear.ComplexSparseMatrix sparseB1 = new IvyFEM.Linear.ComplexSparseMatrix(B);
            IvyFEM.Linear.ComplexSparseMatrix sparseB2;
            System.Numerics.Complex[] dummyVec2;
            int[] indexs;
            {
                bool successOrder = IvyFEM.Linear.Utils.OrderToComplexBandMatrix(
                    out sparseB2, out dummyVec2, out indexs, sparseB1, dummyVec1);
                System.Diagnostics.Debug.Assert(successOrder);
                if (!successOrder)
                {
                    System.Diagnostics.Debug.Assert(false);
                    return;
                }
            }
            // AはBの並び替えにあわせる
            IvyFEM.Linear.ComplexSparseMatrix sparseA2 =
                new IvyFEM.Linear.ComplexSparseMatrix(matLen, matLen);
            for (int row = 0; row < matLen; row++)
            {
                for (int col = 0; col < matLen; col++)
                {
                    sparseA2[row, col] = A[indexs[row], indexs[col]];
                }
            }

            var hermiteBandA = new IvyFEM.Lapack.ComplexHermitianBandMatrix(sparseA2);
            var hermiteBandB = new IvyFEM.Lapack.ComplexHermitianBandMatrix(sparseB2);

            System.Numerics.Complex[][] eVecs2 = null;
            int ret = -1;
            try
            {
                ret = IvyFEM.Lapack.Functions.zhbgv(
                            hermiteBandA.Buffer,
                            hermiteBandA.RowLength, hermiteBandA.ColumnLength, hermiteBandA.SuperdiaLength,
                            hermiteBandB.Buffer,
                            hermiteBandB.RowLength, hermiteBandB.ColumnLength, hermiteBandB.SuperdiaLength,
                            out eVals, out eVecs2);

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
                eVals = new double[n];
                eVecs2 = new System.Numerics.Complex[n][];
                for (int iMode = 0; iMode < n; iMode++)
                {
                    eVecs2[iMode] = new System.Numerics.Complex[n];
                }
            }

            // eVecs2の各ベクトルを元の番号順に戻す
            int modeCnt = eVals.Length;
            eVecs = new System.Numerics.Complex[modeCnt][];
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                System.Numerics.Complex[] eVec2 = eVecs2[iMode];
                System.Numerics.Complex[] eVec = new System.Numerics.Complex[matLen];
                eVecs[iMode] = eVec;
                for (int i = 0; i < matLen; i++)
                {
                    eVec[indexs[i]] = eVec2[i];
                }
            }
        }

        private void TraceMode(double[] defectK0s, System.Numerics.Complex[][] defectEVecs)
        {
            // モード追跡
            System.Numerics.Complex[][] prevModeEVecs = WgPortInfo.PrevModeEVecs;
            if (prevModeEVecs != null)
            {
                int prevModeCnt = prevModeEVecs.Length;
                IList<int> traceModeIndexss = new List<int>();
                for (int iPrevMode = 0; iPrevMode < prevModeCnt; iPrevMode++)
                {
                    System.Numerics.Complex[] prevEVec = prevModeEVecs[iPrevMode];

                    double hitNorm = 0;
                    int traceModeIndex = -1;
                    int modeCnt = defectK0s.Length;
                    for (int iMode = 0; iMode < modeCnt; iMode++)
                    {
                        if (traceModeIndexss.IndexOf(iMode) != -1)
                        {
                            // すでに対応付けている
                            continue;
                        }
                        double k0 = defectK0s[iMode];
                        System.Numerics.Complex[] eVec = defectEVecs[iMode];
                        double norm;
                        bool isHitSameMode = IsSameMode(prevEVec, k0, eVec, out norm);
                        if (isHitSameMode)
                        {
                            // より分布の近いモードを採用する
                            if (Math.Abs(norm - 1.0) < Math.Abs(hitNorm - 1.0))
                            {
                                // 追跡するモードのインデックス退避
                                traceModeIndex = iMode;
                                hitNorm = norm;
                            }
                        }
                    }
                    if (traceModeIndex != -1)
                    {
                        traceModeIndexss.Add(traceModeIndex);
                    }
                }
                // 並び替え
                {
                    int modeCnt = defectK0s.Length;
                    double[] prevK0s = new double[modeCnt];
                    defectK0s.CopyTo(prevK0s, 0);
                    System.Numerics.Complex[][] prevEVecs = new System.Numerics.Complex[modeCnt][];
                    defectEVecs.CopyTo(prevEVecs, 0);
                    int index = 0;
                    for (int i = 0; i < traceModeIndexss.Count; i++)
                    {
                        int traceModeIndex = traceModeIndexss[i];
                        defectK0s[index] = prevK0s[traceModeIndex];
                        defectEVecs[index] = prevEVecs[traceModeIndex];
                        index++;
                    }
                    for (int iMode = 0; iMode < modeCnt; iMode++)
                    {
                        if (!traceModeIndexss.Contains(iMode))
                        {
                            defectK0s[index] = prevK0s[iMode];
                            defectEVecs[index] = prevEVecs[iMode];
                            index++;
                        }
                    }
                    System.Diagnostics.Debug.Assert(index == modeCnt);
                }
                {
                    int modeCnt = defectK0s.Length;
                    WgPortInfo.PrevModeEVecs = new System.Numerics.Complex[modeCnt][];
                    defectEVecs.CopyTo(WgPortInfo.PrevModeEVecs, 0);
                }
            }
            else
            {
                int modeCnt = defectK0s.Length;
                WgPortInfo.PrevModeEVecs = new System.Numerics.Complex[modeCnt][];
                defectEVecs.CopyTo(WgPortInfo.PrevModeEVecs, 0);
            }
        }

        // 欠陥モード一覧の取得
        private void GetDefectModes(
            double minK0, double maxK0,
            double[] k0s, System.Numerics.Complex[][] eVecs,
            out double[] defectK0s, out System.Numerics.Complex[][] defectEVecs)
        {
            defectK0s = null;
            defectEVecs = null;
            IList<IList<int>> channelCoIdss = WgPortInfo.PCChannelCoIds;

            IList<double> retK0s = new List<double>();
            IList<System.Numerics.Complex[]> retEVecs = new List<System.Numerics.Complex[]>();
            int modeCnt = eVecs.Length;
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                double k0 = k0s[iMode];
                System.Numerics.Complex[] eVec = eVecs[iMode];
                bool isDefect = IsDefectMode(minK0, maxK0, channelCoIdss, k0, eVec);
                if (isDefect)
                {
                    retK0s.Add(k0);
                    System.Numerics.Complex[] copyEVec = new System.Numerics.Complex[eVec.Length];
                    eVec.CopyTo(copyEVec, 0);
                    retEVecs.Add(copyEVec);
                }
            }
            defectK0s = retK0s.ToArray();
            defectEVecs = retEVecs.ToArray();
        }

        // 欠陥モード？
        private bool IsDefectMode(
            double minK0, double maxK0,
            IList<IList<int>> channelCoIdss,
            double k0,
            System.Numerics.Complex[] eVec)
        {
            bool isHit = false;
            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            double latticeA = WgPortInfo.LatticeA;

            if (double.IsNaN(k0))
            {
                // 固有値がNaNのモードは除外する
                return isHit;
            }

            // フォトニック結晶導波路の導波モードを判定する
            //
            // 領域内の節点の界の絶対値の２乗の和を計算
            //   要素分割が均一であることが前提。面積を考慮しない。
            double totalPower = 0.0;
            for (int i = 0; i < nodeCnt; i++)
            {
                double fieldAbs = eVec[i].Magnitude;
                double power = fieldAbs * fieldAbs;
                totalPower += power;
            }

            // チャネルの座標IDリストを節点リストに変換
            int channelCnt = channelCoIdss.Count;
            IList<IList<uint>> channelNodess = new List<IList<uint>>();
            for (int channelIndex = 0; channelIndex < channelCnt; channelIndex++)
            {
                IList<int> channelCoIds = channelCoIdss[channelIndex];
                IList<uint> channelNodes = new List<uint>();
                channelNodess.Add(channelNodes);
                foreach (int coId in channelCoIds)
                {
                    int nodeId = World.PortCoord2Node(QuantityId, PortId, coId);
                    if (nodeId == -1)
                    {
                        continue;
                    }
                    channelNodes.Add((uint)nodeId);
                }
            }

            // チャンネル上の節点の界の絶対値の２乗の和を計算
            //   要素分割が均一であることが前提。面積を考慮しない。
            int channelNodeCnt = 0;
            double channelTotalPower = 0.0;
            for (int channelIndex = 0; channelIndex < channelCnt; channelIndex++)
            {
                IList<uint> channelNodes = channelNodess[channelIndex];

                foreach (uint nodeId in channelNodes)
                {
                    System.Numerics.Complex cvalue = eVec[nodeId];
                    double valAbs = cvalue.Magnitude;
                    double channelPower = valAbs * valAbs;
                    channelTotalPower += channelPower;
                    channelNodeCnt++;
                }
            }
            // 総和で比較する
            //const double powerRatioLimit = 0.5;
            const double powerRatioLimit = 0.4;
            if (Math.Abs(totalPower) >= Constants.PrecisionLowerLimit &&
                (channelTotalPower / totalPower) >= powerRatioLimit)
            {
                if (k0 < minK0 || k0 > maxK0)
                {
                    // skip
                    //System.Diagnostics.Debug.WriteLine(
                    //    "Skip defect mode a/λ = " + k0 * latticeA / (2.0 * Math.PI));
                }
                else
                {
                    System.Diagnostics.Debug.WriteLine(
                        "hit defect mode a/λ = " + k0 * latticeA / (2.0 * Math.PI));
                    isHit = true;
                }
            }
            return isHit;
        }

        // 同じモード？
        private bool IsSameMode(
            System.Numerics.Complex[] prevEVec,
            double k0,
            System.Numerics.Complex[] eVec,
            out double retNorm)
        {
            bool isHit = false;
            retNorm = 0;
            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            double latticeA = WgPortInfo.LatticeA;

            if (prevEVec == null)
            {
                return isHit;
            }

            System.Numerics.Complex[] workEVec1 = new System.Numerics.Complex[nodeCnt];
            prevEVec.CopyTo(workEVec1, 0);
            System.Numerics.Complex[] workEVec2 = new System.Numerics.Complex[nodeCnt];
            eVec.CopyTo(workEVec2, 0);
            System.Numerics.Complex norm1 =
                IvyFEM.Lapack.Functions.zdotu(IvyFEM.Lapack.Utils.Conjugate(workEVec1), workEVec1);
            norm1 = System.Numerics.Complex.Sqrt(norm1);
            System.Numerics.Complex norm2 =
                IvyFEM.Lapack.Functions.zdotu(IvyFEM.Lapack.Utils.Conjugate(workEVec2), workEVec2);
            norm2 = System.Numerics.Complex.Sqrt(norm2);
            workEVec1 = IvyFEM.Lapack.Functions.zscal(workEVec1, 1.0 / norm1.Magnitude);
            workEVec2 = IvyFEM.Lapack.Functions.zscal(workEVec2, 1.0 / norm2.Magnitude);
            System.Numerics.Complex norm12 =
                IvyFEM.Lapack.Functions.zdotu(IvyFEM.Lapack.Utils.Conjugate(workEVec1), workEVec2);
            double thLikeMin = 0.8;
            double thLikeMax = 1.0 + 1.0e-6;
            if (norm12.Magnitude >= thLikeMin && norm12.Magnitude < thLikeMax)
            {
                isHit = true;
                retNorm = norm12.Magnitude;
                System.Diagnostics.Debug.WriteLine(
                    "a/λ(eVal)= {0} norm (prev * current)= {1}", (latticeA / (2.0 * Math.PI) * k0), norm12.Magnitude);
            }
            return isHit;
        }

        private void SortEVals(double[] k0s, System.Numerics.Complex[][] eVecs)
        {
            int modeCnt = k0s.Length;
            var eValEVecs = new List<KeyValuePair<double, System.Numerics.Complex[]>>();
            for (int i = 0; i < modeCnt; i++)
            {
                eValEVecs.Add(
                    new KeyValuePair<double, System.Numerics.Complex[]>(k0s[i], eVecs[i]));
            }
            eValEVecs.Sort((a, b) =>
            {
                double k0A = a.Key;
                double k0B = b.Key;
                int cmp = 0;
                double th = 1.0e-12;
                double diff = k0A - k0B;
                if (diff > 0)
                {
                    cmp = 1;
                }
                else if (diff < 0)
                {
                    cmp = -1;
                }
                else
                {
                    cmp = 0;
                }
                return cmp;
            });

            for (int i = 0; i < modeCnt; i++)
            {
                k0s[i] = eValEVecs[i].Key;
                eVecs[i] = eValEVecs[i].Value;
            }
        }

        private void AdjustPhaseEVecs(System.Numerics.Complex[][] eVecs)
        {
            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            int modeCnt = eVecs.Length;
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var eVec = eVecs[iMode];
                int extendedNodeCnt = eVec.Length;
                System.Numerics.Complex maxValue = new System.Numerics.Complex(0, 0);
                double maxAbs = 0;
                for (int iNode = 0; iNode < nodeCnt; iNode++)
                {
                    System.Numerics.Complex value = eVec[iNode];
                    double abs = value.Magnitude;
                    if (abs > maxAbs)
                    {
                        maxAbs = abs;
                        maxValue = value;
                    }
                }
                System.Numerics.Complex phase = maxValue / maxAbs;

                for (int iNode = 0; iNode < extendedNodeCnt; iNode++)
                {
                    eVec[iNode] /= phase;
                }
            }
        }

        // フォトニックバンド固有値問題
        private void SolvePhotonicBand()
        {
            // input
            double betaX = BetaX;
            double betaY = BetaY;

            // 領域
            CalcMatrixs();

            double[] eVals;
            System.Numerics.Complex[][] eVecs;
            // photonic band (Φを直接解く)
            SolvePhotonicBandNonSVEAModeAsQuadraticGeneralizedEigen(betaX, betaY, out eVals, out eVecs);
            SortEVals(eVals, eVecs);
            AdjustPhaseEVecs(eVecs);

            Frequencys = new double[eVals.Length];
            for (int i = 0; i < eVals.Length; i++)
            {
                double k0 = eVals[i];
                double omega = k0 * Constants.C0;
                double freq = omega / (2.0 * Math.PI);
                Frequencys[i] = freq;
            }
            EVecs = eVecs;
        }

        // Φを直接解く
        private void SolvePhotonicBandNonSVEAModeAsQuadraticGeneralizedEigen(
            double betaX, double betaY,
            out double[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            eVals = null;
            eVecs = null;

            int nodeCnt = KMat0.RowLength;
            IList<int> bcNodes1 = WgPortInfo.BcNodess[0];
            IList<int> bcNodes2 = WgPortInfo.BcNodess[1];
            IList<int> bcNodes3 = WgPortInfo.BcNodess[2];
            IList<int> bcNodes4 = WgPortInfo.BcNodess[3];
            int bcNodeCnt1 = bcNodes1.Count;
            int bcNodeCnt3 = bcNodes3.Count;
            System.Diagnostics.Debug.Assert(WgPortInfo.BcNodess[0].Count == WgPortInfo.BcNodess[1].Count);
            System.Diagnostics.Debug.Assert(WgPortInfo.BcNodess[2].Count == WgPortInfo.BcNodess[3].Count);
            // コーナー点を求める(頂点1,2,3,4)
            int[] cornerNodes = new int[4];
            IList<int>[] corner1Bcs = { bcNodes1, bcNodes4 };
            IList<int>[] corner2Bcs = { bcNodes1, bcNodes3 };
            IList<int>[] corner3Bcs = { bcNodes3, bcNodes2 };
            IList<int>[] corner4Bcs = { bcNodes2, bcNodes4 };
            IList<int>[][] cornerBcss = { corner1Bcs, corner2Bcs, corner3Bcs, corner4Bcs};
            for (int cornerIndex = 0; cornerIndex < cornerNodes.Length; cornerIndex++)
            {
                IList<int>[] cornerBc = cornerBcss[cornerIndex];
                IList<int> workBcNodes1 = cornerBc[0];
                IList<int> workBcNodes2 = cornerBc[1];
                int cornerNodeId = -1;
                foreach (int nodeId in workBcNodes1)
                {
                    if (workBcNodes2.Contains(nodeId))
                    {
                        cornerNodeId = nodeId;
                        break;
                    }
                }
                System.Diagnostics.Debug.Assert(cornerNodeId != -1);
                cornerNodes[cornerIndex] = cornerNodeId;
            }
            int cornerNode1 = cornerNodes[0]; // 頂点1 
            int cornerNode2 = cornerNodes[1]; // 頂点2 
            int cornerNode3 = cornerNodes[2]; // 頂点3 
            int cornerNode4 = cornerNodes[3]; // 頂点4

            // 周期境界条件適用後の節点順を求める
            IList<int> nodes1 = new List<int>();
            // 境界1, 3
            uint[] bcIndexs = { 0, 2 };
            foreach (int bcIndex in bcIndexs)
            {
                IList<int> bcNodes = WgPortInfo.BcNodess[bcIndex];
                foreach (int nodeId in bcNodes)
                {
                    if (nodeId == cornerNode1 ||
                        nodeId == cornerNode3 ||
                        nodeId == cornerNode4)
                    {
                        // コーナー頂点2以外は除外する
                        continue;
                    }

                    if (nodes1.IndexOf(nodeId) == -1)
                    {
                        nodes1.Add(nodeId);
                    }
                }
            }
            for (int nodeId = 0; nodeId < nodeCnt; nodeId++)
            {
                if (bcNodes1.IndexOf(nodeId) != -1)
                {
                    continue;
                }
                if (bcNodes2.IndexOf(nodeId) != -1)
                {
                    continue;
                }
                if (bcNodes3.IndexOf(nodeId) != -1)
                {
                    continue;
                }
                if (bcNodes4.IndexOf(nodeId) != -1)
                {
                    continue;
                }

                if (nodes1.IndexOf(nodeId) == -1)
                {
                    // 内部領域
                    nodes1.Add(nodeId);
                }
            }
            int nodeCnt1 = nodes1.Count;
            bool isPortBc2Reverse = WgPortInfo.IsPortBc2Reverse;
            bool isPortBc4Reverse = WgPortInfo.IsPortBc4Reverse;
            double periodicDistanceX = WgPortInfo.PeriodicDistanceX;
            double periodicDistanceY = WgPortInfo.PeriodicDistanceY;

            var KMat = new IvyFEM.Lapack.ComplexMatrix(nodeCnt1, nodeCnt1);
            var MMat = new IvyFEM.Lapack.ComplexMatrix(nodeCnt1, nodeCnt1);

            System.Numerics.Complex expAX =
                System.Numerics.Complex.Exp(-1.0 * System.Numerics.Complex.ImaginaryOne * betaX * periodicDistanceX);
            System.Numerics.Complex expAY =
                System.Numerics.Complex.Exp(-1.0 * System.Numerics.Complex.ImaginaryOne * betaY * periodicDistanceY);

            if (WgPortInfo.IsAslantX)
            {
                // 斜め領域
                // 左下(コーナー頂点2)と左上(コーナー頂点1)のX座標オフセット
                int cornerCoId1 = World.PortNode2Coord(QuantityId, PortId, cornerNode1);
                int cornerCoId2 = World.PortNode2Coord(QuantityId, PortId, cornerNode2);
                double[] cornerCoord1 = World.GetCoord(QuantityId, cornerCoId1);
                double[] cornerCoord2 = World.GetCoord(QuantityId, cornerCoId2);
                double ofsX = cornerCoord1[0] - cornerCoord2[0];
                if (Math.Abs(ofsX) >= Constants.PrecisionLowerLimit)
                {
                    // 斜め領域の場合
                    //  Y方向の周期境界条件にはオフセット分のX方向成分の因子が入る
                    System.Diagnostics.Debug.WriteLine("ofsX: {0}", ofsX);
                    expAY *= System.Numerics.Complex.Exp(-1.0 * System.Numerics.Complex.ImaginaryOne * betaX * ofsX);
                }
            }

            //////////////////////////////////////////////////////
            // 独立節点
            for (int ii = 0; ii < nodeCnt1; ii++)
            {
                int i = nodes1[ii];
                // 独立節点
                for (int jj = 0; jj < nodeCnt1; jj++)
                {
                    int j = nodes1[jj];
                    KMat[ii, jj] += KMat0[i, j];
                    MMat[ii, jj] += MMat0[i, j];
                }
                // 境界2 ->境界1
                for (int j = 0; j < bcNodeCnt1; j++)
                {
                    int jB1 = bcNodes1[j];
                    int jB2 = bcNodes2[isPortBc2Reverse ? (bcNodeCnt1 - 1 - j) : j];
                    int jjB1 = nodes1.IndexOf(jB1);
                    if (jB1 == cornerNode1 || jB1 == cornerNode2 || jB1 == cornerNode3 || jB1 == cornerNode4)
                    {
                        continue;
                    }

                    // [Ki2] ->[Ki1]
                    KMat[ii, jjB1] += expAX * KMat0[i, jB2];
                    MMat[ii, jjB1] += expAX * MMat0[i, jB2];
                }
                // 境界4 ->境界3
                for (int j = 0; j < bcNodeCnt3; j++)
                {
                    int jB3 = bcNodes3[j];
                    int jB4 = bcNodes4[isPortBc4Reverse ? (bcNodeCnt3 - 1 - j) : j];
                    int jjB3 = nodes1.IndexOf(jB3);
                    if (jB3 == cornerNode1 || jB3 == cornerNode2 || jB3 == cornerNode3 || jB3 == cornerNode4)
                    {
                        continue;
                    }
                    // [Ki4] ->[Ki3]
                    KMat[ii, jjB3] += expAY * KMat0[i, jB4];
                    MMat[ii, jjB3] += expAY * MMat0[i, jB4];
                }
                // コーナー頂点1 -> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC1 = cornerNode1;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic1] ->[Kic2]
                    KMat[ii, jjC2] += expAY * KMat0[i, jC1];
                    MMat[ii, jjC2] += expAY * MMat0[i, jC1];
                }
                // コーナー頂点3 -> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC3 = cornerNode3;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic3] ->[Kic2]
                    KMat[ii, jjC2] += expAX * KMat0[i, jC3];
                    MMat[ii, jjC2] += expAX * MMat0[i, jC3];
                }
                // コーナー頂点4-> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC4 = cornerNode4;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic4] ->[Kic2]
                    KMat[ii, jjC2] += expAX * expAY * KMat0[i, jC4];
                    MMat[ii, jjC2] += expAX * expAY * MMat0[i, jC4];
                }
            }
            //////////////////////////////////////////////////////
            // 境界2 ->境界1
            for (int i = 0; i < bcNodeCnt1; i++)
            {
                int iB1 = bcNodes1[i];
                int iB2 = bcNodes2[isPortBc2Reverse ? (bcNodeCnt1 - 1 - i) : i];
                int iiB1 = nodes1.IndexOf(iB1);
                if (iB1 == cornerNode1 || iB1 == cornerNode2 || iB1 == cornerNode3 || iB1 == cornerNode4)
                {
                    continue;
                }
                System.Numerics.Complex phaseI = (1.0 / expAX);
                // 独立節点
                for (int jj = 0; jj < nodeCnt1; jj++)
                {
                    int j = nodes1[jj];
                    KMat[iiB1, jj] += phaseI * KMat0[iB2, j];
                    MMat[iiB1, jj] += phaseI * MMat0[iB2, j];
                }
                // 境界2 ->境界1
                for (int j = 0; j < bcNodeCnt1; j++)
                {
                    int jB1 = bcNodes1[j];
                    int jB2 = bcNodes2[isPortBc2Reverse ? (bcNodeCnt1 - 1 - j) : j];
                    int jjB1 = nodes1.IndexOf(jB1);
                    if (jB1 == cornerNode1 || jB1 == cornerNode2 || jB1 == cornerNode3 || jB1 == cornerNode4)
                    {
                        continue;
                    }

                    // [Ki2] ->[Ki1]
                    KMat[iiB1, jjB1] += phaseI * expAX * KMat0[iB2, jB2];
                    MMat[iiB1, jjB1] += phaseI * expAX * MMat0[iB2, jB2];
                }
                // 境界4 ->境界3
                for (int j = 0; j < bcNodeCnt3; j++)
                {
                    int jB3 = bcNodes3[j];
                    int jB4 = bcNodes4[isPortBc4Reverse ? (bcNodeCnt3 - 1 - j) : j];
                    int jjB3 = nodes1.IndexOf(jB3);
                    if (jB3 == cornerNode1 || jB3 == cornerNode2 || jB3 == cornerNode3 || jB3 == cornerNode4)
                    {
                        continue;
                    }
                    // [Ki4] ->[Ki3]
                    KMat[iiB1, jjB3] += phaseI * expAY * KMat0[iB2, jB4];
                    MMat[iiB1, jjB3] += phaseI * expAY * MMat0[iB2, jB4];
                }
                // コーナー頂点1 -> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC1 = cornerNode1;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic1] ->[Kic2]
                    KMat[iiB1, jjC2] += phaseI * expAY * KMat0[iB2, jC1];
                    MMat[iiB1, jjC2] += phaseI * expAY * MMat0[iB2, jC1];
                }
                // コーナー頂点3 -> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC3 = cornerNode3;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic3] ->[Kic2]
                    KMat[iiB1, jjC2] += phaseI * expAX * KMat0[iB2, jC3];
                    MMat[iiB1, jjC2] += phaseI * expAX * MMat0[iB2, jC3];
                }
                // コーナー頂点4-> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC4 = cornerNode4;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic4] ->[Kic2]
                    KMat[iiB1, jjC2] += phaseI * expAX * expAY * KMat0[iB2, jC4];
                    MMat[iiB1, jjC2] += phaseI * expAX * expAY * MMat0[iB2, jC4];
                }
            }
            //////////////////////////////////////////////////////
            // 境界4 ->境界3
            for (int i = 0; i < bcNodeCnt3; i++)
            {
                int iB3 = bcNodes3[i];
                int iB4 = bcNodes4[isPortBc4Reverse ? (bcNodeCnt3 - 1 - i) : i];
                int iiB3 = nodes1.IndexOf(iB3);
                if (iB3 == cornerNode1 || iB3 == cornerNode2 || iB3 == cornerNode3 || iB3 == cornerNode4)
                {
                    continue;
                }
                System.Numerics.Complex phaseI = (1.0 / expAY);
                // 独立節点
                for (int jj = 0; jj < nodeCnt1; jj++)
                {
                    int j = nodes1[jj];
                    KMat[iiB3, jj] += phaseI * KMat0[iB4, j];
                    MMat[iiB3, jj] += phaseI * MMat0[iB4, j];
                }
                // 境界2 ->境界1
                for (int j = 0; j < bcNodeCnt1; j++)
                {
                    int jB1 = bcNodes1[j];
                    int jB2 = bcNodes2[isPortBc2Reverse ? (bcNodeCnt1 - 1 - j) : j];
                    int jjB1 = nodes1.IndexOf(jB1);
                    if (jB1 == cornerNode1 || jB1 == cornerNode2 || jB1 == cornerNode3 || jB1 == cornerNode4)
                    {
                        continue;
                    }

                    // [Ki2] ->[Ki1]
                    KMat[iiB3, jjB1] += phaseI * expAX * KMat0[iB4, jB2];
                    MMat[iiB3, jjB1] += phaseI * expAX * MMat0[iB4, jB2];
                }
                // 境界4 ->境界3
                for (int j = 0; j < bcNodeCnt3; j++)
                {
                    int jB3 = bcNodes3[j];
                    int jB4 = bcNodes4[isPortBc4Reverse ? (bcNodeCnt3 - 1 - j) : j];
                    int jjB3 = nodes1.IndexOf(jB3);
                    if (jB3 == cornerNode1 || jB3 == cornerNode2 || jB3 == cornerNode3 || jB3 == cornerNode4)
                    {
                        continue;
                    }
                    // [Ki4] ->[Ki3]
                    KMat[iiB3, jjB3] += phaseI * expAY * KMat0[iB4, jB4];
                    MMat[iiB3, jjB3] += phaseI * expAY * MMat0[iB4, jB4];
                }
                // コーナー頂点1 -> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC1 = cornerNode1;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic1] ->[Kic2]
                    KMat[iiB3, jjC2] += phaseI * expAY * KMat0[iB4, jC1];
                    MMat[iiB3, jjC2] += phaseI * expAY * MMat0[iB4, jC1];
                }
                // コーナー頂点3 -> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC3 = cornerNode3;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic3] ->[Kic2]
                    KMat[iiB3, jjC2] += phaseI * expAX * KMat0[iB4, jC3];
                    MMat[iiB3, jjC2] += phaseI * expAX * MMat0[iB4, jC3];
                }
                // コーナー頂点4-> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC4 = cornerNode4;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic4] ->[Kic2]
                    KMat[iiB3, jjC2] += phaseI * expAX * expAY * KMat0[iB4, jC4];
                    MMat[iiB3, jjC2] += phaseI * expAX * expAY * MMat0[iB4, jC4];
                }
            }
            //////////////////////////////////////////////////////
            // コーナー頂点1 -> コーナー頂点2
            {
                int iC2 = cornerNode2;
                int iC1 = cornerNode1;
                int iiC2 = nodes1.IndexOf(iC2);
                System.Numerics.Complex phaseI = (1.0 / expAY);
                // 独立節点
                for (int jj = 0; jj < nodeCnt1; jj++)
                {
                    int j = nodes1[jj];
                    KMat[iiC2, jj] += phaseI * KMat0[iC1, j];
                    MMat[iiC2, jj] += phaseI * MMat0[iC1, j];
                }
                // 境界2 ->境界1
                for (int j = 0; j < bcNodeCnt1; j++)
                {
                    int jB1 = bcNodes1[j];
                    int jB2 = bcNodes2[isPortBc2Reverse ? (bcNodeCnt1 - 1 - j) : j];
                    int jjB1 = nodes1.IndexOf(jB1);
                    if (jB1 == cornerNode1 || jB1 == cornerNode2 || jB1 == cornerNode3 || jB1 == cornerNode4)
                    {
                        continue;
                    }

                    // [Ki2] ->[Ki1]
                    KMat[iiC2, jjB1] += phaseI * expAX * KMat0[iC1, jB2];
                    MMat[iiC2, jjB1] += phaseI * expAX * MMat0[iC1, jB2];
                }
                // 境界4 ->境界3
                for (int j = 0; j < bcNodeCnt3; j++)
                {
                    int jB3 = bcNodes3[j];
                    int jB4 = bcNodes4[isPortBc4Reverse ? (bcNodeCnt3 - 1 - j) : j];
                    int jjB3 = nodes1.IndexOf(jB3);
                    if (jB3 == cornerNode1 || jB3 == cornerNode2 || jB3 == cornerNode3 || jB3 == cornerNode4)
                    {
                        continue;
                    }
                    // [Ki4] ->[Ki3]
                    KMat[iiC2, jjB3] += phaseI * expAY * KMat0[iC1, jB4];
                    MMat[iiC2, jjB3] += phaseI * expAY * MMat0[iC1, jB4];
                }
                // コーナー頂点1 -> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC1 = cornerNode1;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic1] ->[Kic2]
                    KMat[iiC2, jjC2] += phaseI * expAY * KMat0[iC1, jC1];
                    MMat[iiC2, jjC2] += phaseI * expAY * MMat0[iC1, jC1];
                }
                // コーナー頂点3 -> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC3 = cornerNode3;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic3] ->[Kic2]
                    KMat[iiC2, jjC2] += phaseI * expAX * KMat0[iC1, jC3];
                    MMat[iiC2, jjC2] += phaseI * expAX * MMat0[iC1, jC3];
                }
                // コーナー頂点4-> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC4 = cornerNode4;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic4] ->[Kic2]
                    KMat[iiC2, jjC2] += phaseI * expAX * expAY * KMat0[iC1, jC4];
                    MMat[iiC2, jjC2] += phaseI * expAX * expAY * MMat0[iC1, jC4];
                }
            }
            //////////////////////////////////////////////////////
            // コーナー頂点3 -> コーナー頂点2
            {
                int iC2 = cornerNode2;
                int iC3 = cornerNode3;
                int iiC2 = nodes1.IndexOf(iC2);
                System.Numerics.Complex phaseI = (1.0 / expAX);
                // 独立節点
                for (int jj = 0; jj < nodeCnt1; jj++)
                {
                    int j = nodes1[jj];
                    KMat[iiC2, jj] += phaseI * KMat0[iC3, j];
                    MMat[iiC2, jj] += phaseI * MMat0[iC3, j];
                }
                // 境界2 ->境界1
                for (int j = 0; j < bcNodeCnt1; j++)
                {
                    int jB1 = bcNodes1[j];
                    int jB2 = bcNodes2[isPortBc2Reverse ? (bcNodeCnt1 - 1 - j) : j];
                    int jjB1 = nodes1.IndexOf(jB1);
                    if (jB1 == cornerNode1 || jB1 == cornerNode2 || jB1 == cornerNode3 || jB1 == cornerNode4)
                    {
                        continue;
                    }

                    // [Ki2] ->[Ki1]
                    KMat[iiC2, jjB1] += phaseI * expAX * KMat0[iC3, jB2];
                    MMat[iiC2, jjB1] += phaseI * expAX * MMat0[iC3, jB2];
                }
                // 境界4 ->境界3
                for (int j = 0; j < bcNodeCnt3; j++)
                {
                    int jB3 = bcNodes3[j];
                    int jB4 = bcNodes4[isPortBc4Reverse ? (bcNodeCnt3 - 1 - j) : j];
                    int jjB3 = nodes1.IndexOf(jB3);
                    if (jB3 == cornerNode1 || jB3 == cornerNode2 || jB3 == cornerNode3 || jB3 == cornerNode4)
                    {
                        continue;
                    }
                    // [Ki4] ->[Ki3]
                    KMat[iiC2, jjB3] += phaseI * expAY * KMat0[iC3, jB4];
                    MMat[iiC2, jjB3] += phaseI * expAY * MMat0[iC3, jB4];
                }
                // コーナー頂点1 -> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC1 = cornerNode1;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic1] ->[Kic2]
                    KMat[iiC2, jjC2] += phaseI * expAY * KMat0[iC3, jC1];
                    MMat[iiC2, jjC2] += phaseI * expAY * MMat0[iC3, jC1];
                }
                // コーナー頂点3 -> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC3 = cornerNode3;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic3] ->[Kic2]
                    KMat[iiC2, jjC2] += phaseI * expAX * KMat0[iC3, jC3];
                    MMat[iiC2, jjC2] += phaseI * expAX * MMat0[iC3, jC3];
                }
                // コーナー頂点4-> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC4 = cornerNode4;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic4] ->[Kic2]
                    KMat[iiC2, jjC2] += phaseI * expAX * expAY * KMat0[iC3, jC4];
                    MMat[iiC2, jjC2] += phaseI * expAX * expAY * MMat0[iC3, jC4];
                }
            }
            //////////////////////////////////////////////////////
            // コーナー頂点4-> コーナー頂点2
            {
                int iC2 = cornerNode2;
                int iC4 = cornerNode4;
                int iiC2 = nodes1.IndexOf(iC2);
                System.Numerics.Complex phaseI = (1.0 / expAX) * (1.0 /expAY);
                // 独立節点
                for (int jj = 0; jj < nodeCnt1; jj++)
                {
                    int j = nodes1[jj];
                    KMat[iiC2, jj] += phaseI * KMat0[iC4, j];
                    MMat[iiC2, jj] += phaseI * MMat0[iC4, j];
                }
                // 境界2 ->境界1
                for (int j = 0; j < bcNodeCnt1; j++)
                {
                    int jB1 = bcNodes1[j];
                    int jB2 = bcNodes2[isPortBc2Reverse ? (bcNodeCnt1 - 1 - j) : j];
                    int jjB1 = nodes1.IndexOf(jB1);
                    if (jB1 == cornerNode1 || jB1 == cornerNode2 || jB1 == cornerNode3 || jB1 == cornerNode4)
                    {
                        continue;
                    }

                    // [Ki2] ->[Ki1]
                    KMat[iiC2, jjB1] += phaseI * expAX * KMat0[iC4, jB2];
                    MMat[iiC2, jjB1] += phaseI * expAX * MMat0[iC4, jB2];
                }
                // 境界4 ->境界3
                for (int j = 0; j < bcNodeCnt3; j++)
                {
                    int jB3 = bcNodes3[j];
                    int jB4 = bcNodes4[isPortBc4Reverse ? (bcNodeCnt3 - 1 - j) : j];
                    int jjB3 = nodes1.IndexOf(jB3);
                    if (jB3 == cornerNode1 || jB3 == cornerNode2 || jB3 == cornerNode3 || jB3 == cornerNode4)
                    {
                        continue;
                    }
                    // [Ki4] ->[Ki3]
                    KMat[iiC2, jjB3] += phaseI * expAY * KMat0[iC4, jB4];
                    MMat[iiC2, jjB3] += phaseI * expAY * MMat0[iC4, jB4];
                }
                // コーナー頂点1 -> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC1 = cornerNode1;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic1] ->[Kic2]
                    KMat[iiC2, jjC2] += phaseI * expAY * KMat0[iC4, jC1];
                    MMat[iiC2, jjC2] += phaseI * expAY * MMat0[iC4, jC1];
                }
                // コーナー頂点3 -> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC3 = cornerNode3;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic3] ->[Kic2]
                    KMat[iiC2, jjC2] += phaseI * expAX * KMat0[iC4, jC3];
                    MMat[iiC2, jjC2] += phaseI * expAX * MMat0[iC4, jC3];
                }
                // コーナー頂点4-> コーナー頂点2
                {
                    int jC2 = cornerNode2;
                    int jC4 = cornerNode4;
                    int jjC2 = nodes1.IndexOf(jC2);
                    // [Kic4] ->[Kic2]
                    KMat[iiC2, jjC2] += phaseI * expAX * expAY * KMat0[iC4, jC4];
                    MMat[iiC2, jjC2] += phaseI * expAX * expAY * MMat0[iC4, jC4];
                }
            }

            System.Diagnostics.Debug.WriteLine("solve eigen");
            double[] eVals0;
            System.Numerics.Complex[][] eVecs0;
            SolveComplexHermiteBandGeneralizedEigen(KMat, MMat, out eVals0, out eVecs0);

            int modeCnt = eVals0.Length;
            eVecs = new System.Numerics.Complex[modeCnt][];
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                System.Numerics.Complex[] eVec = new System.Numerics.Complex[nodeCnt];
                eVecs[iMode] = eVec;

                System.Numerics.Complex[] eVec0 = eVecs0[iMode];
                System.Diagnostics.Debug.Assert(eVec0.Length == nodeCnt1);

                // 独立節点
                for (int ii = 0; ii < nodeCnt1; ii++)
                {
                    int i = nodes1[ii];
                    eVec[i] = eVec0[ii];
                }
                // 境界1 ->境界2
                for (int i = 0; i < bcNodeCnt1; i++)
                {
                    int iB1 = bcNodes1[i];
                    int iB2 = bcNodes2[isPortBc2Reverse ? (bcNodeCnt1 - 1 - i) : i];
                    int iiB1 = nodes1.IndexOf(iB1);
                    if (iB1 == cornerNode1 || iB1 == cornerNode2 || iB1 == cornerNode3 || iB1 == cornerNode4)
                    {
                        continue;
                    }
                    eVec[iB2] = expAX * eVec0[iiB1];
                }
                // 境界3 ->境界4
                for (int i = 0; i < bcNodeCnt3; i++)
                {
                    int iB3 = bcNodes3[i];
                    int iB4 = bcNodes4[isPortBc4Reverse ? (bcNodeCnt3 - 1 - i) : i];
                    int iiB3 = nodes1.IndexOf(iB3);
                    if (iB3 == cornerNode1 || iB3 == cornerNode2 || iB3 == cornerNode3 || iB3 == cornerNode4)
                    {
                        continue;
                    }
                    eVec[iB4] = expAY * eVec0[iiB3];
                }
                // コーナー頂点2 -> コーナー頂点1
                {
                    int iC2 = cornerNode2;
                    int iC1 = cornerNode1;
                    int iiC2 = nodes1.IndexOf(iC2);
                    eVec[iC1] = expAY * eVec0[iiC2];
                }
                // コーナー頂点2 -> コーナー頂点3
                {
                    int iC2 = cornerNode2;
                    int iC3 = cornerNode3;
                    int iiC2 = nodes1.IndexOf(iC2);
                    eVec[iC3] = expAX * eVec0[iiC2];
                }
                // コーナー頂点2-> コーナー頂点4
                {
                    int iC2 = cornerNode2;
                    int iC4 = cornerNode4;
                    int iiC2 = nodes1.IndexOf(iC2);
                    eVec[iC4] = expAX * expAY * eVec0[iiC2];
                }
            }
            // k0に変換
            {
                eVals = new double[modeCnt];
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    // λ = k0^2と置いた場合( k0 = sqrt(λ) )
                    eVals[iMode] = Math.Sqrt(eVals0[iMode]);
                }
            }
        }

        // フォトニックバンドギャップを計算する
        public static void GetPBG(
            IList<double[]> frequencyss,
            double minFreq,
            double maxFreq,
            out double gapMinFreq,
            out double gapMaxFreq)
        {
            gapMinFreq = 0.0;
            gapMaxFreq = 0.0;
            List<double> freqs = new List<double>();
            foreach (double[] values in frequencyss)
            {
                foreach (double value in values)
                {
                    freqs.Add(value);
                }
            }
            freqs.Sort();
            double maxGap = double.MinValue;
            for (int i = 1; i < freqs.Count; i++)
            {
                if (freqs[i - 1] < minFreq)
                {
                    continue;
                }
                if (freqs[i] > maxFreq)
                {
                    break;
                }
                double gap = freqs[i] - freqs[i - 1];
                if (gap > maxGap)
                {
                    gapMinFreq = freqs[i - 1];
                    gapMaxFreq = freqs[i];
                    maxGap = gap;
                }
            }
        }
    }
}