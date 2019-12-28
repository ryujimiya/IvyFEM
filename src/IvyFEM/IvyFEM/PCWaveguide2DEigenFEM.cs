using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class PCWaveguide2DEigenFEM : FEM
    {
        public uint QuantityId { get; private set; } = 0;
        public uint PortId { get; private set; } = 0;
        public bool IsTMMode { get; set; } = false;
        public PCWaveguidePortInfo WgPortInfo { get; private set; }
        private IvyFEM.Lapack.DoubleMatrix KMat = null;
        private IvyFEM.Lapack.DoubleMatrix CMat = null;
        private IvyFEM.Lapack.DoubleMatrix MMat = null;

        private IvyFEM.Lapack.DoubleMatrix KMat0 = null;
        private IvyFEM.Lapack.DoubleMatrix CMat0 = null;
        private IvyFEM.Lapack.DoubleMatrix MMat0 = null;

        // 境界1積分
        public IvyFEM.Lapack.DoubleMatrix TxxB1 { get; private set; } = null;
        public IvyFEM.Lapack.DoubleMatrix RyyB1 { get; private set; } = null;
        public IvyFEM.Lapack.DoubleMatrix UzzB1 { get; private set; } = null;

        // Solve
        // Input
        public double Frequency { get; set; }
        // Output
        public System.Numerics.Complex[] Betas { get; private set; }
        public System.Numerics.Complex[][] EVecs { get; private set; }
        public System.Numerics.Complex[][] FxEVecs { get; private set; }
        public System.Numerics.Complex[][] BcEVecs { get; private set; }
        public System.Numerics.Complex[][] BcFxEVecs { get; private set; }

        public PCWaveguide2DEigenFEM(
            FEWorld world, uint quantityId, uint portId, PCWaveguidePortInfo wgPortInfo)
        {
            World = world;
            QuantityId = quantityId;
            PortId = portId;
            WgPortInfo = wgPortInfo;
        }

        private void CalcMatrixs(double k0)
        {
            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            IList<uint> feIds = World.GetPeriodicPortTriangleFEIds(QuantityId, PortId);

            KMat0 = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            CMat0 = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
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

                        // 要素剛性行列
                        //  要素質量行列を正定値行列にするために剛性行列、結合行列の符号を反転する
                        double kVal = -maPxx * sNyNy[row, col] - maPyy * sNxNx[row, col]
                                             + k0 * k0 * maQzz * sNN[row, col];
                        double cVal = 0;
                        if (WgPortInfo.IsYDirectionPeriodic)
                        {
                            cVal = -maPyy * (sNyN[row, col] - sNyN[col, row]);
                        }
                        else
                        {
                            cVal = -maPyy * (sNxN[row, col] - sNxN[col, row]);
                        }
                        // 要素質量行列
                        double mVal = maPyy * sNN[row, col];

                        KMat0[rowNodeId, colNodeId] += kVal;
                        CMat0[rowNodeId, colNodeId] += cVal;
                        MMat0[rowNodeId, colNodeId] += mVal;
                    }
                }
            }

            /////////////////////////////////////////////////////
            // 周期構造条件の適用
            int bcNodeCnt = WgPortInfo.BcNodess[0].Count;
            System.Diagnostics.Debug.Assert(WgPortInfo.BcNodess[0].Count == WgPortInfo.BcNodess[1].Count);
            int nodeCnt1 = nodeCnt - bcNodeCnt; // 境界2を除く

            bool isSVEA = WgPortInfo.IsSVEA;
            bool isPortBc2Reverse = WgPortInfo.IsPortBc2Reverse;
            if (isSVEA)
            {
                KMat = new IvyFEM.Lapack.DoubleMatrix(nodeCnt1, nodeCnt1);
                CMat = new IvyFEM.Lapack.DoubleMatrix(nodeCnt1, nodeCnt1);
                MMat = new IvyFEM.Lapack.DoubleMatrix(nodeCnt1, nodeCnt1);
                // 境界1
                for (int i = 0; i < bcNodeCnt; i++)
                {
                    for (int j = 0; j < bcNodeCnt; j++)
                    {
                        KMat[i, j] = KMat0[i, j];
                        CMat[i, j] = CMat0[i, j];
                        MMat[i, j] = MMat0[i, j];
                    }
                    for (int j = bcNodeCnt * 2; j < nodeCnt; j++)
                    {
                        int j1 = j - bcNodeCnt;
                        KMat[i, j1] = KMat0[i, j];
                        CMat[i, j1] = CMat0[i, j];
                        MMat[i, j1] = MMat0[i, j];
                    }
                }
                // 内部
                for (int i = bcNodeCnt * 2; i < nodeCnt; i++)
                {
                    int i1 = i - bcNodeCnt;
                    for (int j = 0; j < bcNodeCnt; j++)
                    {
                        KMat[i1, j] = KMat0[i, j];
                        CMat[i1, j] = CMat0[i, j];
                        MMat[i1, j] = MMat0[i, j];
                    }
                    for (int j = bcNodeCnt * 2; j < nodeCnt; j++)
                    {
                        int j1 = j - bcNodeCnt;
                        KMat[i1, j1] = KMat0[i, j];
                        CMat[i1, j1] = CMat0[i, j];
                        MMat[i1, j1] = MMat0[i, j];
                    }
                }
                // 境界条件
                for (int i = 0; i < bcNodeCnt; i++)
                {
                    for (int j = bcNodeCnt; j < bcNodeCnt * 2; j++)
                    {
                        int j1 = isPortBc2Reverse ?
                            (bcNodeCnt * 2 - 1 - j) :
                            (j - bcNodeCnt);
                        KMat[i, j1] += KMat0[i, j];
                        CMat[i, j1] += CMat0[i, j];
                        MMat[i, j1] += MMat0[i, j];
                    }
                }
                for (int i = bcNodeCnt * 2; i < nodeCnt; i++)
                {
                    int i1 = i - bcNodeCnt;
                    for (int j = bcNodeCnt; j < bcNodeCnt * 2; j++)
                    {
                        int j1 = isPortBc2Reverse ?
                            (bcNodeCnt * 2 - 1 - j) :
                            (j - bcNodeCnt);
                        KMat[i1, j1] += KMat0[i, j];
                        CMat[i1, j1] += CMat0[i, j];
                        MMat[i1, j1] += MMat0[i, j];
                    }
                }
                for (int i = bcNodeCnt; i < bcNodeCnt * 2; i++)
                {
                    int i1 = isPortBc2Reverse ?
                        (bcNodeCnt * 2 - 1 - i) :
                        (i - bcNodeCnt);
                    for (int j = 0; j < bcNodeCnt; j++)
                    {
                        KMat[i1, j] += KMat0[i, j];
                        CMat[i1, j] += CMat0[i, j];
                        MMat[i1, j] += MMat0[i, j];
                    }
                    for (int j = bcNodeCnt * 2; j < nodeCnt; j++)
                    {
                        int j1 = j - bcNodeCnt;
                        KMat[i1, j1] += KMat0[i, j];
                        CMat[i1, j1] += CMat0[i, j];
                        MMat[i1, j1] += MMat0[i, j];
                    }
                    for (int j = bcNodeCnt; j < bcNodeCnt * 2; j++)
                    {
                        int j1 = isPortBc2Reverse ?
                            (bcNodeCnt * 2 - 1 - j) :
                            (j - bcNodeCnt);
                        KMat[i1, j1] += KMat0[i, j];
                        CMat[i1, j1] += CMat0[i, j];
                        MMat[i1, j1] += MMat0[i, j];
                    }
                }
                // 行列要素check
                {
                    double th = 1.0e-12;
                    for (int i = 0; i < nodeCnt1; i++)
                    {
                        for (int j = i; j < nodeCnt1; j++)
                        {
                            // [K]は対称行列
                            System.Diagnostics.Debug.Assert(
                                Math.Abs(KMat[i, j] - KMat[j, i]) < th);
                            // [M]は対称行列
                            System.Diagnostics.Debug.Assert(
                                Math.Abs(MMat[i, j] - MMat[j, i]) < th);
                            // [C]は反対称行列
                            System.Diagnostics.Debug.Assert(
                                Math.Abs((-CMat[i, j]) - CMat[j, i]) < th);
                        }
                    }
                }
            }
            else
            {
                // KMatは使わない
                KMat = null;
                MMat = null;
                CMat = null;
            }
        }

        private void CalcB1Matrix()
        {
            // 境界1
            uint bcIndex = 0;
            IList<uint> feIds = World.GetPeriodicPortLineFEIds(QuantityId, PortId, bcIndex);

            IList<int> bcNodes = WgPortInfo.BcNodess[(int)bcIndex];
            int bcNodeCnt = bcNodes.Count;

            TxxB1 = new IvyFEM.Lapack.DoubleMatrix(bcNodeCnt, bcNodeCnt);
            RyyB1 = new IvyFEM.Lapack.DoubleMatrix(bcNodeCnt, bcNodeCnt);
            UzzB1 = new IvyFEM.Lapack.DoubleMatrix(bcNodeCnt, bcNodeCnt);

            foreach (uint feId in feIds)
            {
                LineFE lineFE = World.GetLineFE(QuantityId, feId);
                uint elemNodeCnt = lineFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = lineFE.NodeCoordIds[iNode];
                    int nodeId = World.PortCoord2Node(QuantityId, PortId, coId);
                    nodes[iNode] = nodeId;
                }

                Material ma0 = World.GetMaterial(lineFE.MaterialId);
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

                double[,] sNN = lineFE.CalcSNN();
                double[,] sNyNy = lineFE.CalcSNxNx();
                for (int row = 0; row < elemNodeCnt; row++)
                {
                    int rowNodeId = nodes[row];
                    if (rowNodeId == -1)
                    {
                        continue;
                    }
                    int rowNodeIdB1 = bcNodes.IndexOf(rowNodeId);
                    for (int col = 0; col < elemNodeCnt; col++)
                    {
                        int colNodeId = nodes[col];
                        if (colNodeId == -1)
                        {
                            continue;
                        }
                        int colNodeIdB1 = bcNodes.IndexOf(colNodeId);
                        double txxVal = maPxx * sNyNy[row, col];
                        double ryyVal = maPyy * sNN[row, col];
                        double uzzVal = maQzz * sNN[row, col];

                        TxxB1[rowNodeIdB1, colNodeIdB1] += txxVal;
                        RyyB1[rowNodeIdB1, colNodeIdB1] += ryyVal;
                        UzzB1[rowNodeIdB1, colNodeIdB1] += uzzVal;
                    }
                }
            }
        }

        public override void Solve()
        {
            Betas = null;
            EVecs = null;
            FxEVecs = null;
            BcEVecs = null;
            BcFxEVecs = null;

            // 角周波数
            double omega = 2.0 * Math.PI * Frequency;
            // 波数
            double k0 = omega / Constants.C0;

            double periodicDistance = WgPortInfo.PeriodicDistanceX;

            bool isYDirectionPeriodic = WgPortInfo.IsYDirectionPeriodic;
            if (isYDirectionPeriodic)
            {
                // 回転移動させない
                World.RotAngle = 0;
                World.RotOrigin = null;
            }
            else
            {
                // 回転移動
                double rotAngle;
                double[] rotOrigin;
                PCWaveguideUtils.GetRotOriginRotAngleFromY(
                    World, WgPortInfo.BcEdgeIds1.ToArray(), out rotAngle, out rotOrigin);

                World.RotAngle = rotAngle;
                World.RotOrigin = rotOrigin;
            }

            // 領域
            CalcMatrixs(k0);

            // 境界1
            CalcB1Matrix();

            double minEffN = WgPortInfo.MinEffN;
            double maxEffN = WgPortInfo.MaxEffN;
            double minBeta = minEffN;
            double maxBeta = maxEffN;
            // PC導波路の場合は、波数が[0, π]の領域から探索する
            {
                double minWaveNum = WgPortInfo.MinWaveNum;
                double maxWaveNum = WgPortInfo.MaxWaveNum;
                double minBetaBZ = minWaveNum * (2.0 * Math.PI / periodicDistance) / k0;
                double maxBetaBZ = maxWaveNum * (2.0 * Math.PI / periodicDistance) / k0;
                if (minBetaBZ > minBeta)
                {
                    minBeta = minBetaBZ;
                }
                if (maxBetaBZ < maxBeta)
                {
                    maxBeta = maxBetaBZ;
                }
                //System.Diagnostics.Debug.WriteLine("minWaveNum:{0}, maxWaveNum: {1}", minWaveNum, maxWaveNum);
                //System.Diagnostics.Debug.WriteLine("minBeta: {0}, maxBeta: {1}", minBeta, maxBeta);
            }

            bool isSVEA = WgPortInfo.IsSVEA;
            System.Numerics.Complex[] eVals;
            System.Numerics.Complex[][] eVecs;
            if (isSVEA)
            {
                SolveAsStandardEigenWithRealMat(out eVals, out eVecs);
                //SolveAsGeneralizedEigenWithRealMat(out eVals, out eVecs);
                //SolveIterAsGeneralizedEigen(k0, minBeta, maxBeta, out eVals, out eVecs);
            }
            else
            {
                SolveNonSVEAModeAsQuadraticGeneralizedEigenWithRealMat(
                    minBeta,
                    maxBeta,
                    k0,
                    out eVals, out eVecs);
            }

            SortEVals(eVals, eVecs);
            AdjustPhaseEVecs(eVecs);

            // 伝搬モードに限定
            System.Numerics.Complex[] propBetas;
            System.Numerics.Complex[][] propEVecs;
            GetPropagationModes(eVals, eVecs, out propBetas, out propEVecs);

            System.Numerics.Complex[] defectBetas;
            System.Numerics.Complex[][] defectEVecs;
            GetDefectModes(
                k0, minBeta, maxBeta,
                propBetas, propEVecs,
                out defectBetas, out defectEVecs);
            // モード追跡
            TraceMode(defectBetas, defectEVecs);

            System.Numerics.Complex[][] defectFxEVecs0;
            System.Numerics.Complex[][] defectFyEVecs0;
            CalcFuValues(defectBetas, defectEVecs, out defectFxEVecs0, out defectFyEVecs0);
            System.Numerics.Complex[][] defectFxEVecs = null;
            if (isYDirectionPeriodic)
            {
                defectFxEVecs = defectFyEVecs0;
            }
            else
            {
                defectFxEVecs = defectFxEVecs0;
            }

            // 境界1上のみのベクトル
            int bcNodeCnt = WgPortInfo.BcNodess[0].Count;
            int modeCnt = defectBetas.Length;
            System.Numerics.Complex[][] defectBcEVecs = new System.Numerics.Complex[modeCnt][];
            System.Numerics.Complex[][] defectBcFxEVecs = new System.Numerics.Complex[modeCnt][];
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                System.Numerics.Complex[] defectEVec = defectEVecs[iMode];
                System.Numerics.Complex[] defectFxEVec = defectFxEVecs[iMode];
                System.Numerics.Complex[] defectBcEVec = new System.Numerics.Complex[bcNodeCnt];
                System.Numerics.Complex[] defectBcFxEVec = new System.Numerics.Complex[bcNodeCnt];
                defectBcEVecs[iMode] = defectBcEVec;
                defectBcFxEVecs[iMode] = defectBcFxEVec;
                for (int i = 0; i < bcNodeCnt; i++)
                {
                    defectBcEVec[i] = defectEVec[i];
                    defectBcFxEVec[i] = defectFxEVec[i];
                }
            }

            // 規格化定数を計算する
            var RyyZ = new IvyFEM.Lapack.ComplexMatrix(RyyB1);
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                System.Numerics.Complex beta = defectBetas[iMode];
                System.Numerics.Complex[] defectBcEVec = defectBcEVecs[iMode];
                System.Numerics.Complex[] defectBcFxEVec = defectBcFxEVecs[iMode];
                System.Numerics.Complex betaPeriodic = ToBetaPeriodic(beta, periodicDistance);

                System.Numerics.Complex[] eVecModify = new System.Numerics.Complex[bcNodeCnt];
                for (int i = 0; i < bcNodeCnt; i++)
                {
                    eVecModify[i] = defectBcEVec[i] - defectBcFxEVec[i] / 
                        (System.Numerics.Complex.ImaginaryOne * betaPeriodic);
                }
                var work = RyyZ * defectBcEVec;
                System.Numerics.Complex work2 = 
                    IvyFEM.Lapack.Functions.zdotu(IvyFEM.Lapack.Utils.Conjugate(eVecModify), work);
                double replacedMu0 = IsTMMode ? Constants.Ep0 : Constants.Mu0;
                System.Numerics.Complex d = System.Numerics.Complex.Sqrt(
                    omega * replacedMu0 * 
                    (1.0 / beta.Magnitude) * 
                    (System.Numerics.Complex.Conjugate(beta) /
                    System.Numerics.Complex.Conjugate(betaPeriodic)) *
                    (1.0 / work2));
                defectBcEVec = IvyFEM.Lapack.Functions.zscal(defectBcEVec, d);
                defectBcFxEVec = IvyFEM.Lapack.Functions.zscal(defectBcFxEVec, d);
                defectBcEVecs[iMode] = defectBcEVec;
                defectBcFxEVecs[iMode] = defectBcFxEVec;
            }

            // 後片付け
            // 回転移動設定をクリアする
            World.RotAngle = 0;
            World.RotOrigin = null;

            Betas = defectBetas;
            EVecs = defectEVecs;
            FxEVecs = defectFxEVecs;
            BcEVecs = defectBcEVecs;
            BcFxEVecs = defectBcFxEVecs;
        }

        // 反復で一般化固有値問題として解く
        private void SolveIterAsGeneralizedEigen(
            double k0, double minBeta, double maxBeta,
            out System.Numerics.Complex[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            System.Diagnostics.Debug.Assert(WgPortInfo.IsSolveIter);
            System.Diagnostics.Debug.Assert(WgPortInfo.ModeCountToSolveIter > 0);
            // 前の固有ベクトルを退避する
            System.Numerics.Complex[][] savePrevModeEVecs = null;
            if (WgPortInfo.PrevModeEVecs != null && WgPortInfo.PrevModeEVecs.Length > 0)
            {
                System.Numerics.Complex[][] prevModeEVecs = WgPortInfo.PrevModeEVecs;
                int prevModeCnt = prevModeEVecs.Length;
                savePrevModeEVecs = new System.Numerics.Complex[prevModeCnt][];
                for (int iMode = 0; iMode < prevModeCnt; iMode++)
                {
                    System.Numerics.Complex[] prevModeEVec = prevModeEVecs[iMode];
                    int nodeCnt = prevModeEVec.Length;
                    System.Numerics.Complex[] savePrevEVec = new System.Numerics.Complex[nodeCnt];
                    prevModeEVec.CopyTo(savePrevEVec, 0);
                }
            }

            int modeCntToSolve = WgPortInfo.ModeCountToSolveIter;
            IList<System.Numerics.Complex> eVals1 = new List<System.Numerics.Complex>();
            IList<System.Numerics.Complex[]> eVecs1 = new List<System.Numerics.Complex[]>();
            for (int iMode = 0; iMode < modeCntToSolve; iMode++)
            {
                System.Numerics.Complex[] eVals0;
                System.Numerics.Complex[][] eVecs0;
                SolveIterEachModeAsGeneralizedEigen(iMode, k0, minBeta, maxBeta, out eVals0, out eVecs0);
                if (eVals0.Length > 0)
                {
                    System.Diagnostics.Debug.Assert(eVals0.Length == 1);
                    eVals1.Add(eVals0[0]);
                    eVecs1.Add(eVecs0[0]);
                }
            }
            eVals = eVals1.ToArray();
            eVecs = eVecs1.ToArray();

            // 前の固有ベクトルを復活する
            WgPortInfo.PrevModeEVecs = savePrevModeEVecs;
        }

        // 反復で一般化固有値問題として解く(モード単位)
        private void SolveIterEachModeAsGeneralizedEigen(
            int targetModeIndex, double k0, double minBeta, double maxBeta,
            out System.Numerics.Complex[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            eVals = null;
            eVecs = null;

            // 境界1 + 内部領域
            int nodeCnt1 = KMat.RowLength;

            System.Numerics.Complex resultBeta = 0;
            System.Numerics.Complex[] resultEVec = null;

            int itrMax = 400;
            //double resMin = 1.0e-12;
            double resMin = 1.0e-6;
            double prevRes = double.MaxValue;
            int itr = 0;
            for (itr = 0; itr < itrMax; itr++)
            {
                System.Numerics.Complex[] workEVals = null;
                System.Numerics.Complex[,] workEVecs = null;
                //  { [K] - jβ[C] - β^2[M] }{Φ}= {0}
                var A = new IvyFEM.Lapack.ComplexMatrix(nodeCnt1, nodeCnt1);
                var B = new IvyFEM.Lapack.ComplexMatrix(nodeCnt1, nodeCnt1);
                for (int i = 0; i < nodeCnt1; i++)
                {
                    for (int j = 0; j < nodeCnt1; j++)
                    {
                        A[i, j] = KMat[i, j] -
                            System.Numerics.Complex.ImaginaryOne * resultBeta * CMat[i, j];
                        B[i, j] = MMat[i, j];
                    }
                }
                // 伝搬定数が実数の時のみに限定
                if (Math.Abs(resultBeta.Imaginary) >= Constants.PrecisionLowerLimit)
                {
                    System.Diagnostics.Debug.WriteLine(
                        "!!!!!!!!Not propagation mode. Skip calculate: betaToSolve: {0} + {1}i",
                        resultBeta.Real, resultBeta.Imaginary);
                    break;
                }

                int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
                int bcNodeCnt = WgPortInfo.BcNodess[0].Count;

                double[] eVals0;
                System.Numerics.Complex[][] eVecs0;
                System.Diagnostics.Debug.WriteLine("solve eigen");
                {
                    // 固有値問題を解く
                    System.Numerics.Complex[][] tmpEVecs0;
                    SolveComplexHermiteBandGeneralizedEigen(A, B, out eVals0, out tmpEVecs0);

                    // eVecs0は節点番号順に並んでいないので並び替える
                    bool isPortBc2Reverse = WgPortInfo.IsPortBc2Reverse;
                    int modeCnt = tmpEVecs0.Length;
                    eVecs0 = new System.Numerics.Complex[modeCnt][];
                    for (int iMode = 0; iMode < modeCnt; iMode++)
                    {
                        System.Numerics.Complex[] tmpEVec0 = tmpEVecs0[iMode];
                        System.Numerics.Complex[] eVec0 = new System.Numerics.Complex[nodeCnt];
                        eVecs0[iMode] = eVec0;
                        for (int i = 0; i < bcNodeCnt; i++)
                        {
                            eVec0[i] = tmpEVec0[i];
                        }
                        for (int i = 0; i < bcNodeCnt; i++)
                        {
                            int i1 = isPortBc2Reverse ?
                                (bcNodeCnt * 2 - 1 - i) :
                                (i + bcNodeCnt);
                            eVec0[i1] = tmpEVec0[i];
                        }
                        for (int i = bcNodeCnt; i < nodeCnt1; i++)
                        {
                            int i1 = i + bcNodeCnt;
                            eVec0[i1] = tmpEVec0[i];
                        }
                    }
                }

                System.Numerics.Complex[] betas1;
                {
                    int modeCnt = eVals0.Length;
                    betas1 = new System.Numerics.Complex[modeCnt];
                    for (int iMode = 0; iMode < modeCnt; iMode++)
                    {
                        // λ = β^2と置いた場合( β = sqrt(λ) )
                        betas1[iMode] = System.Numerics.Complex.Sqrt(eVals0[iMode]);
                    }
                }
                System.Numerics.Complex[][] eVecs1 = new System.Numerics.Complex[eVecs0.Length][];
                for (int i = 0; i < eVecs0.Length; i++)
                {
                    System.Numerics.Complex[] eVec0 = eVecs0[i];
                    System.Numerics.Complex[] eVec1 = new System.Numerics.Complex[eVec0.Length];
                    eVecs1[i] = eVec1;
                    eVec0.CopyTo(eVec1, 0);
                }
                SortEVals(betas1, eVecs1);
                System.Numerics.Complex[] defectBetas;
                System.Numerics.Complex[][] defectEVecs;
                GetDefectModes(
                    k0, minBeta, maxBeta,
                    betas1, eVecs1,
                    out defectBetas, out defectEVecs);

                System.Numerics.Complex prevBeta = resultBeta;
                if (targetModeIndex < defectBetas.Length)
                {
                    resultBeta = defectBetas[targetModeIndex];
                    resultEVec = defectEVecs[targetModeIndex];
                }
                else
                {
                    // 欠陥モードがない
                    System.Diagnostics.Debug.WriteLine("!!!!!!!! no defect mode for targetModeIndex=" + targetModeIndex);
                    resultBeta = 0.0;
                    resultEVec = null;
                    break;
                }

                // 収束判定
                double res = (prevBeta - resultBeta).Magnitude;
                if (res < resMin)
                {
                    System.Diagnostics.Debug.WriteLine(
                        "converged itr:{0} beta:{1} + {2} i", itr, resultBeta.Real, resultBeta.Imaginary);
                    break;
                }
                // 発散判定
                if (itr >= 20 && Math.Abs(res) > Math.Abs(prevRes))
                {
                    System.Diagnostics.Debug.WriteLine(
                        "!!!!!!!! Not converged : prevRes = {0} res = {1}", prevRes, res);
                    resultBeta = 0.0;
                    resultEVec = null;
                    break;
                }
                prevRes = res;
                // check
                if (itr % 20 == 0 && itr >= 20)
                {
                    System.Diagnostics.Debug.WriteLine("itr: {0}", itr);
                }
            }
            if (itr == itrMax)
            {
                System.Diagnostics.Debug.WriteLine(
                    "!!!!!!!! Not converged itr:{0} beta:{1} + {2} i", itr, resultBeta.Real, resultBeta.Imaginary);
            }
            else if (itr >= 20 && (Math.Abs(resultBeta.Real) >= Constants.PrecisionLowerLimit ||
                Math.Abs(resultBeta.Imaginary) >= Constants.PrecisionLowerLimit))
            {
                System.Diagnostics.Debug.WriteLine("converged but too slow!!!: itr: {0}", itr);
            }
            else
            {

            }

            if (resultEVec == null)
            {
                eVals = new System.Numerics.Complex[0];
                eVecs = new System.Numerics.Complex[0][];
            }
            else
            {
                eVals = new System.Numerics.Complex[1] { resultBeta };
                eVecs = new System.Numerics.Complex[1][] { resultEVec };
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

            // 行列要素check
            {
                double th = 1.0e-12;
                for (int i = 0; i < matLen; i++)
                {
                    // [B]の正定値行列チェック
                    System.Diagnostics.Debug.Assert(sparseB2[i, i].Real > 0);
                    for (int j = i; j < matLen; j++)
                    {
                        // [A]はエルミート行列
                        System.Diagnostics.Debug.Assert(
                            (sparseA2[i, j] - System.Numerics.Complex.Conjugate(sparseA2[j, i])).Magnitude < th);
                        // [B]はエルミート行列
                        System.Diagnostics.Debug.Assert(
                            (sparseB2[i, j] - System.Numerics.Complex.Conjugate(sparseB2[j, i])).Magnitude < th);
                    }
                }
            }
            var hermiteBandA = new IvyFEM.Lapack.ComplexHermitianBandMatrix(sparseA2);
            var hermiteBandB = new IvyFEM.Lapack.ComplexHermitianBandMatrix(sparseB2);

            System.Numerics.Complex[][] eVecs2;
            /*
            int ret = IvyFEM.Lapack.Functions.zhbgv(
                hermiteBandA.Buffer,
                hermiteBandA.RowLength, hermiteBandA.ColumnLength, hermiteBandA.SuperdiaLength,
                hermiteBandB.Buffer,
                hermiteBandB.RowLength, hermiteBandB.ColumnLength, hermiteBandB.SuperdiaLength,
                out eVals, out eVecs2);
            */
            int ret = IvyFEM.Lapack.Functions.zhbgv_dirty(
                hermiteBandA.Buffer,
                hermiteBandA.RowLength, hermiteBandA.ColumnLength, hermiteBandA.SuperdiaLength,
                hermiteBandB.Buffer,
                hermiteBandB.RowLength, hermiteBandB.ColumnLength, hermiteBandB.SuperdiaLength,
                out eVals, out eVecs2);

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

        private void SolveAsGeneralizedEigenWithRealMat(
            out System.Numerics.Complex[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            eVals = null;
            eVecs = null;

            // 境界1 + 内部領域
            int nodeCnt1 = KMat.RowLength;

            // 一般化固有値解析(実行列として解く)
            var A = new IvyFEM.Lapack.DoubleMatrix(nodeCnt1 * 2, nodeCnt1 * 2);
            var B = new IvyFEM.Lapack.DoubleMatrix(nodeCnt1 * 2, nodeCnt1 * 2);
            for (int i = 0; i < nodeCnt1; i++)
            {
                for (int j = 0; j < nodeCnt1; j++)
                {
                    A[i, j] = 0.0;
                    A[i, j + nodeCnt1] = (i == j) ? 1.0 : 0.0;
                    //  { [K] - jβ[C] - β^2[M] }{Φ}= {0}
                    // λ = -jβと置いた場合
                    A[i + nodeCnt1, j] = -1.0 * KMat[i, j];
                    A[i + nodeCnt1, j + nodeCnt1] = -1.0 * CMat[i, j];

                    B[i, j] = (i == j) ? 1.0 : 0.0;
                    B[i, j + nodeCnt1] = 0.0;
                    B[i + nodeCnt1, j] = 0.0;
                    // λ = -jβと置いた場合
                    B[i + nodeCnt1, j + nodeCnt1] = MMat[i, j];
                }
            }

            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            int bcNodeCnt = WgPortInfo.BcNodess[0].Count;

            System.Diagnostics.Debug.WriteLine("solve eigen");
            {
                // 固有値問題を解く
                System.Numerics.Complex[][] eVecs0 = null;
                int ret = -1;
                try
                {
                    /*
                    ret = IvyFEM.Lapack.Functions.dggev(A.Buffer, A.RowLength, A.ColumnLength,
                        B.Buffer, B.RowLength, B.ColumnLength,
                        out eVals, out eVecs0);
                    */
                    ret = IvyFEM.Lapack.Functions.dggev_dirty(A.Buffer, A.RowLength, A.ColumnLength,
                        B.Buffer, B.RowLength, B.ColumnLength,
                        out eVals, out eVecs0);
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
                    eVecs0 = new System.Numerics.Complex[n][];
                    for (int iMode = 0; iMode < n; iMode++)
                    {
                        eVecs0[iMode] = new System.Numerics.Complex[n];
                    }
                }

                // eVecs0は節点番号順に並んでいないので並び替える
                bool isPortBc2Reverse = WgPortInfo.IsPortBc2Reverse;
                int modeCnt = eVecs0.Length;
                eVecs = new System.Numerics.Complex[modeCnt][];
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    System.Numerics.Complex[] eVec0 = eVecs0[iMode];
                    System.Numerics.Complex[] eVec = new System.Numerics.Complex[nodeCnt];
                    eVecs[iMode] = eVec;
                    for (int i = 0; i < bcNodeCnt; i++)
                    {
                        eVec[i] = eVec0[i];
                    }
                    for (int i = 0; i < bcNodeCnt; i++)
                    {
                        int i1 = isPortBc2Reverse ?
                            (bcNodeCnt * 2 - 1 - i) :
                            (i + bcNodeCnt);
                        eVec[i1] = eVec0[i];
                    }
                    for (int i = bcNodeCnt; i < nodeCnt1; i++)
                    {
                        int i1 = i + bcNodeCnt;
                        eVec[i1] = eVec0[i];
                    }
                }
            }

            // βに変換
            {
                int modeCnt = eVals.Length;
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    System.Numerics.Complex eVal = eVals[iMode];
                    // λ = -jβと置いた場合(β = jλ)
                    System.Numerics.Complex beta = eVal * System.Numerics.Complex.ImaginaryOne;
                    if (beta.Imaginary > 0)
                    {
                        beta = System.Numerics.Complex.Conjugate(beta);
                    }
                    eVals[iMode] = beta;
                }
            }
        }

        // 標準固有値問題として解く
        private void SolveAsStandardEigenWithRealMat(
            out System.Numerics.Complex[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            eVals = null;
            eVecs = null;

            // 境界1 + 内部領域
            int nodeCnt1 = KMat.RowLength;
            // 非線形固有値問題
            //  { [K] - jβ[C] - β^2[M] }{Φ}= {0}
            //  λ= - jβとおくと
            //  [K] + λ[C] + λ^2[M]{Φ}= {0}

            // [M]の逆行列が存在する緩慢変化包絡線近似の場合のみ有効な方法
            //   Φを直接解く場合は使えない
            System.Diagnostics.Debug.WriteLine("calc [M]-1");
            // [M]の逆行列を求める
            var invMMat = IvyFEM.Lapack.DoubleMatrix.Inverse(MMat);
            System.Diagnostics.Debug.WriteLine("calc [M]-1[K]");
            // [M]-1[K]
            var invMKMat = invMMat * KMat;
            System.Diagnostics.Debug.WriteLine("calc [M]-1[C]");
            // [M]-1[C]
            var invMCMat = invMMat * CMat;

            var A = new IvyFEM.Lapack.DoubleMatrix(nodeCnt1 * 2, nodeCnt1 * 2);
            System.Diagnostics.Debug.WriteLine("set [A]");
            for (int i = 0; i < nodeCnt1; i++)
            {
                for (int j = 0; j < nodeCnt1; j++)
                {
                    A[i, j] = 0.0;
                    A[i, j + nodeCnt1] = i == j ? 1.0 : 0.0;
                    //  { [K] - jβ[C] - β^2[M] }{Φ}= {0}
                    // λ = -jβと置いた場合
                    A[i + nodeCnt1, j] = -1.0 * invMKMat[i, j];
                    A[i + nodeCnt1, j + nodeCnt1] = -1.0 * invMCMat[i, j];
                }
            }

            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            int bcNodeCnt = WgPortInfo.BcNodess[0].Count;

            System.Diagnostics.Debug.WriteLine("solve eigen");
            {
                // 固有値問題を解く
                System.Numerics.Complex[][] eVecs0 = null;
                int ret = -1;
                try
                {
                    /*
                    ret = IvyFEM.Lapack.Functions.dgeev(
                        A.Buffer, A.RowLength, A.ColumnLength,
                        out eVals, out eVecs0);
                    */
                    ret = IvyFEM.Lapack.Functions.dgeev_dirty(
                        A.Buffer, A.RowLength, A.ColumnLength,
                        out eVals, out eVecs0);
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
                    eVecs0 = new System.Numerics.Complex[n][];
                    for (int iMode = 0; iMode < n; iMode++)
                    {
                        eVecs0[iMode] = new System.Numerics.Complex[n];
                    }
                }

                // 縮小前の固有ベクトルを求める
                bool isPortBc2Reverse = WgPortInfo.IsPortBc2Reverse;
                int modeCnt = eVecs0.Length;
                eVecs = new System.Numerics.Complex[modeCnt][];
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    System.Numerics.Complex[] eVec0 = eVecs0[iMode];
                    System.Numerics.Complex[] eVec = new System.Numerics.Complex[nodeCnt];
                    eVecs[iMode] = eVec;
                    for (int i = 0; i < bcNodeCnt; i++)
                    {
                        eVec[i] = eVec0[i];
                    }
                    for (int i = 0; i < bcNodeCnt; i++)
                    {
                        int i1 = isPortBc2Reverse ?
                            (bcNodeCnt * 2 - 1 - i) :
                            (i + bcNodeCnt);
                        eVec[i1] = eVec0[i];
                    }
                    for (int i = bcNodeCnt; i < nodeCnt1; i++)
                    {
                        int i1 = i + bcNodeCnt;
                        eVec[i1] = eVec0[i];
                    }
                }
            }

            // βに変換
            {
                int modeCnt = eVals.Length;
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    System.Numerics.Complex eVal = eVals[iMode];
                    // λ = -jβと置いた場合(β = jλ)
                    System.Numerics.Complex beta = eVal * System.Numerics.Complex.ImaginaryOne;
                    if (beta.Imaginary > 0)
                    {
                        beta = System.Numerics.Complex.Conjugate(beta);
                    }
                    eVals[iMode] = beta;
                }
            }
        }

        // Φを直接解く
        private void SolveNonSVEAModeAsQuadraticGeneralizedEigenWithRealMat(
            double minBeta,
            double maxBeta,
            double k0,
            out System.Numerics.Complex[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            eVals = null;
            eVecs = null;

            int nodeCnt = KMat0.RowLength;
            int bcNodeCnt = WgPortInfo.BcNodess[0].Count;
            System.Diagnostics.Debug.Assert(WgPortInfo.BcNodess[0].Count == WgPortInfo.BcNodess[1].Count);
            int nodeCnt1 = nodeCnt - bcNodeCnt; // 境界2を除く
            bool isPortBc2Reverse = WgPortInfo.IsPortBc2Reverse;
            double periodicDistance = WgPortInfo.PeriodicDistanceX;

            // 複素モード、エバネセントモードの固有ベクトル計算をスキップする？ (計算時間短縮するため)
            bool isSkipCalcComplexAndEvanescentModeVec = true;
            // 境界1のみの式に変換
            int innerNodeCnt = nodeCnt - bcNodeCnt * 2;
            var P11 = new IvyFEM.Lapack.DoubleMatrix(bcNodeCnt, bcNodeCnt);
            var P10 = new IvyFEM.Lapack.DoubleMatrix(bcNodeCnt, innerNodeCnt);
            var P12 = new IvyFEM.Lapack.DoubleMatrix(bcNodeCnt, bcNodeCnt);
            var P01 = new IvyFEM.Lapack.DoubleMatrix(innerNodeCnt, bcNodeCnt);
            var P00 = new IvyFEM.Lapack.DoubleMatrix(innerNodeCnt, innerNodeCnt);
            var P02 = new IvyFEM.Lapack.DoubleMatrix(innerNodeCnt, bcNodeCnt);
            var P21 = new IvyFEM.Lapack.DoubleMatrix(bcNodeCnt, bcNodeCnt);
            var P20 = new IvyFEM.Lapack.DoubleMatrix(bcNodeCnt, innerNodeCnt);
            var P22 = new IvyFEM.Lapack.DoubleMatrix(bcNodeCnt, bcNodeCnt);

            for (int i = 0; i < bcNodeCnt; i++)
            {
                int iB2 = isPortBc2Reverse ? (bcNodeCnt * 2 - 1 - i) : (bcNodeCnt + i);
                for (int j = 0; j < bcNodeCnt; j++)
                {
                    int jB2 = isPortBc2Reverse ? (bcNodeCnt * 2 - 1 - j) : (bcNodeCnt + j);
                    // [K11]
                    P11[i, j] = KMat0[i, j];
                    // [K12]
                    P12[i, j] = KMat0[i, jB2];
                    // [K21]
                    P21[i, j] = KMat0[iB2, j];
                    // [K22]
                    P22[i, j] = KMat0[iB2, jB2];
                }
                for (int j = 0; j < innerNodeCnt; j++)
                {
                    // [K10]
                    P10[i, j] = KMat0[i, j + bcNodeCnt * 2];
                    // [K20]
                    P20[i, j] = KMat0[iB2, j + bcNodeCnt * 2];
                }
            }
            for (int i = 0; i < innerNodeCnt; i++)
            {
                for (int j = 0; j < bcNodeCnt; j++)
                {
                    int jB2 = isPortBc2Reverse ? (bcNodeCnt * 2 - 1 - j) : (bcNodeCnt + j);
                    // [K01]
                    P01[i, j] = KMat0[i + bcNodeCnt * 2, j];
                    // [K02]
                    P02[i, j] = KMat0[i + bcNodeCnt * 2, jB2];
                }
                for (int j = 0; j < innerNodeCnt; j++)
                {
                    // [K00]
                    P00[i, j] = KMat0[i + bcNodeCnt * 2, j + bcNodeCnt * 2];
                }
            }

            System.Diagnostics.Debug.WriteLine("setup [K]B [C]B [M]B");
            var invP00 = IvyFEM.Lapack.DoubleMatrix.Inverse(P00);
            var P10xinvP00 = P10 * invP00;
            var P20xinvP00 = P20 * invP00;
            // for [C]B
            var P10xinvP00xP01 = P10xinvP00 * P01;
            var P20xinvP00xP02 = P20xinvP00 * P02;
            // for [M]B
            var P10xinvP00xP02 = P10xinvP00 * P02;
            // for [K]B
            var P20xinvP00xP01 = P20xinvP00 * P01;
            // [C]B
            var CMatB = new IvyFEM.Lapack.DoubleMatrix(bcNodeCnt, bcNodeCnt);
            // [M]B
            var MMatB = new IvyFEM.Lapack.DoubleMatrix(bcNodeCnt, bcNodeCnt);
            // [K]B
            var KMatB = new IvyFEM.Lapack.DoubleMatrix(bcNodeCnt, bcNodeCnt);
            for (int i = 0; i < bcNodeCnt; i++)
            {
                for (int j = 0; j < bcNodeCnt; j++)
                {
                    CMatB[i, j] =
                        -P10xinvP00xP01[i, j]
                        + P11[i, j]
                        - P20xinvP00xP02[i, j]
                        + P22[i, j];
                    MMatB[i, j] =
                        -P10xinvP00xP02[i, j]
                        + P12[i, j];
                    KMatB[i, j] =
                        -P20xinvP00xP01[i, j]
                        + P21[i, j];
                }
            }

            // 非線形固有値問題
            //  [K] + λ[C] + λ^2[M]{Φ}= {0}
            // 一般化固有値解析(実行列として解く)
            var A = new IvyFEM.Lapack.DoubleMatrix(bcNodeCnt * 2, bcNodeCnt * 2);
            var B = new IvyFEM.Lapack.DoubleMatrix(bcNodeCnt * 2, bcNodeCnt * 2);
            for (int i = 0; i < bcNodeCnt; i++)
            {
                for (int j = 0; j < bcNodeCnt; j++)
                {
                    A[i, j] = 0.0;
                    A[i, j + bcNodeCnt] = (i == j) ? 1.0 : 0.0;
                    A[i + bcNodeCnt, j] = -1.0 * KMatB[i, j];
                    A[i + bcNodeCnt, j + bcNodeCnt] = -1.0 * CMatB[i, j];
                }
            }
            for (int i = 0; i < bcNodeCnt; i++)
            {
                for (int j = 0; j < bcNodeCnt; j++)
                {
                    B[i, j] = (i == j) ? 1.0 : 0.0;
                    B[i , j + bcNodeCnt] = 0.0;
                    B[i + bcNodeCnt, j] = 0.0;
                    B[i + bcNodeCnt, j + bcNodeCnt] = MMatB[i, j];
                }
            }

            System.Diagnostics.Debug.WriteLine("solve eigen");
            {
                // 固有値問題を解く
                System.Numerics.Complex[][] eVecs0 = null;
                int ret = -1;
                try
                {
                    /*
                    ret = IvyFEM.Lapack.Functions.dggev(A.Buffer, A.RowLength, A.ColumnLength,
                        B.Buffer, B.RowLength, B.ColumnLength,
                        out eVals, out eVecs0);
                    */
                    ret = IvyFEM.Lapack.Functions.dggev_dirty(A.Buffer, A.RowLength, A.ColumnLength,
                        B.Buffer, B.RowLength, B.ColumnLength,
                        out eVals, out eVecs0);
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
                    eVecs0 = new System.Numerics.Complex[n][];
                    for (int iMode = 0; iMode < n; iMode++)
                    {
                        eVecs0[iMode] = new System.Numerics.Complex[n];
                    }
                }

                int modeCnt = eVals.Length;
                // βに変換
                double th = 1.0e-12;
                for (int i = 0; i < eVals.Length; i++)
                {
                    System.Numerics.Complex eVal = eVals[i];
                    if ((Math.Abs(eVal.Real) < th &&
                        Math.Abs(eVal.Imaginary) < th)
                        || double.IsInfinity(eVal.Real) || double.IsInfinity(eVal.Imaginary)
                        || double.IsNaN(eVal.Real) || double.IsNaN(eVal.Imaginary)
                        )
                    {
                        // 無効な固有値
                        eVals[i] = System.Numerics.Complex.ImaginaryOne * double.MaxValue;
                    }
                    else
                    {
                        System.Numerics.Complex beta = 
                            -1.0 * System.Numerics.Complex.Log(eVal) / 
                            (System.Numerics.Complex.ImaginaryOne * periodicDistance);
                        eVals[i] = beta;
                    }
                }

                // 2次元配列に格納する ({Φ}のみ格納)
                eVecs = new System.Numerics.Complex[modeCnt][];
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    eVecs[iMode] = new System.Numerics.Complex[nodeCnt];
                }

                System.Diagnostics.Debug.WriteLine("calc {Φ}0");
                var invP00xP01 = invP00 * P01;
                var invP00xP02 = invP00 * P02;
                var transMat = new IvyFEM.Lapack.ComplexMatrix(innerNodeCnt, bcNodeCnt);
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    System.Numerics.Complex beta = eVals[iMode];
                    // 複素モード、エバネセントモードの固有モード計算をスキップする？
                    if (isSkipCalcComplexAndEvanescentModeVec)
                    {
                        if (Math.Abs(beta.Imaginary) >= th)
                        {
                            // 複素モード、エバネセントモードの固有モード計算をスキップする
                            continue;
                        }
                    }

                    System.Numerics.Complex expA = 
                        System.Numerics.Complex.Exp(
                            -1.0 * System.Numerics.Complex.ImaginaryOne * beta * periodicDistance);
                    System.Numerics.Complex[] eVec = eVecs[iMode];
                    System.Numerics.Complex[] eVec0 = eVecs0[iMode];
                    System.Diagnostics.Debug.Assert(eVec0.Length == bcNodeCnt * 2);
                    System.Numerics.Complex[] fVecB = new System.Numerics.Complex[bcNodeCnt];
                    ///////////////////////////////
                    // {Φ}Bを格納
                    for (int i = 0; i < bcNodeCnt; i++)
                    {
                        System.Numerics.Complex value = eVec0[i];
                        eVec[i] = value;
                        fVecB[i] = value;
                    }
                    for (int i = 0; i < bcNodeCnt; i++)
                    {
                        int i1 = isPortBc2Reverse ?
                            (bcNodeCnt * 2 - 1 - i) :
                            (i + bcNodeCnt);
                        System.Numerics.Complex value = eVec0[i];
                        eVec[i1] = value * expA;
                    }

                    ///////////////////////////////
                    // {Φ}0を計算
                    //   変換行列を計算
                    for (int i = 0; i < innerNodeCnt; i++)
                    {
                        for (int j = 0; j < bcNodeCnt; j++)
                        {
                            transMat[i, j] = -1.0 * (invP00xP01[i, j] + expA * invP00xP02[i, j]);
                        }
                    }
                    //   {Φ}0を計算
                    System.Numerics.Complex[] fVecInner = transMat * fVecB;
                    //   {Φ}0を格納
                    for (int i = 0; i < innerNodeCnt; i++)
                    {
                        eVec[i + bcNodeCnt * 2] = fVecInner[i];
                    }
                }
            }

            ////////////////////////////////////////////////////////////////////
            {
                for (int iMode = 0; iMode < eVals.Length; iMode++)
                {
                    // 伝搬定数
                    System.Numerics.Complex beta = eVals[iMode];

                    System.Numerics.Complex betad = beta * periodicDistance;
                    double th = 1.0e-12;
                    if (
                        // [-π, 0]の解を[π, 2π]に移動する
                        ((minBeta * k0 * periodicDistance / (2.0 * Math.PI)) >= (0.5 - th)
                             && (minBeta * k0 * periodicDistance / (2.0 * Math.PI)) < 1.0
                             && Math.Abs(beta.Real) >= th
                             && Math.Abs(beta.Imaginary) < th
                             && (betad.Real / (2.0 * Math.PI)) >= (-0.5 - th)
                             && (betad.Real / (2.0 * Math.PI)) < 0.0)
                        // [0, π]の解を[2π, 3π]に移動する
                        || ((minBeta * k0 * periodicDistance / (2.0 * Math.PI)) >= (1.0 - th)
                                && (minBeta * k0 * periodicDistance / (2.0 * Math.PI)) < 1.5
                                && Math.Abs(beta.Real) >= th
                                && Math.Abs(beta.Imaginary) < th
                                && (betad.Real / (2.0 * Math.PI)) >= (0.0 - th)
                                && (betad.Real / (2.0 * Math.PI)) < 0.5)
                        )
                    {
                        // [0, π]の解を2πだけ移動する
                        double deltaPhase = 2.0 * Math.PI;
                        betad += deltaPhase; // 実部に加算する
                        beta = betad / periodicDistance;
                        //check
                        System.Diagnostics.Debug.WriteLine("shift beta * d / (2π): {0} + {1} i to {2} + {3} i",
                            eVals[iMode].Real * periodicDistance / (2.0 * Math.PI),
                            eVals[iMode].Imaginary * periodicDistance / (2.0 * Math.PI),
                            betad.Real / (2.0 * Math.PI),
                            betad.Imaginary / (2.0 * Math.PI));
                        // 再設定
                        eVals[iMode] = beta;
                    }
                }
            }
        }

        private void TraceMode(System.Numerics.Complex[] defectBetas, System.Numerics.Complex[][] defectEVecs)
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
                    int modeCnt = defectBetas.Length;
                    for (int iMode = 0; iMode < modeCnt; iMode++)
                    {
                        if (traceModeIndexss.IndexOf(iMode) != -1)
                        {
                            // すでに対応付けている
                            continue;
                        }
                        System.Numerics.Complex beta = defectBetas[iMode];
                        System.Numerics.Complex[] eVec = defectEVecs[iMode];
                        double norm;
                        bool isHitSameMode = IsSameMode(prevEVec, beta, eVec, out norm);
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
                    int modeCnt = defectBetas.Length;
                    System.Numerics.Complex[] prevBetas = new System.Numerics.Complex[modeCnt];
                    defectBetas.CopyTo(prevBetas, 0);
                    System.Numerics.Complex[][] prevEVecs = new System.Numerics.Complex[modeCnt][];
                    defectEVecs.CopyTo(prevEVecs, 0);
                    int index = 0;
                    for (int i = 0; i < traceModeIndexss.Count; i++)
                    {
                        int traceModeIndex = traceModeIndexss[i];
                        defectBetas[index] = prevBetas[traceModeIndex];
                        defectEVecs[index] = prevEVecs[traceModeIndex];
                        index++;
                    }
                    for (int iMode = 0; iMode < modeCnt; iMode++)
                    {
                        if (!traceModeIndexss.Contains(iMode))
                        {
                            defectBetas[index] = prevBetas[iMode];
                            defectEVecs[index] = prevEVecs[iMode];
                            index++;
                        }
                    }
                    System.Diagnostics.Debug.Assert(index == modeCnt);
                }
                {
                    int modeCnt = defectBetas.Length;
                    WgPortInfo.PrevModeEVecs = new System.Numerics.Complex[modeCnt][];
                    defectEVecs.CopyTo(WgPortInfo.PrevModeEVecs, 0);
                }
            }
            else
            {
                int modeCnt = defectBetas.Length;
                WgPortInfo.PrevModeEVecs = new System.Numerics.Complex[modeCnt][];
                defectEVecs.CopyTo(WgPortInfo.PrevModeEVecs, 0);
            }
        }

        // 伝搬モード一覧の取得
        private void GetPropagationModes(
            System.Numerics.Complex[] betas, System.Numerics.Complex[][] eVecs,
            out System.Numerics.Complex[] propBetas, out System.Numerics.Complex[][] propEVecs)
        {
            IList<int> indexs = new List<int>();
            int modeCnt = betas.Length;
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var beta = betas[iMode];
                if (Math.Abs(beta.Imaginary) < 1.0e-12)
                {
                    indexs.Add(iMode);
                }
            }

            int propCnt = indexs.Count;
            propBetas = new System.Numerics.Complex[propCnt];
            propEVecs = new System.Numerics.Complex[propCnt][];
            for (int i = 0; i < propCnt; i++)
            {
                int iMode = indexs[i];
                propBetas[i] = betas[iMode];
                propEVecs[i] = eVecs[iMode];
            }
        }

        // 欠陥モード一覧の取得
        private void GetDefectModes(
            double k0, double minBeta, double maxBeta,
            System.Numerics.Complex[] betas, System.Numerics.Complex[][] eVecs,
            out System.Numerics.Complex[] defectBetas, out System.Numerics.Complex[][] defectEVecs)
        {
            defectBetas = null;
            defectEVecs = null;
            IList<IList<int>> channelCoIdss = WgPortInfo.PCChannelCoIds;

            IList<System.Numerics.Complex> retBetas = new List<System.Numerics.Complex>();
            IList<System.Numerics.Complex[]> retEVecs = new List<System.Numerics.Complex[]>();
            int modeCnt = eVecs.Length;
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                System.Numerics.Complex beta = betas[iMode];
                System.Numerics.Complex[] eVec = eVecs[iMode];
                bool isDefect = IsDefectMode(k0, channelCoIdss, minBeta, maxBeta, beta, eVec);
                if (isDefect)
                {
                    retBetas.Add(beta);
                    System.Numerics.Complex[] copyEVec = new System.Numerics.Complex[eVec.Length];
                    eVec.CopyTo(copyEVec, 0);
                    retEVecs.Add(copyEVec);
                }
            }
            defectBetas = retBetas.ToArray();
            defectEVecs = retEVecs.ToArray();
        }

        // 欠陥モード？
        private bool IsDefectMode(
            double k0,
            IList<IList<int>> channelCoIdss,
            double minBeta,
            double maxBeta,
            System.Numerics.Complex beta,
            System.Numerics.Complex[] eVec)
        {
            bool isHit = false;
            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);

            double th = 1.0e-12;
            if (double.IsNaN(beta.Real) || double.IsNaN(beta.Imaginary))
            {
                // 固有値がNaNのモードは除外する
                return isHit;
            }
            else if (Math.Abs(beta.Real) >= th &&
                Math.Abs(beta.Imaginary) >= th)
            {
                // 複素モードは除外する
                return isHit;
            }
            else if (Math.Abs(beta.Real) >= th &&
                Math.Abs(beta.Imaginary) < th)
            {
                // 伝搬モード
                // 後進波は除外する
                if (beta.Real < 0)
                {
                    return isHit;
                }
            }
            else if (Math.Abs(beta.Real) < th &&
                Math.Abs(beta.Imaginary) >= th)
            {
                // 減衰モード
                //  利得のある波は除外する
                if (beta.Imaginary > 0)
                {
                    return isHit;
                }
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
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
                if (Math.Abs(beta.Real / k0) < minBeta || Math.Abs(beta.Real / k0) > maxBeta)
                {
                    System.Diagnostics.Debug.WriteLine(
                        "PCWaveguideMode: beta is invalid skip: β/k0 = {0} at k0 = {1}", beta.Real / k0, k0);
                }
                else
                {
                    System.Diagnostics.Debug.WriteLine(
                        "hit defect mode β/k0 = (" + beta.Real / k0 + ", " + beta.Imaginary + ")");
                    isHit = true;
                }
            }
            return isHit;
        }

        // 同じモード？
        private bool IsSameMode(
            System.Numerics.Complex[] prevEVec,
            System.Numerics.Complex beta,
            System.Numerics.Complex[] eVec,
            out double retNorm)
        {
            bool isHit = false;
            retNorm = 0;
            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            double freq = Frequency;
            double omega = 2.0 * Math.PI * freq;
            double k0 = omega / Constants.C0;

            if (prevEVec == null)
            {
                return isHit;
            }

            double th = 1.0e-12;
            if (Math.Abs(beta.Real) >= th &&
                Math.Abs(beta.Imaginary) >= th)
            {
                // 複素モードは除外する
                return isHit;
            }
            else if (Math.Abs(beta.Real) >= th &&
                Math.Abs(beta.Imaginary) < th)
            {
                // 伝搬モード
                // 後進波は除外する
                if (beta.Real < 0)
                {
                    return isHit;
                }
            }
            else if (Math.Abs(beta.Real) < th &&
                Math.Abs(beta.Imaginary) >= th)
            {
                // 減衰モード
                //  利得のある波は除外する
                if (beta.Imaginary > 0)
                {
                    return isHit;
                }
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
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
                    "β/k0(eVal)= {0} norm (prev * current)= {1}", (beta / k0), norm12.Magnitude);
            }
            return isHit;
        }

        // 周期構造導波路の伝搬定数に変換する(βdが[-π, π]に収まるように変換)
        private System.Numerics.Complex ToBetaPeriodic(System.Numerics.Complex beta, double periodicDistance)
        {
            // βの再変換
            System.Numerics.Complex expA = 
                System.Numerics.Complex.Exp(-1.0 * System.Numerics.Complex.ImaginaryOne * beta * periodicDistance);
            System.Numerics.Complex betaPeriodic = 
                -1.0 * System.Numerics.Complex.Log(expA) / (System.Numerics.Complex.ImaginaryOne * periodicDistance);
            return betaPeriodic;
        }

        private void NonSVEAEVecToSVEAEVec(
            System.Numerics.Complex beta,
            System.Numerics.Complex[] fVec, System.Numerics.Complex[] fxVec, System.Numerics.Complex[] fyVec)
        {
            System.Diagnostics.Debug.Assert(!WgPortInfo.IsSVEA);
            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            System.Diagnostics.Debug.Assert(nodeCnt == fVec.Length);
            System.Diagnostics.Debug.Assert(nodeCnt == fxVec.Length);
            System.Diagnostics.Debug.Assert(nodeCnt == fyVec.Length);

            IList<uint> bcEIds = WgPortInfo.BcEdgeIds1;
            uint eId1 = bcEIds[0];
            uint eId2 = bcEIds[bcEIds.Count - 1];
            Edge2D e1 = World.Mesh.Cad.GetEdge(eId1);
            Edge2D e2 = World.Mesh.Cad.GetEdge(eId2);
            OpenTK.Vector2d firstPt = e1.GetVertexCoord(true);
            OpenTK.Vector2d lastPt = e2.GetVertexCoord(false);
            double[] firstCoord = { firstPt.X, firstPt.Y };
            double[] lasCoord = { lastPt.X, lastPt.Y };

            double periodicDistance = WgPortInfo.PeriodicDistanceX;
            bool isYDirectionPeriodic = WgPortInfo.IsYDirectionPeriodic;
            System.Numerics.Complex betaPeriodic = ToBetaPeriodic(beta, periodicDistance);
            for (int nodeId = 0; nodeId < nodeCnt; nodeId++)
            {
                int coId = World.PortNode2Coord(QuantityId, PortId, nodeId);
                double[] coord = World.GetCoord(QuantityId, coId);
                double ptX = 0;
                if (isYDirectionPeriodic)
                {
                    // Y方向
                    ptX = coord[1] - firstCoord[1];
                }
                else
                {
                    // X方向
                    ptX = coord[0] - firstCoord[0];
                }

                System.Numerics.Complex expX =
                    System.Numerics.Complex.Exp(-1.0 * System.Numerics.Complex.ImaginaryOne * beta * ptX);

                System.Numerics.Complex fVal = fVec[nodeId];
                System.Numerics.Complex fxVal = fxVec[nodeId];
                System.Numerics.Complex fyVal = fyVec[nodeId];
                // ΦのSVEA(φ)
                System.Numerics.Complex fValSVEA = fVal / expX;
                fVec[nodeId] = fValSVEA;

                // SVEAの微分( exp(-jβx)dφ/dx = dΦ/dx + jβφexp(-jβx))
                if (isYDirectionPeriodic)
                {
                    // Y方向周期構造
                    fxVec[nodeId] = fxVal / expX;
                    fyVec[nodeId] = fyVal / expX + System.Numerics.Complex.ImaginaryOne * betaPeriodic * fValSVEA;
                }
                else
                {
                    // X方向周期構造
                    fxVec[nodeId] = fxVal / expX + System.Numerics.Complex.ImaginaryOne * betaPeriodic * fValSVEA;
                    fyVec[nodeId] = fyVal / expX;
                }
            }
        }

        private void SortEVals(System.Numerics.Complex[] betas, System.Numerics.Complex[][] eVecs)
        {
            int modeCnt = betas.Length;
            var eValEVecs = new List<KeyValuePair<System.Numerics.Complex, System.Numerics.Complex[]>>();
            for (int i = 0; i < modeCnt; i++)
            {
                eValEVecs.Add(
                    new KeyValuePair<System.Numerics.Complex, System.Numerics.Complex[]>(betas[i], eVecs[i]));
            }
            eValEVecs.Sort((a, b) =>
            {
                System.Numerics.Complex betaA = a.Key;
                System.Numerics.Complex betaB = b.Key;
                int cmp = 0;
                double th = 1.0e-12;
                if (Math.Abs(betaA.Real) < th && Math.Abs(betaB.Real) < th)
                {
                    double cmpfi = Math.Abs(betaA.Imaginary) - Math.Abs(betaB.Imaginary);
                    // 昇順
                    if (Math.Abs(cmpfi) < 1.0e-15)
                    {
                        cmp = 0;
                    }
                    else if (cmpfi > 0)
                    {
                        cmp = 1;
                    }
                    else
                    {
                        cmp = -1;
                    }
                }
                else if (Math.Abs(betaA.Real) < th && Math.Abs(betaB.Real) >= th)
                {
                    if (betaB.Real >= 0)
                    {
                        cmp = 1;
                    }
                    else
                    {
                        cmp = -1;
                    }
                }
                else if (Math.Abs(betaA.Real) >= th && Math.Abs(betaB.Real) < th)
                {
                    if (betaA.Real >= 0)
                    {
                        cmp = -1;
                    }
                    else
                    {
                        cmp = 1;
                    }
                }
                else
                {
                    double cmpf = betaA.Real - betaB.Real;
                    if (Math.Abs(cmpf) < 1.0e-15)
                    {
                        cmp = 0;
                    }
                    else
                    {
                        if (cmpf > 0)
                        {
                            cmp = -1;
                        }
                        else
                        {
                            cmp = 1;
                        }
                    }
                }
                return cmp;
            });

            for (int i = 0; i < modeCnt; i++)
            {
                betas[i] = eValEVecs[i].Key;
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

        private void CalcFuValues(
            System.Numerics.Complex[] betas,
            System.Numerics.Complex[][] fVecs,
            out System.Numerics.Complex[][] fxVecs, out System.Numerics.Complex[][] fyVecs)
        {
            int modeCnt = fVecs.Length;

            fxVecs = new System.Numerics.Complex[modeCnt][];
            fyVecs = new System.Numerics.Complex[modeCnt][];

            bool isSVEA = WgPortInfo.IsSVEA;
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                System.Numerics.Complex beta = betas[iMode];
                System.Numerics.Complex[] fVec = fVecs[iMode];
                System.Numerics.Complex[] fxVec;
                System.Numerics.Complex[] fyVec;
                CalcModeFuValues(fVec, out fxVec, out fyVec);

                // 境界1と境界2の節点の微分値は同じという条件を弱形式で課している為、微分値は同じにならない。
                // 加えて、CalcModeFuValuesは内部節点からの寄与を片側のみしか計算していない。
                // →境界の両側からの寄与を考慮する為に境界1の微分値と境界2の微分値を平均してみる
                var bcNodes1 = WgPortInfo.BcNodess[0];
                var bcNodes2 = WgPortInfo.BcNodess[1];
                int bcNodeCnt = bcNodes1.Count;
                bool isPortBc2Reverse = WgPortInfo.IsPortBc2Reverse;
                double periodicDistance = WgPortInfo.PeriodicDistanceX;
                for (int i = 0; i < bcNodeCnt; i++)
                {
                    int portNodeId1 = bcNodes1[i];
                    int i2 = isPortBc2Reverse ? (bcNodeCnt - 1 - i) : i;
                    int portNodeId2 = bcNodes2[i2];
                    if (isSVEA)
                    {
                        // 緩慢変化包絡線近似の場合は、Φ2 = Φ1
                        var fx = (fxVec[portNodeId1] + fxVec[portNodeId2]) / 2.0;
                        var fy = (fyVec[portNodeId1] + fyVec[portNodeId2]) / 2.0;
                        fxVec[portNodeId1] = fx;
                        fxVec[portNodeId2] = fx;
                        fyVec[portNodeId1] = fy;
                        fyVec[portNodeId2] = fy;
                    }
                    else
                    {
                        // 緩慢変化包絡線近似でない場合は、Φ2 = expA * Φ1
                        System.Numerics.Complex expA = 
                            System.Numerics.Complex.Exp(
                                -1.0 * System.Numerics.Complex.ImaginaryOne * beta * periodicDistance);
                        var fx = (fxVec[portNodeId1] + fxVec[portNodeId2] / expA) / 2.0;
                        var fy = (fyVec[portNodeId1] + fyVec[portNodeId2] / expA) / 2.0;
                        fxVec[portNodeId1] = fx;
                        fxVec[portNodeId2] = fx * expA;
                        fyVec[portNodeId1] = fy;
                        fyVec[portNodeId2] = fy * expA;
                    }
                }

                if (!isSVEA)
                {
                    NonSVEAEVecToSVEAEVec(beta, fVec, fxVec, fyVec);
                }

                fxVecs[iMode] = fxVec;
                fyVecs[iMode] = fyVec;
            }
        }

        // dF/dx dF/dy
        private void CalcModeFuValues(
            System.Numerics.Complex[] fVec,
            out System.Numerics.Complex[] fxVec, out System.Numerics.Complex[] fyVec)
        {
            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            System.Diagnostics.Debug.Assert(fVec.Length == nodeCnt);

            fxVec = new System.Numerics.Complex[nodeCnt];
            fyVec = new System.Numerics.Complex[nodeCnt];

            IList<uint> feIds = World.GetPeriodicPortTriangleFEIds(QuantityId, PortId);

            Dictionary<int, double> nodeAreaSum = new Dictionary<int, double>();

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

                double A = triFE.GetArea();
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int iNodeId = nodes[iNode];
                    if (iNodeId == -1)
                    {
                        continue;
                    }
                    double[] L = triFE.GetNodeL(iNode);
                    double[][] Nu = triFE.CalcNu(L);
                    double[] Nx = Nu[0];
                    double[] Ny = Nu[1];

                    System.Numerics.Complex fx = 0;
                    System.Numerics.Complex fy = 0;
                    for (int kNode = 0; kNode < elemNodeCnt; kNode++)
                    {
                        int kNodeId = nodes[kNode];
                        if (kNodeId == -1)
                        {
                            continue;
                        }
                        System.Numerics.Complex f = fVec[kNodeId];
                        fx += f * Nx[kNode];
                        fy += f * Ny[kNode];
                    }
                    fxVec[iNodeId] += fx * A;
                    fyVec[iNodeId] += fy * A;

                    if (nodeAreaSum.ContainsKey(iNodeId))
                    {
                        nodeAreaSum[iNodeId] += A;
                    }
                    else
                    {
                        nodeAreaSum[iNodeId] = A;
                    }
                }
            }

            System.Diagnostics.Debug.Assert(nodeAreaSum.Count == nodeCnt);
            for (int i = 0; i < nodeCnt; i++)
            {
                double areaSum = nodeAreaSum[i];
                fxVec[i] /= areaSum;
                fyVec[i] /= areaSum;
            }
        }

        public IvyFEM.Lapack.ComplexMatrix CalcBoundaryMatrix(
            double omega, System.Numerics.Complex[] betas,
            System.Numerics.Complex[][] ezEVecs, System.Numerics.Complex[][] ezXEVecs)
        {
            int bcNodeCnt = WgPortInfo.BcNodess[0].Count;
            double periodicDistance = WgPortInfo.PeriodicDistanceX;

            IvyFEM.Lapack.ComplexMatrix X = new Lapack.ComplexMatrix(bcNodeCnt, bcNodeCnt);

            int modeCnt = betas.Length;
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var beta = betas[iMode];
                var ezEVec = ezEVecs[iMode];
                var ezXEVec = ezXEVecs[iMode];

                var RyyZ = new IvyFEM.Lapack.ComplexMatrix(RyyB1);
                System.Diagnostics.Debug.Assert(ezEVec.Length == bcNodeCnt);
                System.Diagnostics.Debug.Assert(ezXEVec.Length == bcNodeCnt);
                System.Diagnostics.Debug.Assert(RyyZ.RowLength == bcNodeCnt);
                System.Diagnostics.Debug.Assert(RyyZ.ColumnLength == bcNodeCnt);

                var betaPeriodic = ToBetaPeriodic(beta, periodicDistance);
                var ezEVecModify = new System.Numerics.Complex[bcNodeCnt];
                for (int i = 0; i < bcNodeCnt; i++)
                {
                    ezEVecModify[i] = 
                        ezEVec[i] - ezXEVec[i] / (System.Numerics.Complex.ImaginaryOne * betaPeriodic);
                }
                var vec1 = RyyZ * ezEVecModify;
                var vec2 = RyyZ * IvyFEM.Lapack.Utils.Conjugate(ezEVecModify);

                double replacedMu0 = IsTMMode ? Constants.Ep0 : Constants.Mu0;
                for (int row = 0; row < bcNodeCnt; row++)
                {
                    for (int col = 0; col < bcNodeCnt; col++)
                    {
                        System.Numerics.Complex value = 
                            (System.Numerics.Complex.ImaginaryOne / (omega * replacedMu0)) *
                            betaPeriodic * beta.Magnitude *
                            (System.Numerics.Complex.Conjugate(betaPeriodic) /
                            System.Numerics.Complex.Conjugate(beta)) *
                            vec1[row] * vec2[col];

                        X[row, col] += value;
                    }
                }
            }
            return X;
        }

        public System.Numerics.Complex[] CalcIncidentVec(
            System.Numerics.Complex beta0,
            System.Numerics.Complex[] ezEVec0, System.Numerics.Complex[] ezXEVec0)
        {
            double periodicDistance = WgPortInfo.PeriodicDistanceX;

            System.Numerics.Complex[] I = null;

            int bcNodeCnt = WgPortInfo.BcNodess[0].Count;
            System.Diagnostics.Debug.Assert(ezEVec0.Length == bcNodeCnt);

            var RyyZ = new IvyFEM.Lapack.ComplexMatrix(RyyB1);

            var betaPeriodic0 = ToBetaPeriodic(beta0, periodicDistance);
            var ezEVecModify0 = new System.Numerics.Complex[bcNodeCnt];
            for (int i = 0; i < bcNodeCnt; i++)
            {
                ezEVecModify0[i] =
                    ezEVec0[i] - ezXEVec0[i] / (System.Numerics.Complex.ImaginaryOne * betaPeriodic0);
            }
            var vec1 = RyyZ * ezEVecModify0;
            var a1 = System.Numerics.Complex.ImaginaryOne * 2.0 * betaPeriodic0;
            vec1 = IvyFEM.Lapack.Functions.zscal(vec1, a1);
            I = vec1;
            return I;
        }

        public System.Numerics.Complex[] CalcSMatrix(double omega, int incidentModeId,
            System.Numerics.Complex[] betas,
            System.Numerics.Complex[][] ezEVecs, System.Numerics.Complex[][] ezXEVecs,
            System.Numerics.Complex[] Ez)
        {
            double periodicDistance = WgPortInfo.PeriodicDistanceX;

            int bcNodeCnt = WgPortInfo.BcNodess[0].Count;
            int modeCnt = betas.Length;
            System.Numerics.Complex[] S = new System.Numerics.Complex[modeCnt];

            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var beta = betas[iMode];
                var ezEVec = ezEVecs[iMode];
                var ezXEVec = ezXEVecs[iMode];
                var RyyZ = new IvyFEM.Lapack.ComplexMatrix(RyyB1);
                System.Diagnostics.Debug.Assert(ezEVec.Length == bcNodeCnt);
                System.Diagnostics.Debug.Assert(ezXEVec.Length == bcNodeCnt);
                System.Diagnostics.Debug.Assert(RyyZ.RowLength == bcNodeCnt);
                System.Diagnostics.Debug.Assert(RyyZ.ColumnLength == bcNodeCnt);

                var betaPeriodic = ToBetaPeriodic(beta, periodicDistance);
                var ezEVecModify = new System.Numerics.Complex[bcNodeCnt];
                for (int i = 0; i < bcNodeCnt; i++)
                {
                    ezEVecModify[i] =
                        ezEVec[i] - ezXEVec[i] / (System.Numerics.Complex.ImaginaryOne * betaPeriodic);
                }
                var vec1 = RyyZ * IvyFEM.Lapack.Utils.Conjugate(ezEVecModify);
                double replacedMu0 = IsTMMode ? Constants.Ep0 : Constants.Mu0;
                System.Numerics.Complex work1 = IvyFEM.Lapack.Functions.zdotu(vec1, Ez);
                var b = beta.Magnitude *
                    (System.Numerics.Complex.Conjugate(betaPeriodic) / 
                    System.Numerics.Complex.Conjugate(beta)) *
                    (1.0 / (omega * replacedMu0)) * work1;
                if (incidentModeId != -1 && incidentModeId == iMode)
                {
                    b = (-1.0) + b;
                }
                S[iMode] = b;
            }

            return S;
        }

        public System.Numerics.Complex CalcModeAmp(double omega, int modeIndex,
            System.Numerics.Complex[] betas,
            System.Numerics.Complex[][] ezEVecs, System.Numerics.Complex[][] ezXEVecs,
            System.Numerics.Complex[] Ez)
        {
            System.Numerics.Complex b = 0;

            double periodicDistance = WgPortInfo.PeriodicDistanceX;

            int bcNodeCnt = WgPortInfo.BcNodess[0].Count;
            int modeCnt = betas.Length;
            if (modeIndex >= modeCnt)
            {
                return b;
            }

            {
                int iMode = modeIndex;
                var beta = betas[iMode];
                var ezEVec = ezEVecs[iMode];
                var ezXEVec = ezXEVecs[iMode];
                var RyyZ = new IvyFEM.Lapack.ComplexMatrix(RyyB1);
                System.Diagnostics.Debug.Assert(ezEVec.Length == bcNodeCnt);
                System.Diagnostics.Debug.Assert(ezXEVec.Length == bcNodeCnt);
                System.Diagnostics.Debug.Assert(RyyZ.RowLength == bcNodeCnt);
                System.Diagnostics.Debug.Assert(RyyZ.ColumnLength == bcNodeCnt);

                var betaPeriodic = ToBetaPeriodic(beta, periodicDistance);
                var ezEVecModify = new System.Numerics.Complex[bcNodeCnt];
                for (int i = 0; i < bcNodeCnt; i++)
                {
                    ezEVecModify[i] =
                        ezEVec[i] - ezXEVec[i] / (System.Numerics.Complex.ImaginaryOne * betaPeriodic);
                }
                var vec1 = RyyZ * IvyFEM.Lapack.Utils.Conjugate(ezEVecModify);
                double replacedMu0 = IsTMMode ? Constants.Ep0 : Constants.Mu0;
                System.Numerics.Complex work1 = IvyFEM.Lapack.Functions.zdotu(vec1, Ez);
                b = beta.Magnitude *
                    (System.Numerics.Complex.Conjugate(betaPeriodic) /
                    System.Numerics.Complex.Conjugate(beta)) *
                    (1.0 / (omega * replacedMu0)) * work1;
            }
            return b;
        }
    }
}