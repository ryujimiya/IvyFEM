using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class PCWaveguide2DModalABCZTDFEM : FEM
    {
        //---------------------------------
        // ※ポートの設定順序
        // ABC境界、参照(観測)面、励振源の順
        //---------------------------------

        public uint QuantityId { get; private set; } = 0;

        /// <summary>
        /// Newmarkのβ法
        /// </summary>
        private double NewmarkBeta = 1.0 / 4.0;
        /// <summary>
        /// 時刻ステップ数
        /// </summary>
        public int TimeLoopCnt { get; set; } = 0;
        /// <summary>
        /// 計算時刻インデックス
        /// </summary>
        public int TimeIndex { get; set; } = 0;
        /// <summary>
        /// 時間ステップ幅
        /// </summary>
        public double TimeDelta { get; set; } = 0.0;
        /// <summary>
        /// ガウシアンパルス？
        /// true: ガウシアンパルス
        /// false: 正弦波
        /// </summary>
        public bool IsGaussian { get; set; } = true;
        /// <summary>
        /// ガウシアンパルスの種別
        /// </summary>
        public GaussianType GaussianType { get; set; } = GaussianType.Normal;
        /// <summary>
        /// ガウシアンパルスの遅延時間
        /// </summary>
        public double GaussianT0 { get; set; } = 0.0;
        /// <summary>
        /// ガウシアンパルスの時間幅
        /// </summary>
        public double GaussianTp { get; set; } = 0.0;
        /// <summary>
        /// 励振源のモード分布の周波数
        /// </summary>
        public double SrcFrequency { get; set; } = 0.0;
        /// <summary>
        /// Sマトリクスを求める開始周波数
        /// </summary>
        public double StartFrequencyForSMatrix { get; set; } = 0.0;
        /// <summary>
        /// Sマトリクスを求める開始周波数
        /// </summary>
        public double EndFrequencyForSMatrix { get; set; } = double.MaxValue;

        public IList<PCWaveguidePortInfo> WgPortInfos { get; set; } = null;
        /// <summary>
        /// TMモード？
        /// </summary>
        public bool IsTMMode { get; set; } = false;

        /// <summary>
        /// 観測点の頂点ID
        /// </summary>
        public IList<uint> RefVIds { get; set; } = new List<uint>();
        /// <summary>
        /// 観測点ポート数
        /// </summary>
        public int RefPortCount { get; set; } = 0;

        /// <summary>
        /// [A]
        /// </summary>
        private IvyFEM.Linear.ComplexSparseMatrix A = null;
        /// <summary>
        /// {b}
        /// </summary>
        private System.Numerics.Complex[] B = null;
        /// <summary>
        /// 剛性行列
        /// </summary>
        private IvyFEM.Linear.ComplexSparseMatrix K = null;
        // Note: defaultのPresicionLowerLimitだと正しい結果が得られない
        /// <summary>
        /// 質量行列
        /// </summary>
        private IvyFEM.Linear.ComplexSparseMatrix M = null;
        // Note: defaultのPresicionLowerLimitだと正しい結果が得られない
        /// <summary>
        /// 境界質量行列リスト(ポート単位)  p∫NiNj dy
        /// </summary>
        private IList<IvyFEM.Lapack.DoubleMatrix> Qbs = null;
        /// <summary>
        /// 境界剛性行列リスト(ポート単位)  p∫dNi/dy dNj/dy dy
        /// </summary>
        private IList<IvyFEM.Lapack.DoubleMatrix> Rbs = null;
        /// <summary>
        /// 境界質量行列リスト(ポート単位) q∫NiNj dy
        /// </summary>
        private IList<IvyFEM.Lapack.DoubleMatrix> Tbs = null;
        /// <summary>
        /// 境界積分
        /// </summary>
        private IList<IvyFEM.Lapack.ComplexMatrix> Bbs = null;

        /// <summary>
        /// 境界の界の伝搬定数(ポート単位)
        /// </summary>
        public IList<double> SrcBetaXs { get; private set; } = null;
        /// <summary>
        /// 境界の界のモード分布(ポート単位)
        /// </summary>
        public IList<System.Numerics.Complex[]> SrcProfiles { get; private set; } = null;
        /// <summary>
        /// 境界の界（周期構造用に修正)のモード分布(ポート単位)
        /// </summary>
        public IList<System.Numerics.Complex[]> SrcModifyProfiles { get; private set; } = null;

        /// <summary>
        /// 電界（現在値)
        /// </summary>
        public double[] Ez { get; private set; } = null;
        /// <summary>
        /// 電界（現在値)
        /// </summary>
        private System.Numerics.Complex[] EzZ = null;
        /// <summary>
        /// 電界(1つ前)
        /// </summary>
        private System.Numerics.Complex[] EzZPrev = null;
        /// <summary>
        /// 電界(2つ前)
        /// </summary>
        private System.Numerics.Complex[] EzZPrev2 = null;
        /// <summary>
        /// 観測点の電界(時間変化リスト)
        /// </summary>
        public IList<IList<double[]>> RefTimeEzsss { get; private set; } = new List<IList<double[]>>();

        public PCWaveguide2DModalABCZTDFEM(FEWorld world)
        {
            World = world;
        }

        public override void Solve()
        {
            int t;

            if (TimeIndex == -1)
            {
                return;
            }
            if (TimeIndex >= TimeLoopCnt)
            {
                TimeIndex = -1;
                return;
            }
            if (TimeIndex == 0)
            {
                int cnt = 0;
                if (RefVIds.Count > 0)
                {
                    cnt = RefVIds.Count;
                }
                else if (RefPortCount > 0)
                {
                    cnt = RefPortCount;
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
                RefTimeEzsss.Clear();
                for (int i = 0; i < cnt; i++)
                {
                    var refTimeEzss = new List<double[]>();
                    RefTimeEzsss.Add(refTimeEzss);
                }
            }

            if (TimeIndex == 0)
            {
                //--------------------------------------------------------------
                // 全体行列
                //--------------------------------------------------------------
                //t = System.Environment.TickCount;
                CalcA();
                //System.Diagnostics.Debug.WriteLine("CalcA t = " + (System.Environment.TickCount - t));
            }

            int nodeCnt = A.RowLength;

            //--------------------------------------------------------------
            // 電界
            //--------------------------------------------------------------
            // 以下TimeIndexに対する計算
            double time = TimeIndex * TimeDelta;
            double dt = TimeDelta;
            EzZPrev.CopyTo(EzZPrev2, 0);
            EzZ.CopyTo(EzZPrev, 0);
            EzZ = null;

            //--------------------------------------------------------------
            // 残差
            //--------------------------------------------------------------
            //t = System.Environment.TickCount;
            CalcB();
            //System.Diagnostics.Debug.WriteLine("CalcB t = " + (System.Environment.TickCount - t));

            //------------------------------------------------------------------
            // Ezを求める
            //------------------------------------------------------------------
            {
                System.Numerics.Complex[] X;
                Solver.ComplexSolve(out X, A, B);
                EzZ = X;
            }

            // 電界の実部を取得
            Ez = new double[nodeCnt];
            for (int nodeId = 0; nodeId < nodeCnt; nodeId++)
            {
                Ez[nodeId] = EzZ[nodeId].Real;
            }

            if (RefVIds.Count > 0)
            {
                // 観測点
                for (int refIndex = 0; refIndex < RefVIds.Count; refIndex++)
                {
                    uint vId = RefVIds[refIndex];
                    IList<int> coIds = World.GetCoordIdsFromCadId(QuantityId, vId, CadElementType.Vertex);
                    System.Diagnostics.Debug.Assert(coIds.Count == 1);
                    int coId = coIds[0];
                    int nodeId = World.Coord2Node(QuantityId, coId);
                    double fValue = Ez[nodeId];
                    double[] fValues = new double[] { fValue };
                    var refTimeEzss = RefTimeEzsss[refIndex];
                    refTimeEzss.Add(fValues);
                }
            }
            else if (RefPortCount > 0)
            {
                int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1; // 参照面と励振源を除く
                for (int refIndex = 0; refIndex < RefPortCount; refIndex++)
                {
                    int portId = portCnt + refIndex;
                    int nodeCntB;
                    {
                        var wgPortInfo = WgPortInfos[portId];
                        IList<int> bcNodes = wgPortInfo.BcNodess[0]; // 境界1
                        nodeCntB = bcNodes.Count;
                    }
                    double[] fValues = new double[nodeCntB];
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);
                        fValues[nodeIdB] = Ez[nodeId];
                    }
                    var refTimeEzss = RefTimeEzsss[refIndex];
                    refTimeEzss.Add(fValues);
                }
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }

            if ((TimeIndex + 1) % 50 == 0)
            {
                System.Diagnostics.Debug.WriteLine("timeIndex: {0}", (TimeIndex + 1));
            }
        }

        private void CalcA()
        {
            System.Diagnostics.Debug.Assert(TimeIndex == 0);

            IvyFEM.Linear.ComplexSparseMatrix _A;
            _A = null;
            Qbs = new List<IvyFEM.Lapack.DoubleMatrix>();
            Rbs = new List<IvyFEM.Lapack.DoubleMatrix>();
            Tbs = new List<IvyFEM.Lapack.DoubleMatrix>();
            Bbs = new List<IvyFEM.Lapack.ComplexMatrix>();
            SrcBetaXs = new List<double>();
            SrcProfiles = new List<System.Numerics.Complex[]>();
            SrcModifyProfiles = new List<System.Numerics.Complex[]>();
            EzZ = null;
            EzZPrev = null;
            EzZPrev2 = null;

            //------------------------------------------------------
            // 剛性行列、質量行列を作成
            //------------------------------------------------------
            CalcKM();

            //------------------------------------------------------
            // モード分布計算
            //------------------------------------------------------
            // 周波数
            double srcFreq = SrcFrequency;
            // 角周波数
            double srcOmega = 2.0 * Math.PI * srcFreq;
            // 波長
            double srcWaveLength = Constants.C0 / srcFreq;
            // 波数
            double srcK0 = 2.0 * Math.PI / srcWaveLength;

            int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1; // 参照面と励振源を除く
            var portConditions = World.GetPortConditions(QuantityId);

            for (int portId = 0; portId < (portCnt + RefPortCount + 1); portId++)
            {
                // 周期構造導波路
                var wgPortInfo = WgPortInfos[portId];
                wgPortInfo.PrevModeEVecs = null; // モード追跡初期化
                IvyFEM.Lapack.DoubleMatrix ryy1D;
                IvyFEM.Lapack.DoubleMatrix txx1D;
                IvyFEM.Lapack.DoubleMatrix uzz1D;
                System.Numerics.Complex[] betas;
                System.Numerics.Complex[][] eVecs;
                System.Numerics.Complex[][] fxEVecs;
                PCWaveguide2DEigenFEM eigenFEM;
                CalcEigen(
                    portId, srcFreq,
                    out ryy1D, out txx1D, out uzz1D, out betas, out eVecs, out fxEVecs, out eigenFEM);

                int nodeCntB = ryy1D.RowLength;
                Qbs.Add(ryy1D);
                Rbs.Add(txx1D);
                Tbs.Add(uzz1D);

                System.Diagnostics.Debug.WriteLine("port = {0} mode Count = {1}", portId, betas.Length);
                // 基本モード
                int iMode = 0;
                System.Numerics.Complex beta = betas[iMode];
                System.Numerics.Complex[] fVec = eVecs[iMode];
                System.Numerics.Complex[] fxVec = fxEVecs[iMode];
                System.Numerics.Complex[] fVecModify = new System.Numerics.Complex[nodeCntB];
                System.Diagnostics.Debug.Assert(fVec.Length == nodeCntB);
                for (int i = 0; i < nodeCntB; i++)
                {
                    fVecModify[i] =
                        fVec[i] - fxVec[i] / (System.Numerics.Complex.ImaginaryOne * beta);
                }
                // 実数部を取得する
                double betaReal = beta.Real;
                SrcBetaXs.Add(betaReal);
                SrcProfiles.Add(fVec);
                SrcModifyProfiles.Add(fVecModify);

                {
                    System.Numerics.Complex[] workBetas = { beta };
                    System.Numerics.Complex[][] workFVecs = { fVec };
                    System.Numerics.Complex[][] workFxVecs = { fxVec };
                    IvyFEM.Lapack.ComplexMatrix Bb =
                        eigenFEM.CalcBoundaryMatrix(srcOmega, workBetas, workFVecs, workFxVecs);
                    Bbs.Add(Bb);
                }
            }

            /////////////////////////////////////////////////////////
            //------------------------------------------------------
            // 節点数
            //-----------------------------------------------------
            int nodeCnt = (int)World.GetNodeCount(QuantityId);

            //------------------------------------------------------
            // 電界
            //-----------------------------------------------------
            EzZ = new System.Numerics.Complex[nodeCnt];
            EzZPrev = new System.Numerics.Complex[nodeCnt];
            EzZPrev2 = new System.Numerics.Complex[nodeCnt];

            //------------------------------------------------------
            // 全体係数行列の作成
            //------------------------------------------------------
            _A = new IvyFEM.Linear.ComplexSparseMatrix(nodeCnt, nodeCnt);
            double dt = TimeDelta;
            for (int rowNodeId = 0; rowNodeId < nodeCnt; rowNodeId++)
            {
                for (int colNodeId = 0; colNodeId < nodeCnt; colNodeId++)
                {
                    _A[rowNodeId, colNodeId] =
                        (1.0 / (dt * dt)) * M[rowNodeId, colNodeId] +
                        NewmarkBeta * K[rowNodeId, colNodeId];
                }
            }

            ////////////////////////////////////////////////////////////////////////////////////////
            // 吸収境界
            // Ez
            for (int portId = 0; portId < portCnt; portId++)
            {
                var Bb = Bbs[portId];
                int nodeCntB = Bb.RowLength;
                // Ez - Ez
                for (int rowNodeIdB = 0; rowNodeIdB < nodeCntB; rowNodeIdB++)
                {
                    // Ez
                    int rowCoId = World.PortNode2Coord(QuantityId, (uint)portId, rowNodeIdB);
                    int rowNodeId = World.Coord2Node(QuantityId, rowCoId);

                    for (int colNodeIdB = 0; colNodeIdB < nodeCntB; colNodeIdB++)
                    {
                        // Ez
                        int colCoId = World.PortNode2Coord(QuantityId, (uint)portId, colNodeIdB);
                        int colNodeId = World.Coord2Node(QuantityId, colCoId);
                        System.Numerics.Complex BbVal = Bb[rowNodeIdB, colNodeIdB];
                        _A[rowNodeId, colNodeId] +=
                            (1.0 / (2.0 * dt)) * BbVal / (System.Numerics.Complex.ImaginaryOne * srcOmega);
                    }
                }
            }

            A = _A;
        }

        private void CalcKM()
        {
            int nodeCnt = (int)World.GetNodeCount(QuantityId);
            K = new IvyFEM.Linear.ComplexSparseMatrix(nodeCnt, nodeCnt);
            K.PrecisionLowerLimit = 0.0;
            M = new IvyFEM.Linear.ComplexSparseMatrix(nodeCnt, nodeCnt);
            M.PrecisionLowerLimit = 0.0;

            IList<uint> feIds = World.GetTriangleFEIds(QuantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = World.GetTriangleFE(QuantityId, feId);
                uint elemNodeCnt = triFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = triFE.NodeCoordIds[iNode];
                    int nodeId = World.Coord2Node(QuantityId, coId);
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
                double[,] sNyNx = sNuNv[1, 0];
                double[,] sNxNy = sNuNv[0, 1];
                double[,] sNyNy = sNuNv[1, 1];
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
                        double kValue = maPxx * sNyNy[row, col] + maPyy * sNxNx[row, col];
                        // 要素質量行列
                        double mValue = Constants.Ep0 * Constants.Mu0 * maQzz * sNN[row, col];

                        K[rowNodeId, colNodeId] += kValue;
                        M[rowNodeId, colNodeId] += mValue;
                    }
                }
            }
        }

        private void CalcEigen(
            int portId, double srcFreq,
            out IvyFEM.Lapack.DoubleMatrix ryy1D,
            out IvyFEM.Lapack.DoubleMatrix txx1D,
            out IvyFEM.Lapack.DoubleMatrix uzz1D,
            out System.Numerics.Complex[] betas,
            out System.Numerics.Complex[][] bcEVecs,
            out System.Numerics.Complex[][] bcFxEVecs,
            out PCWaveguide2DEigenFEM eigenFEM)
        {
            var wgPortInfo = WgPortInfos[portId];
            eigenFEM = new PCWaveguide2DEigenFEM(World, QuantityId, (uint)portId, wgPortInfo);
            eigenFEM.IsTMMode = IsTMMode;
            eigenFEM.Frequency = srcFreq;
            eigenFEM.Solve();
            ryy1D = eigenFEM.RyyB1;
            txx1D = eigenFEM.TxxB1;
            uzz1D = eigenFEM.UzzB1;
            betas = eigenFEM.Betas;
            bcEVecs = eigenFEM.BcEVecs;
            bcFxEVecs = eigenFEM.BcFxEVecs;
        }

        private void CalcB()
        {
            double dt = TimeDelta;
            // 周波数
            double srcFreq = SrcFrequency;
            // 角周波数
            double srcOmega = 2.0 * Math.PI * srcFreq;

            int nodeCnt = A.RowLength;
            B = new System.Numerics.Complex[nodeCnt];
            {
                // nodeCntまで
                System.Numerics.Complex[] vecM = new System.Numerics.Complex[nodeCnt];
                System.Numerics.Complex[] vecK = new System.Numerics.Complex[nodeCnt];
                for (int i = 0; i < nodeCnt; i++)
                {
                    // M
                    vecM[i] = (2.0 / (dt * dt)) * EzZPrev[i] - (1.0 / (dt * dt)) * EzZPrev2[i];
                    // K
                    vecK[i] = -(1.0 - 2.0 * NewmarkBeta) * EzZPrev[i] - NewmarkBeta * EzZPrev2[i];
                }
                vecM = M * vecM;
                vecK = K * vecK;
                for (int i = 0; i < nodeCnt; i++)
                {
                    System.Numerics.Complex value = vecM[i] + vecK[i];
                    B[i] = value;
                }
            }

            int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1;  // 参照面、励振源分引く

            ////////////////////////////////////////////////////////////////////////////////
            // 吸収境界            
            // Ez
            for (int portId = 0; portId < portCnt; portId++)
            {
                var Bb = Bbs[portId];
                int nodeCntB = Bb.RowLength;

                System.Numerics.Complex[] workEzPrev = new System.Numerics.Complex[nodeCntB];
                System.Numerics.Complex[] workEzPrev2 = new System.Numerics.Complex[nodeCntB];
                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                {
                    // Ez
                    int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                    int nodeId = World.Coord2Node(QuantityId, coId);
                    workEzPrev[nodeIdB] = EzZPrev[nodeId];
                    workEzPrev2[nodeIdB] = EzZPrev2[nodeId];
                }

                System.Numerics.Complex[] vec = new System.Numerics.Complex[nodeCntB];
                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                {
                    for (int kNodeIdB = 0; kNodeIdB < nodeCntB; kNodeIdB++)
                    {
                        System.Numerics.Complex BbVal = Bb[nodeIdB, kNodeIdB];
                        vec[nodeIdB] +=
                            (1.0 / (2.0 * dt)) * workEzPrev2[kNodeIdB] *
                            BbVal / (System.Numerics.Complex.ImaginaryOne * srcOmega);
                    }
                }

                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                {
                    int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                    int nodeId = World.Coord2Node(QuantityId, coId);
                    B[nodeId] += vec[nodeIdB];
                }
            }

            //--------------------------------------------------------------
            // 励振源
            //--------------------------------------------------------------
            {
                int portId = portCnt + RefPortCount; // ポートリストの最後の要素が励振境界
                var Qb = Qbs[portId];
                int nodeCntB = Qb.RowLength;
                double srcBetaX = SrcBetaXs[portId];
                double vpx = srcOmega / srcBetaX;
                //System.Numerics.Complex[] srcProfile = SrcProfiles[portId];
                System.Numerics.Complex[] srcModifyProfile = SrcModifyProfiles[portId];

                System.Numerics.Complex srcU0 = 0.0;
                System.Numerics.Complex srcU1 = 0.0;
                System.Numerics.Complex srcU2 = 0.0;

                int n = TimeIndex;
                if (IsGaussian)
                {
                    if (GaussianType == GaussianType.Normal)
                    {
                        // ガウシアンパルス
                        if ((n * dt) <= (2.0 * GaussianT0 + Constants.PrecisionLowerLimit))
                        {
                            srcU0 = Math.Exp(-1.0 * ((n) * dt - GaussianT0) * ((n) * dt - GaussianT0) /
                                (2.0 * GaussianTp * GaussianTp));
                            srcU1 = Math.Exp(-1.0 * ((n - 1) * dt - GaussianT0) * ((n - 1) * dt - GaussianT0) /
                                (2.0 * GaussianTp * GaussianTp));
                            srcU2 = Math.Exp(-1.0 * ((n + 1) * dt - GaussianT0) * ((n + 1) * dt - GaussianT0) /
                                (2.0 * GaussianTp * GaussianTp));
                        }
                    }
                    else if (GaussianType == GaussianType.SinModulation)
                    {
                        // 正弦波変調ガウシアンパルス
                        if ((n * dt) <= (2.0 * GaussianT0 + Constants.PrecisionLowerLimit))
                        {
                            srcU0 = System.Numerics.Complex.Exp(
                                System.Numerics.Complex.ImaginaryOne * srcOmega * (n) * dt) *
                                Math.Exp(-1.0 * ((n) * dt - GaussianT0) * ((n) * dt - GaussianT0) /
                                (2.0 * GaussianTp * GaussianTp));
                            srcU1 = System.Numerics.Complex.Exp(
                                System.Numerics.Complex.ImaginaryOne * srcOmega * (n - 1) * dt) *
                                Math.Exp(-1.0 * ((n - 1) * dt - GaussianT0) * ((n - 1) * dt - GaussianT0) /
                                (2.0 * GaussianTp * GaussianTp));
                            srcU2 = System.Numerics.Complex.Exp(
                                System.Numerics.Complex.ImaginaryOne * srcOmega * (n + 1) * dt) *
                                Math.Exp(-1.0 * ((n + 1) * dt - GaussianT0) * ((n + 1) * dt - GaussianT0) /
                                (2.0 * GaussianTp * GaussianTp));
                        }
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                }
                else
                {
                    // 正弦波
                    srcU0 = System.Numerics.Complex.Exp(
                                System.Numerics.Complex.ImaginaryOne * srcOmega * (n) * dt);
                    srcU1 = System.Numerics.Complex.Exp(
                                System.Numerics.Complex.ImaginaryOne * srcOmega * (n - 1) * dt);
                    srcU2 = System.Numerics.Complex.Exp(
                                System.Numerics.Complex.ImaginaryOne * srcOmega * (n + 1) * dt);
                }

                {
                    // dEzinc/dt
                    System.Numerics.Complex[] srcUt = new System.Numerics.Complex[nodeCntB];
                    for (int i = 0; i < nodeCntB; i++)
                    {
                        srcUt[i] = (2.0 / vpx) * (srcModifyProfile[i].Real * (srcU2 - srcU1) / (2.0 * dt) +
                            srcModifyProfile[i].Imaginary * (1.0 / srcOmega) * (1.0 / (dt * dt)) *
                            (srcU2 - 2.0 * srcU0 + srcU1));
                    }
                    var QbZ = new IvyFEM.Lapack.ComplexMatrix(Qb);
                    System.Numerics.Complex[] vecQb = QbZ * srcUt;
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);
                        B[nodeId] += vecQb[nodeIdB];
                    }
                }
                /*
                {
                    // dEzinc/dt
                    System.Numerics.Complex[] srcUt = new System.Numerics.Complex[nodeCntB];
                    for (int i = 0; i < nodeCntB; i++)
                    {
                        srcUt[i] = (2.0 / vpx) * srcModifyProfile[i] * (srcU2 - srcU1) / (2.0 * dt);
                    }
                    var QbZ = new IvyFEM.Lapack.ComplexMatrix(Qb);
                    System.Numerics.Complex[] vecQb = QbZ * srcUt;
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);
                        B[nodeId] += vecQb[nodeIdB];
                    }
                }
                */
            }
        }

        /// <summary>
        /// Sパラメータの計算 freqs[iFreq] Sss[portId][iFreq]
        /// </summary>
        /// <returns></returns>
        public void CalcSParameter(
            System.Numerics.Complex[] freqDomainAmpsInc,
            out double[] freqs,
            out IList<System.Numerics.Complex[]> freqDomainAmpss,
            out IList<System.Numerics.Complex[]> Sss)
        {
            freqs = null;
            freqDomainAmpss = null;
            Sss = null;

            if (RefVIds.Count > 0)
            {
                CalcSParameterAtRefPoints(freqDomainAmpsInc, out freqs, out freqDomainAmpss, out Sss);
            }
            else if (RefPortCount > 0)
            {
                CalcSParameterAtRefPorts(freqDomainAmpsInc, out freqs, out freqDomainAmpss, out Sss);
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
        }

        private void CalcSParameterAtRefPoints(
            System.Numerics.Complex[] freqDomainAmpsInc,
            out double[] freqs,
            out IList<System.Numerics.Complex[]> freqDomainAmpss,
            out IList<System.Numerics.Complex[]> Sss)
        {
            System.Diagnostics.Debug.Assert(RefVIds.Count > 0);
            double[][] datasRefs = new double[RefVIds.Count][];
            for (int refIndex = 0; refIndex < RefVIds.Count; refIndex++)
            {
                int timeCnt = RefTimeEzsss[refIndex].Count;
                double[] datasRef = new double[timeCnt];
                for (int timeIndex = 0; timeIndex < timeCnt; timeIndex++)
                {
                    System.Diagnostics.Debug.Assert(RefTimeEzsss[refIndex][timeIndex].Length == 1);
                    datasRef[timeIndex] = RefTimeEzsss[refIndex][timeIndex][0];
                }
                datasRefs[refIndex] = datasRef;
            }

            int dataCnt = datasRefs[0].Length;
            double dt = TimeDelta;

            double[] times = new double[dataCnt];
            for (int i = 0; i < dataCnt; i++)
            {
                times[i] = i * dt;
            }

            ////////////////////////////////////////////////////////////////////
            // FFT
            freqs = null;
            freqDomainAmpss = new List<System.Numerics.Complex[]>();
            for (int refIndex = 0; refIndex < RefVIds.Count; refIndex++)
            {
                double[] datas = datasRefs[refIndex];
                double[] _freqs;
                System.Numerics.Complex[] freqDomainAmps;
                IvyFEM.FFT.Functions.DoFFT(times, datas, out _freqs, out freqDomainAmps);
                if (refIndex == 0)
                {
                    freqs = _freqs;
                }
                freqDomainAmpss.Add(freqDomainAmps);
            }

            System.Diagnostics.Debug.Assert(freqs != null && freqs.Length == dataCnt);
            System.Numerics.Complex[] _freqDomainAmpsInc =
                freqDomainAmpsInc != null ? freqDomainAmpsInc : freqDomainAmpss[0];

            // 散乱パラメータを計算する
            Sss = new List<System.Numerics.Complex[]>();
            for (int refIndex = 0; refIndex < RefVIds.Count; refIndex++)
            {
                System.Numerics.Complex[] Sp1s = new System.Numerics.Complex[dataCnt];
                // Sss: Sss[port][freq]
                Sss.Add(Sp1s);
                System.Numerics.Complex[] freqDomainDatas = freqDomainAmpss[refIndex];

                for (int i = 0; i < dataCnt; i++)
                {
                    System.Numerics.Complex dataInc = _freqDomainAmpsInc[i];
                    System.Numerics.Complex data = freqDomainDatas[i];
                    System.Numerics.Complex Sp1 = 0;
                    if (refIndex == 0)
                    {
                        Sp1 = (data - dataInc) / dataInc;
                    }
                    else
                    {
                        Sp1 = data / dataInc;
                    }

                    Sp1s[i] = Sp1;
                }
            }
        }

        private void CalcSParameterAtRefPorts(
            System.Numerics.Complex[] freqDomainAmpsInc,
            out double[] freqs,
            out IList<System.Numerics.Complex[]> freqDomainAmpss,
            out IList<System.Numerics.Complex[]> Sss)
        {
            System.Diagnostics.Debug.Assert(RefPortCount > 0);
            double[][][] datasRefs = new double[RefPortCount][][];
            for (int refIndex = 0; refIndex < RefPortCount; refIndex++)
            {
                int timeCnt = RefTimeEzsss[refIndex].Count;
                int nodeCntB = RefTimeEzsss[refIndex][0].Length;
                double[][] datasRef = new double[nodeCntB][];
                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                {
                    datasRef[nodeIdB] = new double[timeCnt];
                    for (int timeIndex = 0; timeIndex < timeCnt; timeIndex++)
                    {
                        datasRef[nodeIdB][timeIndex] = RefTimeEzsss[refIndex][timeIndex][nodeIdB];
                    }
                }
                datasRefs[refIndex] = datasRef;
            }

            int dataCnt = datasRefs[0][0].Length;
            double dt = TimeDelta;

            double[] times = new double[dataCnt];
            for (int i = 0; i < dataCnt; i++)
            {
                times[i] = i * dt;
            }

            ////////////////////////////////////////////////////////////////////
            // FFT
            freqs = null;
            freqDomainAmpss = new List<System.Numerics.Complex[]>();

            for (int refIndex = 0; refIndex < RefPortCount; refIndex++)
            {
                double[][] datass = datasRefs[refIndex];
                int nodeCntB = datass.Length;
                System.Diagnostics.Debug.Assert(datass[0].Length == dataCnt);
                System.Numerics.Complex[][] freqDomainDatass = new System.Numerics.Complex[dataCnt][];
                for (int freqIndex = 0; freqIndex < dataCnt; freqIndex++)
                {
                    freqDomainDatass[freqIndex] = new System.Numerics.Complex[nodeCntB];
                }
                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                {
                    double[] datas = datass[nodeIdB];
                    double[] _freqs;
                    System.Numerics.Complex[] _freqDomainDatas;
                    IvyFEM.FFT.Functions.DoFFT(times, datas, out _freqs, out _freqDomainDatas);
                    System.Diagnostics.Debug.Assert(_freqs.Length == dataCnt);
                    if (refIndex == 0 && nodeIdB == 0)
                    {
                        freqs = _freqs;
                    }
                    for (int freqIndex = 0; freqIndex < dataCnt; freqIndex++)
                    {
                        System.Numerics.Complex data = _freqDomainDatas[freqIndex];
                        freqDomainDatass[freqIndex][nodeIdB] = data;
                    }
                }

                // モード振幅の計算
                int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1; // 参照面、励振源を引く
                System.Numerics.Complex[] freqDomainAmps = new System.Numerics.Complex[dataCnt];
                for (int freqIndex = 0; freqIndex < dataCnt; freqIndex++)
                {
                    // 周波数
                    double freq = freqs[freqIndex];
                    // 角周波数
                    double omega = 2.0 * Math.PI * freq;

                    if (freq < StartFrequencyForSMatrix || freq > EndFrequencyForSMatrix)
                    {
                        freqDomainAmps[freqIndex] = 0;
                        continue;
                    }

                    int portId = portCnt + refIndex;
                    var wgPortInfo = WgPortInfos[portId];
                    if (freqIndex == 0)
                    {
                        wgPortInfo.PrevModeEVecs = null; // モード追跡初期化
                    }
                    // モード
                    IvyFEM.Lapack.DoubleMatrix ryy1D;
                    IvyFEM.Lapack.DoubleMatrix txx1D;
                    IvyFEM.Lapack.DoubleMatrix uzz1D;
                    System.Numerics.Complex[] betas;
                    System.Numerics.Complex[][] eVecs;
                    System.Numerics.Complex[][] fxEVecs;
                    PCWaveguide2DEigenFEM eigenFEM;
                    CalcEigen(
                        portId, freq,
                        out ryy1D, out txx1D, out uzz1D, out betas, out eVecs, out fxEVecs, out eigenFEM);

                    // 振幅分布
                    System.Numerics.Complex[] freqEz = freqDomainDatass[freqIndex];
                    // モード振幅の算出
                    int iMode = 0; // 基本モード
                    System.Numerics.Complex b = eigenFEM.CalcModeAmp(omega, iMode, betas, eVecs, fxEVecs, freqEz);
                    freqDomainAmps[freqIndex] = b;
                }
                freqDomainAmpss.Add(freqDomainAmps);
            }

            System.Diagnostics.Debug.Assert(freqs != null && freqs.Length == dataCnt);
            System.Numerics.Complex[] _freqDomainAmpsInc =
                freqDomainAmpsInc != null ? freqDomainAmpsInc : freqDomainAmpss[0];

            // 散乱パラメータを計算する
            Sss = new List<System.Numerics.Complex[]>();
            for (int refIndex = 0; refIndex < RefPortCount; refIndex++)
            {
                System.Numerics.Complex[] Sp1s = new System.Numerics.Complex[dataCnt];
                // Sss: Sss[port][freq]
                Sss.Add(Sp1s);
                System.Numerics.Complex[] freqDomainDatas = freqDomainAmpss[refIndex];

                for (int i = 0; i < dataCnt; i++)
                {
                    System.Numerics.Complex dataInc = _freqDomainAmpsInc[i];
                    System.Numerics.Complex data = freqDomainDatas[i];
                    System.Numerics.Complex Sp1 = 0;
                    if (refIndex == 0)
                    {
                        Sp1 = (data - dataInc) / dataInc;
                    }
                    else
                    {
                        Sp1 = data / dataInc;
                    }

                    Sp1s[i] = Sp1;
                }
            }
        }
    }
}
