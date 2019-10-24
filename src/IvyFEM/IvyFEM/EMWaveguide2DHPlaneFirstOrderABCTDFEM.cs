using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class EMWaveguide2DHPlaneFirstOrderABCTDFEM : FEM
    {
        //---------------------------------
        // ※ポートの設定順序
        // ABC境界、参照(観測)面、励振源の順
        //---------------------------------

        public uint QuantityId { get; private set; } = 0;
        /// <summary>
        /// 吸収境界条件の次数
        /// </summary>
        public IList<int> ABCOrdersToSet { get; set; } = new List<int>();
        /// <summary>
        /// Evanescent Waveの吸収境界条件の次数(ポート単位)
        /// </summary>
        public IList<int> ABCOrdersForEvanescentToSet { get; set; } = new List<int>();
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
        /// 励振を点波源で行う？
        /// </summary>
        public bool IsPointExcitation { get; set; } = false;
        /// <summary>
        /// ABCで使う速度 (default: リスト無し(defaultの速度), -1を指定してもdefaultの速度になる)
        /// </summary>
        public IList<double> VelocitysToSet { get; set; } = new List<double>();
        /// <summary>
        /// 1D固有値問題で減衰定数を用いる？
        /// </summary>
        public IList<bool> IsEigen1DUseDecayParameters { get; set; } = new List<bool>(); 
        /// <summary>
        /// 1D固有値問題のクラッド比誘電率
        /// </summary>
        public IList<double> Eigen1DCladdingEps { get; set; } = new List<double>();
        /// <summary>
        /// 減衰定数を持ってくる1D固有値問題のポート
        /// </summary>
        public IList<int> DecayParameterEigen1DPortIds { get; set; } = new List<int>();

        /// <summary>
        /// TEモードで実装した式をTMモードに流用するため
        ///   TEモードの場合は μ0
        ///   TMモードの場合は ε0
        /// </summary>
        public double ReplacedMu0 { get; set; } = Constants.Mu0;
        /// <summary>
        /// 観測点の頂点ID
        /// </summary>
        public IList<uint> RefVIds { get; set; } = new List<uint>();
        /// <summary>
        /// 観測点ポート数
        /// </summary>
        public int RefPortCount { get; set; } = 0;

        /// <summary>
        /// 逆行列を用いる？
        /// </summary>
        public bool IsUseInvMatrix { get; set; } = true;
        /// <summary>
        /// [A]
        /// </summary>
        private IvyFEM.Linear.DoubleSparseMatrix A = null;
        /// <summary>
        /// {b}
        /// </summary>
        private double[] B = null;
        /// <summary>
        /// 剛性行列
        /// </summary>
        private IvyFEM.Lapack.DoubleMatrix K = null;
        // Note: IvyFEM.Linear.DoubleSparseMatrixだと正しい結果が得られない
        /// <summary>
        /// 質量行列
        /// </summary>
        private IvyFEM.Lapack.DoubleMatrix M = null;
        // Note: IvyFEM.Linear.DoubleSparseMatrixだと正しい結果が得られない
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
        /// 境界の界の伝搬定数(ポート単位)
        /// </summary>
        public IList<double> SrcBetaXs { get; private set; } = null;
        /// <summary>
        /// 境界の界のモード分布(ポート単位)
        /// </summary>
        public IList<double[]> SrcProfiles { get; private set; } = null;
        /// <summary>
        /// 境界の界の両端の減衰定数
        /// </summary>
        public IList<double> SrcDecayParameters { get; private set; } = null;

        /// <summary>
        /// ABC(Evanescent)の次数
        /// </summary>
        private IList<int> ABCOrdersFor1 = null;
        /// <summary>
        /// ABCの次数
        /// </summary>
        private IList<int> ABCOrdersFor2 = null;

        /// <summary>
        /// Evanescent Waveの減衰定数
        /// </summary>
        private IList<double> AlphasFor1 = null;

        /// <summary>
        /// 吸収境界波の速度(ポート単位)
        /// </summary>
        private IList<double> VelosFor2 = null;

        /// <summary>
        /// Φ1(2)のオフセット(Traveling)
        /// </summary>
        private IList<int> Pz1Offsets = null;

        /// <summary>
        ///  電界（現在値）※ABC分の自由度を除いたもの
        /// </summary>
        public double[] Ez { get; private set; } = null;
        /// <summary>
        /// 電界（現在値)
        /// </summary>
        private double[] EzPz = null;
        /// <summary>
        /// 電界(1つ前)
        /// </summary>
        private double[] EzPzPrev = null;
        /// <summary>
        /// 電界(2つ前)
        /// </summary>
        private double[] EzPzPrev2 = null;
        /// <summary>
        /// 観測点の電界(時間変化リスト)
        /// </summary>
        public IList<IList<double[]>> RefTimeEzsss { get; private set; } = new List<IList<double[]>>();

        public EMWaveguide2DHPlaneFirstOrderABCTDFEM(FEWorld world)
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

            int nodeCnt = K.RowLength;
            int nodeCntPlusABC = A.RowLength;

            //--------------------------------------------------------------
            // 電界
            //--------------------------------------------------------------
            // 以下TimeIndexに対する計算
            double time = TimeIndex * TimeDelta;
            double dt = TimeDelta;
            EzPzPrev.CopyTo(EzPzPrev2, 0);
            EzPz.CopyTo(EzPzPrev, 0);
            EzPz = null;

            //--------------------------------------------------------------
            // 残差
            //--------------------------------------------------------------
            //t = System.Environment.TickCount;
            CalcB();
            //System.Diagnostics.Debug.WriteLine("CalcB t = " + (System.Environment.TickCount - t));

            //------------------------------------------------------------------
            // Ezを求める
            //------------------------------------------------------------------
            if (IsUseInvMatrix)
            {
                // 逆行列を用いる
                EzPz = A * B;
            }
            else
            {
                double[] X;
                Solver.DoubleSolve(out X, A, B);
                EzPz = X;
            }

            // ABC分を除いた電界を取得
            Ez = new double[nodeCnt];
            for (int nodeId = 0; nodeId < nodeCnt; nodeId++)
            {
                Ez[nodeId] = EzPz[nodeId];
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
                    int nodeCntB = (int)World.GetPortNodeCount(QuantityId, (uint)portId);
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

            IvyFEM.Lapack.DoubleMatrix _A;
            _A = null;
            Qbs = new List<IvyFEM.Lapack.DoubleMatrix>();
            Rbs = new List<IvyFEM.Lapack.DoubleMatrix>();
            Tbs = new List<IvyFEM.Lapack.DoubleMatrix>();
            SrcBetaXs = new List<double>();
            SrcProfiles = new List<double[]>();
            SrcDecayParameters = new List<double>();
            EzPz = null;
            EzPzPrev = null;
            EzPzPrev2 = null;

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

            System.Diagnostics.Debug.Assert(
                IsEigen1DUseDecayParameters.Count == 0 ||
                IsEigen1DUseDecayParameters.Count == (portCnt + RefPortCount + 1));
            for (int portId = 0; portId < (portCnt + RefPortCount + 1); portId++)
            {
                IvyFEM.Lapack.DoubleMatrix ryy1D;
                IvyFEM.Lapack.DoubleMatrix txx1D;
                IvyFEM.Lapack.DoubleMatrix uzz1D;
                System.Numerics.Complex[] betas;
                System.Numerics.Complex[][] eVecs;
                double alpha;
                EMWaveguide1DEigenBaseFEM eigen1DBaseFEM;
                CalcEigen(
                    portId, srcFreq,
                    out ryy1D, out txx1D, out uzz1D, out betas, out eVecs, out alpha, out eigen1DBaseFEM);

                int nodeCntB = ryy1D.RowLength;
                Qbs.Add(ryy1D);
                Rbs.Add(txx1D);
                Tbs.Add(uzz1D);

                System.Diagnostics.Debug.WriteLine("port = {0} mode Count = {1}", portId, betas.Length);
                // 基本モード
                int iMode = 0;
                System.Numerics.Complex beta = betas[iMode];
                System.Numerics.Complex[] fVec = eVecs[iMode];
                // 実数部を取得する
                double betaReal = beta.Real;
                double[] fVecReal = new double[nodeCntB];
                for (int i = 0; i < nodeCntB; i++)
                {
                    double value = fVec[i].Real;
                    fVecReal[i] = value;
                }
                SrcBetaXs.Add(betaReal);
                SrcProfiles.Add(fVecReal);
                SrcDecayParameters.Add(alpha);
            }

            /////////////////////////////////////////////////////////
            //------------------------------------------------------
            // ABC次数
            //------------------------------------------------------
            System.Diagnostics.Debug.Assert(
                ABCOrdersToSet.Count == 0 || ABCOrdersToSet.Count == portCnt);
            System.Diagnostics.Debug.Assert(
                ABCOrdersForEvanescentToSet.Count == 0 || ABCOrdersForEvanescentToSet.Count == portCnt);
            ABCOrdersFor1 = new List<int>();
            ABCOrdersFor2 = new List<int>();
            for (int portId = 0; portId < portCnt; portId++)
            {
                int abcOrderFor1 = 0;
                if (ABCOrdersForEvanescentToSet.Count == portCnt)
                {
                    abcOrderFor1 = ABCOrdersForEvanescentToSet[portId];
                }
                int abcOrderFor2 = 0;
                if (ABCOrdersToSet.Count == portCnt)
                {
                    abcOrderFor2 = ABCOrdersToSet[portId];
                }
                System.Diagnostics.Debug.Assert(abcOrderFor1 == 1 || abcOrderFor2 == 1);
                ABCOrdersFor1.Add(abcOrderFor1);
                ABCOrdersFor2.Add(abcOrderFor2);
            }
            //------------------------------------------------------
            // 電界
            //-----------------------------------------------------
            int nodeCnt = (int)World.GetNodeCount(QuantityId);
            int nodeCntPlusABC = nodeCnt;
            for (int portId = 0; portId < portCnt; portId++)
            {
                int nodeCntB = (int)World.GetPortNodeCount(QuantityId, (uint)portId);
                System.Diagnostics.Debug.Assert(Qbs[portId].RowLength == nodeCntB);
                int abcOrderFor1 = ABCOrdersFor1[portId];
                int abcOrderFor2 = ABCOrdersFor2[portId];
                if (abcOrderFor1 == 1 && abcOrderFor2 == 1)
                {
                    nodeCntPlusABC += nodeCntB;
                }
            }
            EzPz = new double[nodeCntPlusABC];
            EzPzPrev = new double[nodeCntPlusABC];
            EzPzPrev2 = new double[nodeCntPlusABC];

            //------------------------------------------------------
            // 全体係数行列の作成
            //------------------------------------------------------
            _A = new IvyFEM.Lapack.DoubleMatrix(nodeCntPlusABC, nodeCntPlusABC);
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
            // 吸収境界パラメータ
            AlphasFor1 = new List<double>();
            for (int portId = 0; portId < portCnt; portId++)
            {
                int portId0 = -1;
                System.Diagnostics.Debug.Assert(
                    DecayParameterEigen1DPortIds.Count == 0 ||
                    DecayParameterEigen1DPortIds.Count == portCnt);
                if (DecayParameterEigen1DPortIds.Count == portCnt)
                {
                    portId0 = DecayParameterEigen1DPortIds[portId];
                }
                else
                {
                    portId0 = portCnt; // 最後のポート(励振源)を採用
                }
                double alpha0 = 0;
                if (portId0 == -1)
                {
                    alpha0 = 0;
                }
                else
                {
                    alpha0 = SrcDecayParameters[portId0];
                }
                AlphasFor1.Add(alpha0);
            }

            VelosFor2 = new List<double>();
            for (int portId = 0; portId < portCnt; portId++)
            {
                double velo0 = -1.0;
                if (VelocitysToSet.Count == portCnt)
                {
                    velo0 = VelocitysToSet[portId];
                }
                if (velo0 < 0) // velo0: -1
                {
                    if (GaussianType == GaussianType.Normal)
                    {
                        velo0 = Constants.C0; // ガウシアンパルスの場合
                    }
                    else if (GaussianType == GaussianType.SinModulation)
                    {
                        double srcBetaX = SrcBetaXs[portId];
                        double[] srcProfile = SrcProfiles[portId];
                        double vpx = srcOmega / srcBetaX;
                        velo0 = vpx; // 正弦波変調ガウシアンパルスの場合
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                }
                VelosFor2.Add(velo0);
            }

            ////////////////////////////////////////////////////////////////////////////////////////
            // Φ1の開始位置
            Pz1Offsets = new List<int>();
            for (int portId = 0; portId < portCnt; portId++)
            {
                int offset = 0;
                if (portId == 0)
                {
                    offset = nodeCnt;
                }
                else
                {
                    int nodeCntBPrevPort = (int)World.GetPortNodeCount(QuantityId, (uint)(portId - 1));
                    int abcOrderFor1PrevPort = ABCOrdersFor1[portId - 1];
                    int abcOrderFor2PrevPort = ABCOrdersFor2[portId - 1];

                    offset = Pz1Offsets[portId - 1];
                    if (abcOrderFor1PrevPort == 1 && abcOrderFor2PrevPort == 1)
                    {
                        offset += nodeCntBPrevPort;
                    }
                }
                Pz1Offsets.Add(offset);
            }

            ////////////////////////////////////////////////////////////////////////////////////////
            // 吸収境界
            // Ez
            for (int portId = 0; portId < portCnt; portId++)
            {
                int abcOrderFor1 = ABCOrdersFor1[portId];
                int abcOrderFor2 = ABCOrdersFor2[portId];
                var Qb = Qbs[portId];
                int nodeCntB = Qb.RowLength;

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

                        if (abcOrderFor1 == 1 && abcOrderFor2 == 0)
                        {
                            // Evanescent
                            double alpha0 = AlphasFor1[portId];
                            _A[rowNodeId, colNodeId] +=
                                NewmarkBeta * alpha0  * Qb[rowNodeIdB, colNodeIdB];
                        }
                        else if (abcOrderFor2 == 1 && abcOrderFor1 == 0)
                        {
                            // Traveling
                            double velo0 = VelosFor2[portId];
                            _A[rowNodeId, colNodeId] +=
                                (1.0 / (2.0 * dt * velo0)) * Qb[rowNodeIdB, colNodeIdB];
                        }
                        else if (abcOrderFor1 == 1 && abcOrderFor2 == 1)
                        {
                            // Evanescent & Traveling
                        }
                        else
                        {
                            System.Diagnostics.Debug.Assert(false);
                        }
                    }
                }

                // Ez - Φ1
                for (int rowNodeIdB = 0; rowNodeIdB < nodeCntB; rowNodeIdB++)
                {
                    int rowCoId = World.PortNode2Coord(QuantityId, (uint)portId, rowNodeIdB);
                    int rowNodeId = World.Coord2Node(QuantityId, rowCoId);

                    if (abcOrderFor1 == 1 && abcOrderFor2 == 1)
                    {
                        for (int colNodeIdB = 0; colNodeIdB < nodeCntB; colNodeIdB++)
                        {
                            // Φ1
                            int colNodeId1 = colNodeIdB + Pz1Offsets[portId];
                            _A[rowNodeId, colNodeId1] += -NewmarkBeta * Qb[rowNodeIdB, colNodeIdB];
                        }
                    }
                }
            }

            // Evanescent & Traveling
            // Φ1
            for (int portId = 0; portId < portCnt; portId++)
            {
                var Qb = Qbs[portId];
                var Rb = Rbs[portId];
                var Tb = Tbs[portId];
                double alpha0 = AlphasFor1[portId];
                double velo0 = VelosFor2[portId];
                int nodeCntB = Qb.RowLength;
                int abcOrderFor1 = ABCOrdersFor1[portId];
                int abcOrderFor2 = ABCOrdersFor2[portId];

                if (abcOrderFor1 == 1 && abcOrderFor2 == 1)
                {
                    for (int rowNodeIdB = 0; rowNodeIdB < nodeCntB; rowNodeIdB++)
                    {
                        // Φ1
                        int rowNodeId1 = rowNodeIdB + Pz1Offsets[portId];

                        for (int colNodeIdB = 0; colNodeIdB < nodeCntB; colNodeIdB++)
                        {
                            int colNodeId1 = 0;
                            // Ez
                            int colCoId = World.PortNode2Coord(QuantityId, (uint)portId, colNodeIdB);
                            colNodeId1 = World.Coord2Node(QuantityId, colCoId);
                            _A[rowNodeId1, colNodeId1] +=
                                (alpha0 / velo0) * (1.0 / (2.0 * dt)) * Qb[rowNodeIdB, colNodeIdB] +
                                NewmarkBeta * Rb[rowNodeIdB, colNodeIdB] +
                                (1.0 / (dt * dt * velo0 * velo0)) * Tb[rowNodeIdB, colNodeIdB];

                            // Φ1
                            colNodeId1 = colNodeIdB + Pz1Offsets[portId];
                            _A[rowNodeId1, colNodeId1] +=
                                (NewmarkBeta * alpha0 +
                                1.0 / (2.0 * dt * velo0)) * Qb[rowNodeIdB, colNodeIdB];

                        }
                    }
                }
            }

            if (IsUseInvMatrix)
            {
                // 逆行列を計算
                System.Diagnostics.Debug.WriteLine("calc [A]-1");
                _A = IvyFEM.Lapack.DoubleMatrix.Inverse(_A);
                System.Diagnostics.Debug.WriteLine("calc [A]-1 done");

            }

            A = new IvyFEM.Linear.DoubleSparseMatrix(_A);
        }

        private void CalcKM()
        {
            int nodeCnt = (int)World.GetNodeCount(QuantityId);
            K = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            M = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);

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
                        double kValue = (1.0 / ma.Muxx) * sNyNy[row, col] + (1.0 / ma.Muyy) * sNxNx[row, col];
                        // 要素質量行列
                        double mValue = Constants.Ep0 * Constants.Mu0 * ma.Epzz * sNN[row, col];

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
            out System.Numerics.Complex[][] eVecs,
            out double alpha,
            out EMWaveguide1DEigenBaseFEM eigen1DBaseFEM)
        {
            int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1; // 参照面と励振源を除く
            alpha = 0;
            if (IsEigen1DUseDecayParameters.Count == (portCnt + RefPortCount + 1) &&
                IsEigen1DUseDecayParameters[portId])
            {
                // 減衰定数を考慮した固有値問題
                var eigen1DFEM = new EMWaveguide1DOpenEigenFEM(World, QuantityId, (uint)portId);
                eigen1DBaseFEM = eigen1DFEM;
                eigen1DFEM.CladdingEp = Eigen1DCladdingEps[portId];
                eigen1DFEM.ReplacedMu0 = ReplacedMu0;
                eigen1DFEM.Frequency = srcFreq;
                eigen1DFEM.Solve();
                ryy1D = eigen1DFEM.Ryy;
                txx1D = eigen1DFEM.Txx;
                uzz1D = eigen1DFEM.Uzz;
                betas = eigen1DFEM.Betas;
                eVecs = eigen1DFEM.EzEVecs;
                alpha = eigen1DFEM.DecayParameter;
            }
            else
            {
                // 通常の固有値問題
                var eigen1DFEM = new EMWaveguide1DEigenFEM(World, QuantityId, (uint)portId);
                eigen1DBaseFEM = eigen1DFEM;
                eigen1DFEM.ReplacedMu0 = ReplacedMu0;
                eigen1DFEM.Frequency = srcFreq;
                eigen1DFEM.Solve();
                ryy1D = eigen1DFEM.Ryy;
                txx1D = eigen1DFEM.Txx;
                uzz1D = eigen1DFEM.Uzz;
                betas = eigen1DFEM.Betas;
                eVecs = eigen1DFEM.EzEVecs;
            }
        }

        private void CalcB()
        {
            double dt = TimeDelta;
            int nodeCntPlusABC = A.RowLength;
            int nodeCnt = K.RowLength;
            B = new double[nodeCntPlusABC];
            {
                // nodeCntまで
                double[] vecM = new double[nodeCnt];
                double[] vecK = new double[nodeCnt];
                for (int i = 0; i < nodeCnt; i++)
                {
                    // M
                    vecM[i] = (2.0 / (dt * dt)) * EzPzPrev[i] - (1.0 / (dt * dt)) * EzPzPrev2[i];
                    // K
                    vecK[i] = -(1.0 - 2.0 * NewmarkBeta) * EzPzPrev[i] - NewmarkBeta * EzPzPrev2[i];
                }
                vecM = M * vecM;
                vecK = K * vecK;
                for (int i = 0; i < nodeCnt; i++)
                {
                    double value = vecM[i] + vecK[i];
                    B[i] = value;
                }
            }
            int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1;  // 参照面、励振源分引く
            ////////////////////////////////////////////////////////////////////////////////
            // 吸収境界            
            // Ez
            for (int portId = 0; portId < portCnt; portId++)
            {
                var Qb = Qbs[portId];
                int nodeCntB = Qb.RowLength;
                int abcOrderFor1 = ABCOrdersFor1[portId];
                int abcOrderFor2 = ABCOrdersFor2[portId];

                if (abcOrderFor1 == 1 && abcOrderFor2 == 0)
                {
                    // Evanescent
                    double alpha0 = AlphasFor1[portId];
                    double[] workEzPrev = new double[nodeCntB];
                    double[] workEzPrev2 = new double[nodeCntB];
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        // Ez
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);
                        workEzPrev[nodeIdB] = EzPzPrev[nodeId];
                        workEzPrev2[nodeIdB] = EzPzPrev2[nodeId];
                    }

                    double[] vecQb = new double[nodeCntB];
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        vecQb[nodeIdB] = -1.0 * alpha0 *
                            ((1.0 - 2.0 * NewmarkBeta) * workEzPrev[nodeIdB] +
                            NewmarkBeta * workEzPrev2[nodeIdB]);
                    }
                    vecQb = Qb * vecQb;

                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);
                        B[nodeId] += vecQb[nodeIdB];
                    }
                }
                else if (abcOrderFor2 == 1 && abcOrderFor1 == 0)
                {
                    // Traveling
                    double velo0 = VelosFor2[portId];
                    double[] workEzPrev2 = new double[nodeCntB];
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        // Ez
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);
                        workEzPrev2[nodeIdB] = EzPzPrev2[nodeId];
                    }

                    double[] vecQb = new double[nodeCntB];
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        vecQb[nodeIdB] = (1.0 / (2.0 * dt * velo0)) * workEzPrev2[nodeIdB];
                    }
                    vecQb = Qb * vecQb;

                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);
                        B[nodeId] += vecQb[nodeIdB];
                    }
                }
                else if (abcOrderFor1 == 1 && abcOrderFor2 == 1)
                {
                    // Traveling & Evanescent
                    double[] workPzPrev = new double[nodeCntB];
                    double[] workPzPrev2 = new double[nodeCntB];
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        // Pz
                        int nodeId = nodeIdB + Pz1Offsets[portId];
                        workPzPrev[nodeIdB] = EzPzPrev[nodeId];
                        workPzPrev2[nodeIdB] = EzPzPrev2[nodeId];
                    }

                    double[] vecQb = new double[nodeCntB];
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        vecQb[nodeIdB] =
                            ((1.0 - 2.0 * NewmarkBeta) * workPzPrev[nodeIdB] +
                            NewmarkBeta * workPzPrev2[nodeIdB]);
                    }
                    vecQb = Qb * vecQb;

                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);
                        B[nodeId] += vecQb[nodeIdB];
                    }

                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
            }
            // Φ1
            for (int portId = 0; portId < portCnt; portId++)
            {
                var Qb = Qbs[portId];
                var Rb = Rbs[portId];
                var Tb = Tbs[portId];
                int nodeCntB = Qb.RowLength;
                int abcOrderFor1 = ABCOrdersFor1[portId];
                int abcOrderFor2 = ABCOrdersFor2[portId];

                if (abcOrderFor1 == 1 && abcOrderFor2 == 1)
                {
                    // Traveling & Evanescent
                    double velo0 = VelosFor2[portId];
                    double alpha0 = AlphasFor1[portId];
                    double[] workEzPrev = new double[nodeCntB];
                    double[] workEzPrev2 = new double[nodeCntB];
                    double[] workPzPrev = new double[nodeCntB];
                    double[] workPzPrev2 = new double[nodeCntB];
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        // Ez
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);
                        workEzPrev[nodeIdB] = EzPzPrev[nodeId];
                        workEzPrev2[nodeIdB] = EzPzPrev2[nodeId];
                    }
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        // Pz
                        int nodeId = nodeIdB + Pz1Offsets[portId];
                        workPzPrev[nodeIdB] = EzPzPrev[nodeId];
                        workPzPrev2[nodeIdB] = EzPzPrev2[nodeId];
                    }

                    double[] vecQb = new double[nodeCntB];
                    double[] vecRb = new double[nodeCntB];
                    double[] vecTb = new double[nodeCntB];
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        vecQb[nodeIdB] =
                            (alpha0 / velo0) * (1.0 / (2.0 * dt)) * workEzPrev2[nodeIdB] +
                            -1.0 * alpha0 * ((1.0 - 2.0 * NewmarkBeta) * workPzPrev[nodeIdB] +
                            NewmarkBeta * workPzPrev2[nodeIdB]) +
                            (1.0 / (2.0 * dt * velo0)) * workPzPrev2[nodeIdB];
                        vecRb[nodeIdB] =
                            -1.0 * ((1.0 - NewmarkBeta) * workEzPrev[nodeIdB] +
                            NewmarkBeta * workEzPrev2[nodeIdB]);
                        vecTb[nodeIdB] =
                            (-1.0 / (dt * dt * velo0 * velo0)) *
                            (-2.0 * workEzPrev[nodeIdB] + workEzPrev2[nodeIdB]);
                    }
                    vecQb = Qb * vecQb;
                    vecRb = Rb * vecRb;
                    vecTb = Tb * vecTb;

                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int nodeId = nodeIdB + Pz1Offsets[portId];
                        B[nodeId] += vecQb[nodeIdB] + vecRb[nodeIdB] + vecTb[nodeIdB];
                    }
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
                double[] srcProfile = null;
                if (IsPointExcitation)
                {
                    // 中央に点波源を置く
                    srcProfile = new double[nodeCntB];
                    // Note: 節点は境界に沿って順に並んでいる(FEWorldでそのようになっている)
                    srcProfile[nodeCntB / 2] = 1.0;

                }
                else
                {
                    srcProfile = SrcProfiles[portId];
                }

                double srcU0 = 0.0;
                double srcU1 = 0.0;
                double srcU2 = 0.0;

                int n = TimeIndex;
                // 周波数
                double srcFreq = SrcFrequency;
                // 角周波数
                double srcOmega = 2.0 * Math.PI * srcFreq;

                if (IsGaussian)
                {
                    if (GaussianType == GaussianType.Normal)
                    {
                        // ガウシアンパルス
                        // Note: veloはC0
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
                        // Note: veloはvpx
                        if ((n * dt) <= (2.0 * GaussianT0 + Constants.PrecisionLowerLimit))
                        {
                            srcU0 = Math.Sin(srcOmega * ((n) * dt - GaussianT0)) *
                                Math.Exp(-1.0 * ((n) * dt - GaussianT0) * ((n) * dt - GaussianT0) /
                                (2.0 * GaussianTp * GaussianTp));
                            srcU1 = Math.Sin(srcOmega * ((n - 1) * dt - GaussianT0)) *
                                Math.Exp(-1.0 * ((n - 1) * dt - GaussianT0) * ((n - 1) * dt - GaussianT0) /
                                (2.0 * GaussianTp * GaussianTp));
                            srcU2 = Math.Sin(srcOmega * ((n + 1) * dt - GaussianT0)) *
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
                    srcU0 = Math.Sin(srcOmega * (n) * dt);
                    srcU1 = Math.Sin(srcOmega * (n - 1) * dt);
                    srcU2 = Math.Sin(srcOmega * (n + 1) * dt);
                }

                {
                    double[] srcUt = new double[nodeCntB];
                    double vpx = srcOmega / srcBetaX;
                    System.Diagnostics.Debug.Assert(srcProfile.Length == nodeCntB);
                    for (int i = 0; i < nodeCntB; i++)
                    {
                        srcUt[i] = srcProfile[i] * (srcU2 - srcU1) / (2.0 * dt);
                    }
                    double[] work = IvyFEM.Lapack.Functions.dscal(srcUt, (-2.0 / vpx));
                    double[] vecQb = Qb * work;
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);
                        B[nodeId] += vecQb[nodeIdB];
                    }
                }
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
                    int portId = portCnt + refIndex;
                    // モード
                    IvyFEM.Lapack.DoubleMatrix ryy1D;
                    IvyFEM.Lapack.DoubleMatrix txx1D;
                    IvyFEM.Lapack.DoubleMatrix uzz1D;
                    System.Numerics.Complex[] betas;
                    System.Numerics.Complex[][] eVecs;
                    double alpha;
                    EMWaveguide1DEigenBaseFEM eigen1DBaseFEM;
                    CalcEigen(
                        portId, freq,
                        out ryy1D, out txx1D, out uzz1D, out betas, out eVecs, out alpha, out eigen1DBaseFEM);

                    // 振幅分布
                    System.Numerics.Complex[] freqEz = freqDomainDatass[freqIndex];
                    // モード振幅の算出
                    int iMode = 0; // 基本モード
                    System.Numerics.Complex b = eigen1DBaseFEM.CalcModeAmp(omega, iMode, betas, eVecs, freqEz);
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
