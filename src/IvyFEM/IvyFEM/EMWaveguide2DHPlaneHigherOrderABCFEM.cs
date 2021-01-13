using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class EMWaveguide2DHPlaneHigherOrderABCFEM : FEM
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
        /// 周波数
        /// </summary>
        public double Frequency { get; set; } = 0.0;
        /// <summary>
        /// ABCで使う位相定数 (default: リスト無し(defaultの位相定数), -1を指定してもdefaultの位相定数になる)
        /// </summary>
        public IList<double> BetasToSet { get; set; } = new List<double>();
        /// <summary>
        /// 1D固有値問題で減衰定数を用いる？
        /// </summary>
        public IList<bool> IsEigen1DUseDecayParameters { get; set; } = new List<bool>(); 
        /// <summary>
        /// 1D固有値問題のクラッド比誘電率
        /// </summary>
        public IList<double> Eigen1DCladdingEp { get; set; } = new List<double>();
        /// <summary>
        /// 減衰定数を持ってくる1D固有値問題のポート
        /// </summary>
        public IList<int> DecayParameterEigen1DPortIds { get; set; } = new List<int>();

        /// <summary>
        /// TMモード？
        /// </summary>
        public bool IsTMMode { get; set; } = false;
        /// <summary>
        /// 観測点ポート数
        /// </summary>
        public int RefPortCount { get; set; } = 0;

        /// <summary>
        /// 節点数(ABCを除く)
        /// </summary>
        private int NodeCount = 0;
        /// <summary>
        /// 節点数
        /// </summary>
        private int NodeCountPlusABC = 0;
        /// <summary>
        /// [A]
        /// </summary>
        private IvyFEM.Linear.ComplexSparseMatrix A = null;
        /// <summary>
        /// {b}
        /// </summary>
        private System.Numerics.Complex[] B = null;
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
        public IList<System.Numerics.Complex> SrcBetaXs { get; private set; } = null;
        /// <summary>
        /// 境界の界のモード分布(ポート単位)
        /// </summary>
        public IList<System.Numerics.Complex[]> SrcFVecs { get; private set; } = null;
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
        private IList<double[]> AlphasFor1 = null;
        /// <summary>
        /// Evanescent Wave b0_abc
        /// </summary>
        private IList<double> B0AbcsFor1 = null;
        /// <summary>
        /// Evanescent Wave b_abc
        /// </summary>
        private IList<double[]> BAbcsFor1 = null;

        /// <summary>
        /// 吸収境界波の位相定数リスト(ポート単位)
        /// </summary>
        private IList<System.Numerics.Complex[]> BetasFor2 = null;
        /// <summary>
        /// 吸収境界 b0_abc
        /// </summary>
        private IList<System.Numerics.Complex> B0AbcsFor2 = null;
        /// <summary>
        /// 吸収境界 b_abc
        /// </summary>
        private IList<System.Numerics.Complex[]> BAbcsFor2 = null;

        /// <summary>
        /// Φ1(1)のオフセット(Evanescent)
        /// </summary>
        private IList<int> Pz1Offsets1 = null;
        /// <summary>
        /// Φ1(2)のオフセット(Traveling)
        /// </summary>
        private IList<int> Pz1Offsets2 = null;

        /// <summary>
        ///  電界※ABC分の自由度を除いたもの
        /// </summary>
        public System.Numerics.Complex[] Ez { get; private set; } = null;
        /// <summary>
        /// 電界
        /// </summary>
        private System.Numerics.Complex[] EzPz = null;
        /// <summary>
        /// Sパラメータ
        /// </summary>
        public System.Numerics.Complex[][] S { get; private set; }
        /// <summary>
        /// 固有値問題 EMWaveguide1DEigenFEM or EMWaveguide1DOepnEigenFEM
        /// </summary>
        public EMWaveguide1DEigenBaseFEM[] EigenBaseFEMs { get; private set; }

        public EMWaveguide2DHPlaneHigherOrderABCFEM(FEWorld world)
        {
            World = world;
        }

        public override void Solve()
        {
            // 周波数
            double freq = Frequency;
            // 角周波数
            double omega = 2.0 * Math.PI * freq;
            // 波長
            double waveLength = Constants.C0 / freq;
            // 波数
            double k0 = 2.0 * Math.PI / waveLength;

            //--------------------------------------------------------------
            // 全体行列
            //--------------------------------------------------------------
            CalcA(freq);

            int nodeCnt = NodeCount;
            int nodeCntPlusABC = NodeCountPlusABC;

            //--------------------------------------------------------------
            // 残差
            //--------------------------------------------------------------
            CalcB();

            //------------------------------------------------------------------
            // Ezを求める
            //------------------------------------------------------------------
            System.Numerics.Complex[] X;
            Solver.ComplexSolve(out X, A, B);
            EzPz = X;

            // ABC分を除いた電界を取得
            Ez = new System.Numerics.Complex[nodeCnt];
            for (int nodeId = 0; nodeId < nodeCnt; nodeId++)
            {
                Ez[nodeId] = EzPz[nodeId];
            }

            //------------------------------------------------------------------
            // Sマトリクスを求める
            //------------------------------------------------------------------
            S = CalcS(omega);
        }

        private void CalcA(double freq)
        {
            A = null;
            Qbs = new List<IvyFEM.Lapack.DoubleMatrix>();
            Rbs = new List<IvyFEM.Lapack.DoubleMatrix>();
            Tbs = new List<IvyFEM.Lapack.DoubleMatrix>();
            SrcBetaXs = new List<System.Numerics.Complex>();
            SrcFVecs = new List<System.Numerics.Complex[]>();
            SrcDecayParameters = new List<double>();

            // 角周波数
            double omega = 2.0 * Math.PI * freq;
            // 波長
            double waveLength = Constants.C0 / freq;
            // 波数
            double k0 = 2.0 * Math.PI / waveLength;

            /////////////////////////////////////////////////////////
            int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1; // 参照面と励振源を除く
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
                System.Diagnostics.Debug.Assert(abcOrderFor1 >= 1 || abcOrderFor2 >= 1);
                ABCOrdersFor1.Add(abcOrderFor1);
                ABCOrdersFor2.Add(abcOrderFor2);
            }
            //------------------------------------------------------
            // 節点数
            //-----------------------------------------------------
            int nodeCnt = (int)World.GetNodeCount(QuantityId);
            int nodeCntPlusABC = nodeCnt;
            for (int portId = 0; portId < portCnt; portId++)
            {
                int nodeCntB = (int)World.GetPortNodeCount(QuantityId, (uint)portId);
                int abcOrderFor1 = ABCOrdersFor1[portId];
                int abcOrderFor2 = ABCOrdersFor2[portId];
                if (abcOrderFor1 >= 1)
                {
                    nodeCntPlusABC += nodeCntB * (abcOrderFor1 - 1);
                }
                if (abcOrderFor1 >= 1 && abcOrderFor2 >= 1)
                {
                    nodeCntPlusABC += nodeCntB;
                }
                if (abcOrderFor2 >= 1)
                {
                    nodeCntPlusABC += nodeCntB * (abcOrderFor2 - 1);
                }
            }
            NodeCount = nodeCnt;
            NodeCountPlusABC = nodeCntPlusABC;

            A = new IvyFEM.Linear.ComplexSparseMatrix(nodeCntPlusABC, nodeCntPlusABC);

            //------------------------------------------------------
            // 剛性行列、質量行列を作成
            //------------------------------------------------------
            CalcKM(k0);

            //------------------------------------------------------
            // モード分布計算
            //------------------------------------------------------
            System.Diagnostics.Debug.Assert(
                IsEigen1DUseDecayParameters.Count == 0 ||
                IsEigen1DUseDecayParameters.Count == (portCnt + RefPortCount + 1));
            EigenBaseFEMs = new EMWaveguide1DEigenBaseFEM[(portCnt + RefPortCount + 1)];
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
                    portId, freq,
                    out ryy1D, out txx1D, out uzz1D, out betas, out eVecs, out alpha,
                    out eigen1DBaseFEM);

                EigenBaseFEMs[portId] = eigen1DBaseFEM;
                int nodeCntB = ryy1D.RowLength;
                Qbs.Add(ryy1D);
                Rbs.Add(txx1D);
                Tbs.Add(uzz1D);

                System.Diagnostics.Debug.WriteLine("port = {0} mode Count = {1}", portId, betas.Length);
                // 基本モード
                int iMode = 0;
                System.Diagnostics.Debug.Assert(World.GetIncidentModeId(QuantityId) == 0);//現状基本モード固定
                System.Numerics.Complex beta = betas[iMode];
                System.Numerics.Complex[] fVec = eVecs[iMode];
                SrcBetaXs.Add(beta);
                SrcFVecs.Add(fVec);
                SrcDecayParameters.Add(alpha);
            }

            ////////////////////////////////////////////////////////////////////////////////////////
            // 吸収境界パラメータ
            AlphasFor1 = new List<double[]>();
            B0AbcsFor1 = new List<double>();
            BAbcsFor1 = new List<double[]>();
            for (int portId = 0; portId < portCnt; portId++)
            {
                int abcOrderFor1 = ABCOrdersFor1[portId];
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
                double[] alpha =  new double[abcOrderFor1];
                for (int order = 0; order < abcOrderFor1; order++)
                {
                    alpha[order] = alpha0;
                }
                AlphasFor1.Add(alpha);
            }
            for (int portId = 0; portId < portCnt; portId++)
            {
                int abcOrderFor1 = ABCOrdersFor1[portId];
                double[] alpha = AlphasFor1[portId];
                double B0Abc = alpha.Length > 0 ? alpha[0] : 0.0;
                double[] BAbc = new double[(abcOrderFor1 >= 1 ? (abcOrderFor1 - 1) : 0)];
                for (int order = 0; order < abcOrderFor1 - 1; order++)
                {
                    BAbc[order] = alpha[order] + alpha[order + 1];
                }
                B0AbcsFor1.Add(B0Abc);
                BAbcsFor1.Add(BAbc);
            }

            BetasFor2 = new List<System.Numerics.Complex[]>();
            B0AbcsFor2 = new List<System.Numerics.Complex>();
            BAbcsFor2 = new List<System.Numerics.Complex[]>();

            System.Diagnostics.Debug.Assert(BetasToSet.Count == 0 || BetasToSet.Count == portCnt);
            for (int portId = 0; portId < portCnt; portId++)
            {
                int abcOrderFor2 = ABCOrdersFor2[portId];
                System.Numerics.Complex[] beta = new System.Numerics.Complex[abcOrderFor2];
                for (int order = 0; order < abcOrderFor2; order++)
                {
                    System.Numerics.Complex beta0 = -1.0;
                    if (BetasToSet.Count == portCnt)
                    {
                        beta0 = BetasToSet[portId];
                    }
                    if (beta0.Real < 0) // beta0: -1
                    {
                        //beta0 = omega / Constants.C0;
                        //beta0 = k0;
                        System.Numerics.Complex srcBetaX = SrcBetaXs[portId];
                        beta0 = srcBetaX;
                    }
                    beta[order] = beta0;
                }
                BetasFor2.Add(beta);
            }
            for (int portId = 0; portId < portCnt; portId++)
            {
                int abcOrderFor2 = ABCOrdersFor2[portId];
                System.Numerics.Complex[] beta = BetasFor2[portId];
                System.Numerics.Complex B0Abc = beta.Length > 0 ?
                    System.Numerics.Complex.ImaginaryOne * beta[0] : 0.0;
                System.Numerics.Complex[] BAbc = 
                    new System.Numerics.Complex[(abcOrderFor2 >= 1 ? (abcOrderFor2 - 1) : 0)];
                for (int order = 0; order < (abcOrderFor2 - 1); order++)
                {
                    BAbc[order] = System.Numerics.Complex.ImaginaryOne * (beta[order] + beta[order + 1]);
                }
                B0AbcsFor2.Add(B0Abc);
                BAbcsFor2.Add(BAbc);
            }

            ////////////////////////////////////////////////////////////////////////////////////////
            // Φ1の開始位置
            Pz1Offsets1 = new List<int>();
            Pz1Offsets2 = new List<int>();
            for (int portId = 0; portId < portCnt; portId++)
            {
                int offset1 = 0;
                if (portId == 0)
                {
                    offset1 = nodeCnt;
                }
                else
                {
                    int nodeCntBPrevPort = (int)World.GetPortNodeCount(QuantityId, (uint)(portId - 1));
                    int abcOrderFor1PrevPort = ABCOrdersFor1[portId - 1];
                    int abcOrderFor2PrevPort = ABCOrdersFor2[portId - 1];

                    offset1 = Pz1Offsets1[portId - 1];
                    if (abcOrderFor1PrevPort >= 1)
                    {
                        offset1 += nodeCntBPrevPort * (abcOrderFor1PrevPort - 1);
                    }
                    if (abcOrderFor2PrevPort >= 1 && abcOrderFor1PrevPort >= 1)
                    {
                        offset1 += nodeCntBPrevPort;
                    }
                    if (abcOrderFor2PrevPort >= 1)
                    {
                        offset1 += nodeCntBPrevPort * (abcOrderFor2PrevPort - 1);
                    }
                }
                int offset2 = offset1;
                {
                    int nodeCntB = (int)World.GetPortNodeCount(QuantityId, (uint)portId);
                    int abcOrderFor1 = ABCOrdersFor1[portId];
                    int abcOrderFor2 = ABCOrdersFor2[portId];
                    if (abcOrderFor1 >= 1)
                    {
                        offset2 += nodeCntB * (abcOrderFor1 - 1);
                    }
                    if (abcOrderFor2 >= 1 && abcOrderFor1 >= 1)
                    {
                        offset2 += nodeCntB;
                    }
                }
                Pz1Offsets1.Add(offset1);
                Pz1Offsets2.Add(offset2);
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

                        if (abcOrderFor1 >= 1)
                        {
                            // Evanescent
                            double B0AbcFor1 = B0AbcsFor1[portId];
                            A[rowNodeId, colNodeId] +=
                                B0AbcFor1  * Qb[rowNodeIdB, colNodeIdB];
                        }
                        else if (abcOrderFor2 >= 1 && abcOrderFor1 == 0)
                        {
                            // Traveling
                            System.Numerics.Complex B0AbcFor2 = B0AbcsFor2[portId];
                            A[rowNodeId, colNodeId] +=
                                B0AbcFor2 * Qb[rowNodeIdB, colNodeIdB];
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

                    for (int colNodeIdB = 0; colNodeIdB < nodeCntB; colNodeIdB++)
                    {
                        if (abcOrderFor1 >= 1)
                        {
                            // Evanescent
                            if (abcOrderFor1 >= 2)
                            {
                                // Φ1(1)
                                int colNodeId1 = colNodeIdB + Pz1Offsets1[portId];
                                A[rowNodeId, colNodeId1] += -1.0 * Qb[rowNodeIdB, colNodeIdB];
                            }
                        }
                        else if (abcOrderFor2 >= 1 && abcOrderFor1 == 0)
                        {
                            // Traveling
                            if (abcOrderFor2 >= 2)
                            {
                                System.Diagnostics.Debug.Assert(Pz1Offsets1[portId] == Pz1Offsets2[portId]);
                                // Φ1(2)
                                int colNodeId1 = colNodeIdB + Pz1Offsets2[portId];
                                A[rowNodeId, colNodeId1] += -1.0 * Qb[rowNodeIdB, colNodeIdB];
                            }
                        }
                        else
                        {
                            System.Diagnostics.Debug.Assert(false);
                        }
                    }
                }
            }

            ////////////////////////////////////////////////////////////////////////////////////////
            // Evanescent Wave
            // Φ1(1) ～ Φ(ABC_order -1)(1)
            for (int portId = 0; portId < portCnt; portId++)
            {
                var Qb = Qbs[portId];
                var Rb = Rbs[portId];
                var Tb = Tbs[portId];
                int abcOrderFor1 = ABCOrdersFor1[portId];
                int abcOrderFor2 = ABCOrdersFor2[portId];
                double[] alphaFor1 = AlphasFor1[portId];
                double[] BAbcFor1 = BAbcsFor1[portId];
                int nodeCntB = Qb.RowLength;

                // Φ(order)に関する式
                for (int order = 0; order < (abcOrderFor1 - 1); order++)
                {
                    for (int rowNodeIdB = 0; rowNodeIdB < nodeCntB; rowNodeIdB++)
                    {
                        // order
                        int rowNodeId1 = rowNodeIdB + nodeCntB * order + Pz1Offsets1[portId];

                        for (int colNodeIdB = 0; colNodeIdB < nodeCntB; colNodeIdB++)
                        {
                            int colNodeId1 = 0;
                            // Φ(order)
                            colNodeId1 = colNodeIdB + nodeCntB * order + Pz1Offsets1[portId];
                            A[rowNodeId1, colNodeId1] += BAbcFor1[order] * Qb[rowNodeIdB, colNodeIdB];
                            // Φ(order-1)
                            if (order == 0)
                            {
                                // Φ0 (Ez)
                                int colCoId = World.PortNode2Coord(QuantityId, (uint)portId, colNodeIdB);
                                colNodeId1 = World.Coord2Node(QuantityId, colCoId);
                            }
                            else
                            {
                                colNodeId1 = colNodeIdB + nodeCntB * (order - 1) + Pz1Offsets1[portId];
                            }
                            A[rowNodeId1, colNodeId1] +=
                                -alphaFor1[order] * alphaFor1[order] * Qb[rowNodeIdB, colNodeIdB] +
                                -k0 * k0 * Tb[rowNodeIdB, colNodeIdB] +
                               Rb[rowNodeIdB, colNodeIdB];
                            // Φ(order + 1)
                            colNodeId1 = 0;
                            if (order == (abcOrderFor1 - 2) && abcOrderFor2 == 0)
                            {
                                // なし
                            }
                            else
                            {
                                colNodeId1 = colNodeIdB + nodeCntB * (order + 1) + Pz1Offsets1[portId];
                                A[rowNodeId1, colNodeId1] +=
                                    -1.0 * Qb[rowNodeIdB, colNodeIdB];
                            }
                        }
                    }
                }
            }

            ///////////////////////////////////////////////////////////////////////////
            // Φ(ABCOrder - 1)(1) = Φ0(2)に関する式
            for (int portId = 0; portId < portCnt; portId++)
            {
                var Qb = Qbs[portId];
                var Rb = Rbs[portId];
                var Tb = Tbs[portId];
                double[] alpha = AlphasFor1[portId];
                System.Numerics.Complex B0ABCFor2 = B0AbcsFor2[portId];
                System.Numerics.Complex[] betaFor2 = BetasFor2[portId];
                int nodeCntB = Qb.RowLength;
                int abcOrderFor1 = ABCOrdersFor1[portId];
                int abcOrderFor2 = ABCOrdersFor2[portId];

                // Φ(ABCOrder - 1)(1) = Φ0(2)
                if (abcOrderFor1 >= 1 && abcOrderFor2 >= 1)
                {
                    int order = abcOrderFor1 - 1;
                    for (int rowNodeIdB = 0; rowNodeIdB < nodeCntB; rowNodeIdB++)
                    {
                        // order
                        int rowNodeId1 = rowNodeIdB + nodeCntB * order + Pz1Offsets1[portId];

                        for (int colNodeIdB = 0; colNodeIdB < nodeCntB; colNodeIdB++)
                        {
                            int colNodeId1 = 0;
                            // Φ(order)
                            colNodeId1 = colNodeIdB + nodeCntB * order + Pz1Offsets1[portId];
                            A[rowNodeId1, colNodeId1] +=
                                B0ABCFor2 * Qb[rowNodeIdB, colNodeIdB] +
                                alpha[order] * Qb[rowNodeIdB, colNodeIdB];

                            // Φ(order-1)
                            colNodeId1 = 0;
                            if (order == 0)
                            {
                                // Φ0 (Ez)
                                int colCoId = World.PortNode2Coord(QuantityId, (uint)portId, colNodeIdB);
                                colNodeId1 = World.Coord2Node(QuantityId, colCoId);
                            }
                            else
                            {
                                colNodeId1 = colNodeIdB + nodeCntB * (order - 1) + Pz1Offsets1[portId];
                            }
                            A[rowNodeId1, colNodeId1] +=
                                -alpha[order] * alpha[order] * Qb[rowNodeIdB, colNodeIdB] +
                                -k0 * k0 * Tb[rowNodeIdB, colNodeIdB] +
                               Rb[rowNodeIdB, colNodeIdB];

                            // Φ(order + 1)
                            if (abcOrderFor2 <= 1)
                            {

                            }
                            else
                            {
                                colNodeId1 = colNodeIdB + nodeCntB * (order + 1) + Pz1Offsets1[portId];
                                A[rowNodeId1, colNodeId1] +=
                                    -1.0 * Qb[rowNodeIdB, colNodeIdB];
                            }
                        }
                    }
                }
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////
            // Traveling Wave
            // Φ1(2) ～ Φ(ABC_order - 1)(2)
            for (int portId = 0; portId < portCnt; portId++)
            {
                var Qb = Qbs[portId];
                var Rb = Rbs[portId];
                var Tb = Tbs[portId];
                System.Numerics.Complex[] betaFor2 = BetasFor2[portId];
                System.Numerics.Complex[] BAbcFor2 = BAbcsFor2[portId];
                int nodeCntB = Qb.RowLength;
                int abcOrderFor1 = ABCOrdersFor1[portId];
                int abcOrderFor2 = ABCOrdersFor2[portId];

                // Φ(order)に関する式
                for (int order = 0; order < (abcOrderFor2 - 1); order++)
                {
                    for (int rowNodeIdB = 0; rowNodeIdB < nodeCntB; rowNodeIdB++)
                    {
                        // order
                        int rowNodeId1 = rowNodeIdB + nodeCntB * order + Pz1Offsets2[portId];

                        for (int colNodeIdB = 0; colNodeIdB < nodeCntB; colNodeIdB++)
                        {
                            int colNodeId1 = 0;
                            // Φ(order)
                            colNodeId1 = colNodeIdB + nodeCntB * order + Pz1Offsets2[portId];
                            A[rowNodeId1, colNodeId1] +=
                                BAbcFor2[order] * Qb[rowNodeIdB, colNodeIdB];

                            // Φ(order-1)
                            colNodeId1 = 0;
                            if (order == 0 && abcOrderFor1 <= 1)
                            {
                                // Φ0 (Ez)
                                int colCoId = World.PortNode2Coord(QuantityId, (uint)portId, colNodeIdB);
                                colNodeId1 = World.Coord2Node(QuantityId, colCoId);
                            }
                            else
                            {
                                colNodeId1 = colNodeIdB + nodeCntB * (order - 1) + Pz1Offsets2[portId];
                            }
                            A[rowNodeId1, colNodeId1] +=
                                betaFor2[order] * betaFor2[order] * Qb[rowNodeIdB, colNodeIdB] +
                                -k0 * k0 * Tb[rowNodeIdB, colNodeIdB] +
                                Rb[rowNodeIdB, colNodeIdB];

                            // Φ(order + 1)
                            colNodeId1 = 0;
                            if (order == (abcOrderFor2 - 2))
                            {

                                // なし
                            }
                            else
                            {
                                colNodeId1 = colNodeIdB + nodeCntB * (order + 1) + Pz1Offsets2[portId];
                                A[rowNodeId1, colNodeId1] +=
                                    -1.0 * Qb[rowNodeIdB, colNodeIdB];
                            }
                        }
                    }
                }
            }

            // バンド幅縮小を座標を用いて行う
            if (Solver is IvyFEM.Linear.LapackEquationSolver)
            {
                var solver = Solver as IvyFEM.Linear.LapackEquationSolver;
                if (solver.IsOrderingToBandMatrix)
                {
                    int[] coIds = new int[nodeCntPlusABC];
                    for (int nodeId = 0; nodeId < nodeCntPlusABC; nodeId++)
                    {
                        if (nodeId < nodeCnt)
                        {
                            int coId = World.Node2Coord(QuantityId, nodeId);
                            coIds[nodeId] = coId;
                        }
                        else
                        {
                            coIds[nodeId] = -1;
                        }
                    }
                    for (int portId = 0; portId < portCnt; portId++)
                    {
                        var Qb = Qbs[portId];
                        int nodeCntB = Qb.RowLength;
                        int abcOrderFor1 = ABCOrdersFor1[portId];
                        int abcOrderFor2 = ABCOrdersFor2[portId];
                        if (abcOrderFor1 >= 2)
                        {
                            for (int order = 0; order < (abcOrderFor1 - 1); order++)
                            {
                                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                                {
                                    int nodeId = nodeIdB + nodeCntB * order + Pz1Offsets1[portId];
                                    int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                                    coIds[nodeId] = coId;
                                }
                            }
                        }
                        if (abcOrderFor1 >= 1 && abcOrderFor2 >= 1)
                        {
                            int order = abcOrderFor1 - 1;
                            for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                            {
                                int nodeId = nodeIdB + nodeCntB * order + Pz1Offsets1[portId];
                                int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                                coIds[nodeId] = coId;
                            }
                        }
                        if (abcOrderFor2 >= 2)
                        {
                            for (int order = 0; order < (abcOrderFor2 - 1); order++)
                            {
                                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                                {
                                    int nodeId = nodeIdB + nodeCntB * order + Pz1Offsets2[portId];
                                    int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                                    coIds[nodeId] = coId;
                                }
                            }
                        }
                    }
                    double[][] coords = new double[nodeCntPlusABC][];
                    for (int nodeId = 0; nodeId < nodeCntPlusABC; nodeId++)
                    {
                        int coId = coIds[nodeId];
                        double[] coord = World.GetCoord(QuantityId, coId);
                        coords[nodeId] = coord;
                    }

                    solver.CoordsForBandMatrix = coords;
                }
            }
        }

        private void CalcKM(double k0)
        {
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

                        double a = 
                            maPxx * sNyNy[row, col] + maPyy * sNxNx[row, col] -
                            k0 * k0 * maQzz * sNN[row, col];
                        A[rowNodeId, colNodeId] += a;
                    }
                }
            }
        }

        private void CalcEigen(
            int portId, double freq,
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
                eigen1DFEM.CladdingEp = Eigen1DCladdingEp[portId];
                eigen1DFEM.IsTMMode = IsTMMode;
                eigen1DFEM.Frequency = freq;
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
                eigen1DFEM.IsTMMode = IsTMMode;
                eigen1DFEM.Frequency = freq;
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
            int nodeCntPlusABC = NodeCountPlusABC;
            int nodeCnt = NodeCount;
            B = new System.Numerics.Complex[nodeCntPlusABC];

            int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1;  // 参照面、励振源分引く

            //--------------------------------------------------------------
            // 励振源
            //--------------------------------------------------------------
            {
                int portId = portCnt + RefPortCount; // ポートリストの最後の要素が励振境界
                var Qb = Qbs[portId];
                int nodeCntB = Qb.RowLength;
                System.Numerics.Complex srcBetaX = SrcBetaXs[portId];
                System.Numerics.Complex[] srcEVec = SrcFVecs[portId];

                {
                    System.Numerics.Complex betaX = srcBetaX;
                    System.Diagnostics.Debug.Assert(srcEVec.Length == nodeCntB);
                    System.Diagnostics.Debug.Assert(betaX.Real >= 0);
                    System.Numerics.Complex[] work = new System.Numerics.Complex[nodeCntB];
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        work[nodeIdB] = (2.0 * System.Numerics.Complex.ImaginaryOne * betaX) * srcEVec[nodeIdB];
                    }
                    var QbZ = new IvyFEM.Lapack.ComplexMatrix(Qb);
                    System.Numerics.Complex[] vecQb = QbZ * work;
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);
                        B[nodeId] += vecQb[nodeIdB];
                    }
                }
            }
        }

        private System.Numerics.Complex[][] CalcS(double omega)
        {
            System.Diagnostics.Debug.Assert(RefPortCount > 0);
            int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1; // 参照面と励振源を除く
            int incidentPortId = World.GetIncidentPortId(QuantityId);
            int excitationPortId = portCnt + RefPortCount;

            System.Diagnostics.Debug.Assert(World.Mesh is Mesher2D);
            Mesher2D mesh = World.Mesh as Mesher2D;

            // 励振面から入射参照面までの距離と位相差の計算
            var portConditions = World.GetPortConditions(QuantityId);
            PortCondition[] tagtPortConditions = { portConditions[excitationPortId], portConditions[incidentPortId] };
            IList<OpenTK.Vector2d[]> portSEPts = new List<OpenTK.Vector2d[]>();
            foreach (PortCondition portCondition in tagtPortConditions)
            {
                OpenTK.Vector2d[] sePt = new OpenTK.Vector2d[2]; 
                IList<uint> eIds = portCondition.EIds;
                uint eId1 = eIds[0];
                Edge2D e1 = mesh.Cad.GetEdge(eId1);
                sePt[0] = e1.GetVertexCoord(true);
                uint eId2 = eIds[eIds.Count - 1];
                Edge2D e2 = mesh.Cad.GetEdge(eId2);
                sePt[1] = e2.GetVertexCoord(false);
                portSEPts.Add(sePt);
            }
            System.Numerics.Complex a;
            {
                // 励振面と入射面の距離を算出
                // (両面は平行であるものとする)
                OpenTK.Vector2d v1 = portSEPts[1][0];
                OpenTK.Vector2d v2 = portSEPts[0][0];
                OpenTK.Vector2d v3 = portSEPts[0][1];
                double distanceX = Math.Abs(IvyFEM.CadUtils2D.TriHeight(v1, v2, v3));

                System.Numerics.Complex betaX = SrcBetaXs[excitationPortId];
                // 入射面（ポート1)における振幅を計算
                a = 1.0 * System.Numerics.Complex.Exp(-1.0 * System.Numerics.Complex.ImaginaryOne * betaX * distanceX);
            }

            // Sマトリクスの計算
            var S = new System.Numerics.Complex[RefPortCount][];
            
            for (int refIndex = 0; refIndex < RefPortCount; refIndex++)
            {
                int portId = portCnt + refIndex;
                int nodeCntB = (int)World.GetPortNodeCount(QuantityId, (uint)portId);
                System.Numerics.Complex[] portEz = GetPortEz((uint)portId, Ez);
                int incidentModeId = -1;
                if (incidentPortId == portId)
                {
                    incidentModeId = (int)World.GetIncidentModeId(QuantityId);
                    // 現状0固定
                    System.Diagnostics.Debug.Assert(incidentModeId == 0);
                }

                EMWaveguide1DEigenBaseFEM eigenBaseFEM = EigenBaseFEMs[portId];
                System.Numerics.Complex[] betas = eigenBaseFEM.Betas;
                System.Numerics.Complex[][] ezEVecs = eigenBaseFEM.EzEVecs;
                int modeCnt = betas.Length;
                System.Numerics.Complex[] S1 = new System.Numerics.Complex[modeCnt];
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    System.Numerics.Complex b = eigenBaseFEM.CalcModeAmp(omega, iMode, betas, ezEVecs, portEz);
                    if (incidentModeId == iMode)
                    {
                        b += -a;
                    }
                    S1[iMode] = b / a;
                }
                S[refIndex] = S1;
            }
            return S;
        }

        private System.Numerics.Complex[] GetPortEz(uint portId, System.Numerics.Complex[] Ez)
        {
            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, portId);
            System.Numerics.Complex[] portEz = new System.Numerics.Complex[nodeCnt];
            for (int row = 0; row < nodeCnt; row++)
            {
                int coId = World.PortNode2Coord(QuantityId, portId, row);
                int nodeId = World.Coord2Node(QuantityId, coId);
                portEz[row] = Ez[nodeId];
            }
            return portEz;
        }
    }
}
