using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class EMWaveguide3DFirstOrderABCTDFEM : FEM
    {
        //---------------------------------
        // ※ポートの設定順序
        // ABC境界、参照(観測)面、励振源の順
        //---------------------------------
        // 磁界？
        public bool IsMagneticField { get; set; } = false;

        public uint QuantityId { get; private set; } = 0;

        private double NewmarkBeta = 1.0 / 4.0;

        public double TimeStep { get; set; } = 0.0;

        /// <summary>
        /// 時刻ステップ数
        /// </summary>
        public int TimeLoopCnt { get; set; } = 0;
        /// <summary>
        /// 計算時刻インデックス
        /// </summary>
        public int TimeIndex { get; set; } = 0;

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
        /// <summary>
        /// ガウシアンパルスの振幅
        /// </summary>
        public double GaussianAmp { get; set; } = 1.0;

        /// <summary>
        /// 観測点ポート数
        /// </summary>
        public int RefPortCount { get; set; } = 0;

        /// <summary>
        /// [A]
        /// </summary>
        private IvyFEM.Linear.DoubleSparseMatrix A = null;
        /// <summary>
        /// {b}
        /// </summary>
        private double[] B = null;

        /// <summary>
        ///  電界（現在値）
        /// </summary>
        public double[] E { get; private set; } = null;
        public double[] CoordExyz { get; private set; }

        /// <summary>
        /// <summary>
        /// 電界(1つ前)
        /// </summary>
        private double[] EPrev = null;
        /// <summary>
        /// 電界(2つ前)
        /// </summary>
        private double[] EPrev2 = null;

        public EMWaveguide2DEigenFEMForPort[] EigenFEMs { get; private set; }

        /// <summary>
        /// 観測点の電界(時間変化リスト)
        /// </summary>
        public IList<IList<double[]>> RefTimeEsss { get; private set; } = new List<IList<double[]>>();

        public EMWaveguide3DFirstOrderABCTDFEM(FEWorld world)
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
                if (RefPortCount > 0)
                {
                    cnt = RefPortCount; 
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
                RefTimeEsss.Clear();
                for (int i = 0; i < cnt; i++)
                {
                    var refTimeEss = new List<double[]>();
                    RefTimeEsss.Add(refTimeEss);
                }
            }

            if (TimeIndex == 0)
            {
                //--------------------------------------------------------------
                // 全体行列
                //--------------------------------------------------------------
                CalcA();
            }

            //--------------------------------------------------------------
            // 電界
            //--------------------------------------------------------------
            // 以下TimeIndexに対する計算
            double time = TimeIndex * TimeStep;
            double dt = TimeStep;
            EPrev.CopyTo(EPrev2, 0);
            E.CopyTo(EPrev, 0);
            E = null;

            //--------------------------------------------------------------
            // 残差
            //--------------------------------------------------------------
            CalcB();

            //------------------------------------------------------------------
            // Eを求める
            //------------------------------------------------------------------
            {
                double[] X;
                Solver.DoubleSolve(out X, A, B);
                E = X;
            }

            // 節点におけるEを求める
            double[] coordExyz;
            CalcCoordExyz(E, out coordExyz);
            CoordExyz = coordExyz;

            // 参照面の界の格納
            if (RefPortCount > 0)
            {
                int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1; // 参照面と励振源を除く
                for (int refIndex = 0; refIndex < RefPortCount; refIndex++)
                {
                    int portId = portCnt + refIndex;
                    int edgeNodeCntB = (int)World.GetPortEdgeNodeCount(QuantityId, (uint)portId);
                    double[] fValues = new double[edgeNodeCntB];
                    for (int edgeNodeIdB = 0; edgeNodeIdB < edgeNodeCntB; edgeNodeIdB++)
                    {
                        int edgeId = World.PortEdgeNode2Edge(QuantityId, (uint)portId, edgeNodeIdB);
                        int edgeNodeId = World.Edge2EdgeNode(QuantityId, edgeId);
                        fValues[edgeNodeIdB] = E[edgeNodeId];
                    }
                    var refTimeEss = RefTimeEsss[refIndex];
                    refTimeEss.Add(fValues);
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

            double dt = TimeStep;
            double beta = NewmarkBeta;

            int edgeNodeCnt = (int)World.GetEdgeNodeCount(QuantityId);

            //------------------------------------------------------
            // 電界
            //-----------------------------------------------------
            E = new double[edgeNodeCnt];
            EPrev = new double[edgeNodeCnt];
            EPrev2 = new double[edgeNodeCnt];

            //------------------------------------------------------
            // 剛性行列、質量行列を作成
            //------------------------------------------------------
            A = new IvyFEM.Linear.DoubleSparseMatrix(edgeNodeCnt, edgeNodeCnt);
            _CalcA();

            //------------------------------------------------------
            // モード分布計算
            //------------------------------------------------------
            // 周波数
            double srcFreq = SrcFrequency;
            // 角周波数
            double srcOmega = 2.0 * Math.PI * srcFreq;

            int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1; // 参照面と励振源を除く

            EigenFEMs = new EMWaveguide2DEigenFEMForPort[(portCnt + RefPortCount + 1)];
            for (int portId = 0; portId < (portCnt + RefPortCount + 1); portId++)
            {
                System.Numerics.Complex[] betas;
                System.Numerics.Complex[][] etEVecs;
                System.Numerics.Complex[][] ezEVecs;
                EMWaveguide2DEigenFEMForPort eigenFEM;
                CalcEigen(
                    portId, srcFreq,
                    out betas, out etEVecs, out ezEVecs, out eigenFEM);
                EigenFEMs[portId] = eigenFEM;

                System.Diagnostics.Debug.WriteLine("port = {0} mode Count = {1}", portId, betas.Length);
            }

            ////////////////////////////////////////////////////////////////////////////////////////
            // 吸収境界
            for (int portId = 0; portId < portCnt; portId++)
            {
                int edgeNodeCntB = (int)World.GetPortEdgeNodeCount(QuantityId, (uint)portId);
                var eigenFEM = EigenFEMs[portId];
                var Rtt = eigenFEM.Rtt;
                // 基本モード
                double srcBeta = eigenFEM.Betas[0].Real;
                double srcVelo = srcOmega / srcBeta;

                for (int rowEdgeNodeIdB = 0; rowEdgeNodeIdB < edgeNodeCntB; rowEdgeNodeIdB++)
                {
                    // E
                    int rowEdgeId = World.PortEdgeNode2Edge(QuantityId, (uint)portId, rowEdgeNodeIdB);
                    int rowEdgeNodeId = World.Edge2EdgeNode(QuantityId, rowEdgeId);

                    for (int colEdgeNodeIdB = 0; colEdgeNodeIdB < edgeNodeCntB; colEdgeNodeIdB++)
                    {
                        // E
                        int colEdgeId = World.PortEdgeNode2Edge(QuantityId, (uint)portId, colEdgeNodeIdB);
                        int colEdgeNodeId = World.Edge2EdgeNode(QuantityId, colEdgeId);

                        // Traveling
                        double velo = srcVelo;
                        //double velo = Constants.C0;
                        double matBValue = -1.0 * (1.0 / velo) * Rtt[rowEdgeNodeIdB, colEdgeNodeIdB];
                        A[rowEdgeNodeId, colEdgeNodeId] +=
                            (1.0 / (2.0 * dt)) * (-1.0 * matBValue);
                    }
                }
            }
        }

        private void CalcB()
        {
            double dt = TimeStep;
            double beta = NewmarkBeta;

            int edgeNodeCnt = (int)World.GetEdgeNodeCount(QuantityId);

            B = new double[edgeNodeCnt];
            _CalcB();

            int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1; // 参照面と励振源を除く

            // 周波数
            double srcFreq = SrcFrequency;
            // 角周波数
            double srcOmega = 2.0 * Math.PI * srcFreq;

            ////////////////////////////////////////////////////////////////////////////////////////
            // 吸収境界
            for (int portId = 0; portId < portCnt; portId++)
            {
                int edgeNodeCntB = (int)World.GetPortEdgeNodeCount(QuantityId, (uint)portId);
                var eigenFEM = EigenFEMs[portId];
                var Rtt = eigenFEM.Rtt;
                // 基本モード
                double srcBeta = eigenFEM.Betas[0].Real;
                double srcVelo = srcOmega / srcBeta;

                for (int rowEdgeNodeIdB = 0; rowEdgeNodeIdB < edgeNodeCntB; rowEdgeNodeIdB++)
                {
                    // E
                    int rowEdgeId = World.PortEdgeNode2Edge(QuantityId, (uint)portId, rowEdgeNodeIdB);
                    int rowEdgeNodeId = World.Edge2EdgeNode(QuantityId, rowEdgeId);

                    for (int colEdgeNodeIdB = 0; colEdgeNodeIdB < edgeNodeCntB; colEdgeNodeIdB++)
                    {
                        // E
                        int colEdgeId = World.PortEdgeNode2Edge(QuantityId, (uint)portId, colEdgeNodeIdB);
                        int colEdgeNodeId = World.Edge2EdgeNode(QuantityId, colEdgeId);

                        // Traveling
                        double velo = srcVelo;
                        //double velo = Constants.C0;
                        double matBValue = -1.0 * (1.0 / velo) * Rtt[rowEdgeNodeIdB, colEdgeNodeIdB];
                        B[rowEdgeNodeId] +=
                            (-1.0 * matBValue) * (1.0 / (2.0 * dt)) * EPrev2[colEdgeNodeId];
                    }
                }
            }

            //---------------------------------------------------
            // 励振源
            //---------------------------------------------------
            _CalcBSrc();
        }

        private void _CalcA()
        {
            double dt = TimeStep;
            double beta = NewmarkBeta;

            IList<uint> feIds = World.GetTetrahedronFEIds(QuantityId);
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = World.GetTetrahedronFE(QuantityId, feId);
                uint elemNodeCnt = tetFE.NodeCount;
                uint elemEdgeNodeCnt = tetFE.EdgeCount;
                int[] edgeIds = new int[elemEdgeNodeCnt];
                bool[] isReverses = new bool[elemEdgeNodeCnt];
                int[] edgeNodes = new int[elemEdgeNodeCnt];
                for (int iENode = 0; iENode < elemEdgeNodeCnt; iENode++)
                {
                    int[] coIds = tetFE.EdgeCoordIdss[iENode];
                    bool isReverse;
                    int edgeId = World.GetEdgeIdFromCoords(QuantityId, coIds[0], coIds[1], out isReverse);
                    int edgeNodeId = World.Edge2EdgeNode(QuantityId, edgeId);

                    edgeIds[iENode] = edgeId;
                    isReverses[iENode] = isReverse;
                    edgeNodes[iENode] = edgeNodeId;
                }

                Material ma0 = World.GetMaterial(tetFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is DielectricMaterial);
                var ma = ma0 as DielectricMaterial;
                // TODO: 複素誘電率対応する
                double maPxx = 0.0;
                double maPyy = 0.0;
                double maPzz = 0.0;
                double maQxx = 0.0;
                double maQyy = 0.0;
                double maQzz = 0.0;
                if (IsMagneticField)
                {
                    // 磁界
                    maPxx = 1.0 / ma.Epxx;
                    maPyy = 1.0 / ma.Epyy;
                    maPzz = 1.0 / ma.Epzz;
                    maQxx = ma.Muxx;
                    maQyy = ma.Muyy;
                    maQzz = ma.Muzz;
                }
                else
                {
                    // 電界
                    maPxx = 1.0 / ma.Muxx;
                    maPyy = 1.0 / ma.Muyy;
                    maPzz = 1.0 / ma.Muzz;
                    maQxx = ma.Epxx;
                    maQyy = ma.Epyy;
                    maQzz = ma.Epzz;
                }

                IntegrationPoints ip = TetrahedronFE.GetIntegrationPoints(World.TetIntegrationPointCount);//Point5
                for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                {
                    double[] L = ip.Ls[ipPt];
                    double[][] N = tetFE.CalcEdgeN(L);
                    double[][] rotN = tetFE.CalcRotEdgeN(L);

                    double detJ = tetFE.GetDetJacobian(L);
                    double weight = ip.Weights[ipPt];
                    double detJWeight = (1.0 / 6.0) * weight * detJ;

                    for (int row = 0; row < elemEdgeNodeCnt; row++)
                    {
                        int rowEdgeNodeId = edgeNodes[row];
                        if (rowEdgeNodeId == -1)
                        {
                            continue;
                        }
                        double rowEdgeSgn = isReverses[row] ? -1.0 : 1.0;
                        for (int col = 0; col < elemEdgeNodeCnt; col++)
                        {
                            int colEdgeNodeId = edgeNodes[col];
                            if (colEdgeNodeId == -1)
                            {
                                continue;
                            }
                            double colEdgeSgn = isReverses[col] ? -1.0 : 1.0;

                            double kValue = detJWeight * (
                                rotN[row][0] * maPxx * rotN[col][0] +
                                rotN[row][1] * maPyy * rotN[col][1] +
                                rotN[row][2] * maPzz * rotN[col][2]);
                            double c0 = Constants.C0;
                            double mValue = detJWeight * (1.0 / (c0 * c0)) * (
                                N[row][0] * maQxx * N[col][0] +
                                N[row][1] * maQyy * N[col][1] +
                                N[row][2] * maQzz * N[col][2]);
                            double aValue =
                                (1.0 / (dt * dt)) * mValue +
                                beta * kValue;
                            A[rowEdgeNodeId, colEdgeNodeId] +=
                                rowEdgeSgn * colEdgeSgn * aValue;
                        }
                    }
                }
            }
        }

        private void _CalcB()
        {
            double dt = TimeStep;
            double beta = NewmarkBeta;

            IList<uint> feIds = World.GetTetrahedronFEIds(QuantityId);
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = World.GetTetrahedronFE(QuantityId, feId);
                uint elemNodeCnt = tetFE.NodeCount;
                uint elemEdgeNodeCnt = tetFE.EdgeCount;
                int[] edgeIds = new int[elemEdgeNodeCnt];
                bool[] isReverses = new bool[elemEdgeNodeCnt];
                int[] edgeNodes = new int[elemEdgeNodeCnt];
                for (int iENode = 0; iENode < elemEdgeNodeCnt; iENode++)
                {
                    int[] coIds = tetFE.EdgeCoordIdss[iENode];
                    bool isReverse;
                    int edgeId = World.GetEdgeIdFromCoords(QuantityId, coIds[0], coIds[1], out isReverse);
                    int edgeNodeId = World.Edge2EdgeNode(QuantityId, edgeId);

                    edgeIds[iENode] = edgeId;
                    isReverses[iENode] = isReverse;
                    edgeNodes[iENode] = edgeNodeId;
                }

                Material ma0 = World.GetMaterial(tetFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is DielectricMaterial);
                var ma = ma0 as DielectricMaterial;
                // TODO: 複素誘電率対応する
                double maPxx = 0.0;
                double maPyy = 0.0;
                double maPzz = 0.0;
                double maQxx = 0.0;
                double maQyy = 0.0;
                double maQzz = 0.0;
                if (IsMagneticField)
                {
                    // 磁界
                    maPxx = 1.0 / ma.Epxx;
                    maPyy = 1.0 / ma.Epyy;
                    maPzz = 1.0 / ma.Epzz;
                    maQxx = ma.Muxx;
                    maQyy = ma.Muyy;
                    maQzz = ma.Muzz;
                }
                else
                {
                    // 電界
                    maPxx = 1.0 / ma.Muxx;
                    maPyy = 1.0 / ma.Muyy;
                    maPzz = 1.0 / ma.Muzz;
                    maQxx = ma.Epxx;
                    maQyy = ma.Epyy;
                    maQzz = ma.Epzz;
                }

                IntegrationPoints ip = TetrahedronFE.GetIntegrationPoints(World.TetIntegrationPointCount);//Point5
                for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                {
                    double[] L = ip.Ls[ipPt];
                    double[][] N = tetFE.CalcEdgeN(L);
                    double[][] rotN = tetFE.CalcRotEdgeN(L);

                    double detJ = tetFE.GetDetJacobian(L);
                    double weight = ip.Weights[ipPt];
                    double detJWeight = (1.0 / 6.0) * weight * detJ;

                    for (int row = 0; row < elemEdgeNodeCnt; row++)
                    {
                        int rowEdgeNodeId = edgeNodes[row];
                        if (rowEdgeNodeId == -1)
                        {
                            continue;
                        }
                        double rowEdgeSgn = isReverses[row] ? -1.0 : 1.0;
                        for (int col = 0; col < elemEdgeNodeCnt; col++)
                        {
                            int colEdgeNodeId = edgeNodes[col];
                            if (colEdgeNodeId == -1)
                            {
                                continue;
                            }
                            double colEdgeSgn = isReverses[col] ? -1.0 : 1.0;

                            double kValue = detJWeight * (
                                rotN[row][0] * maPxx * rotN[col][0] +
                                rotN[row][1] * maPyy * rotN[col][1] +
                                rotN[row][2] * maPzz * rotN[col][2]);
                            double c0 = Constants.C0;
                            double mValue = detJWeight * (1.0 / (c0 * c0)) * (
                                N[row][0] * maQxx * N[col][0] +
                                N[row][1] * maQyy * N[col][1] +
                                N[row][2] * maQzz * N[col][2]);
                            double bValue =
                                mValue * (1.0 / (dt * dt)) * (
                                2.0 * EPrev[colEdgeNodeId] - EPrev2[colEdgeNodeId]) +
                                kValue * (
                                (2.0 * beta - 1.0) * EPrev[colEdgeNodeId] - beta * EPrev2[colEdgeNodeId]);
                            B[rowEdgeNodeId] +=
                                rowEdgeSgn * colEdgeSgn * bValue;
                        }
                    }
                }
            }
        }

        private void CalcEigen(
            int portId, double srcFreq,
            out System.Numerics.Complex[] betas,
            out System.Numerics.Complex[][] etEVecs,
            out System.Numerics.Complex[][] ezEVecs,
            out EMWaveguide2DEigenFEMForPort eigenFEM)
        {
            eigenFEM = new EMWaveguide2DEigenFEMForPort(World, (uint)portId);

            eigenFEM.IsMagneticField = IsMagneticField;
            eigenFEM.Frequency = srcFreq;
            eigenFEM.Solve();
            betas = eigenFEM.Betas;
            etEVecs = eigenFEM.EtEVecs;
            ezEVecs = eigenFEM.EzEVecs;
        }

        private void _CalcBSrc()
        {
            double dt = TimeStep;
            // 周波数
            double srcFreq = SrcFrequency;
            // 角周波数
            double srcOmega = 2.0 * Math.PI * srcFreq;

            int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1;  // 参照面、励振源分引く

            //--------------------------------------------------------------
            // 励振源
            //--------------------------------------------------------------
            {
                int portId = portCnt + RefPortCount; // ポートリストの最後の要素が励振境界
                var eigenFEM = EigenFEMs[portId];
                var etVec = eigenFEM.EtEVecs[0];
                var ezVec = eigenFEM.EzEVecs[0];
                int edgeNodeCntB = etVec.Length;
                int nodeCntB = ezVec.Length;
                var Stz = eigenFEM.Stz;
                var Rtt = eigenFEM.Rtt;
                double mur = 1.0;//!!!!! FIXME
                System.Numerics.Complex srcBeta = eigenFEM.Betas[0];

                double srcU0 = 0.0;
                double srcU1 = 0.0;
                double srcU2 = 0.0;

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
                            srcU0 = Math.Sin(srcOmega * (n) * dt) *
                                Math.Exp(-1.0 * ((n) * dt - GaussianT0) * ((n) * dt - GaussianT0) /
                                (2.0 * GaussianTp * GaussianTp));
                            srcU1 = Math.Sin(srcOmega * (n - 1) * dt) *
                                Math.Exp(-1.0 * ((n - 1) * dt - GaussianT0) * ((n - 1) * dt - GaussianT0) /
                                (2.0 * GaussianTp * GaussianTp));
                            srcU2 = Math.Sin(srcOmega * (n + 1) * dt) *
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

                var vec1 = new System.Numerics.Complex[edgeNodeCntB];
                var vec2 = new System.Numerics.Complex[edgeNodeCntB];
                for (int i = 0; i < edgeNodeCntB; i++)
                {
                    for (int j = 0; j < nodeCntB; j++)
                    {
                        vec1[i] += Stz[i, j] * ezVec[j];
                    }
                }
                for (int i = 0; i < edgeNodeCntB; i++)
                {
                    for (int j = 0; j < edgeNodeCntB; j++)
                    {
                        vec2[i] += Rtt[i, j] * etVec[j];
                    }
                }
                var vec3 = new System.Numerics.Complex[edgeNodeCntB];
                for (int i = 0; i < edgeNodeCntB; i++)
                {
                    vec3[i] =
                        (-2.0 * Constants.Mu0) / (-Constants.Mu0 * mur) *
                        (vec1[i] + System.Numerics.Complex.ImaginaryOne * srcBeta * vec2[i]);
                }

                {
                    double[] srcVec = new double[edgeNodeCntB];
                    for (int i = 0; i < edgeNodeCntB; i++)
                    {
                        srcVec[i] =
                            vec3[i].Real * srcU0 +
                            (1.0 / srcOmega) * vec3[i].Imaginary * (srcU2 - srcU1) / (2.0 * dt);
                    }
                    for (int edgeNodeIdB = 0; edgeNodeIdB < edgeNodeCntB; edgeNodeIdB++)
                    {
                        int edgeId = World.PortEdgeNode2Edge(QuantityId, (uint)portId, edgeNodeIdB);
                        int edgeNodeId = World.Edge2EdgeNode(QuantityId, edgeId);
                        B[edgeNodeId] += srcVec[edgeNodeIdB];
                    }
                }
            }
        }

        private void CalcCoordExyz(
            double[] E,
            out double[] coordExyz)
        {
            uint quantityId = 0;
            int edgeNodeCnt = (int)World.GetEdgeNodeCount(quantityId);
            int nodeCnt = (int)World.GetNodeCount(quantityId);
            int coCnt = (int)World.GetCoordCount(quantityId);

            int dof = 3; // x, y, z成分
            coordExyz = new double[coCnt * dof];
            int[] coordValueCnt = new int[coCnt];

            IList<uint> feIds = World.GetTetrahedronFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = World.GetTetrahedronFE(quantityId, feId);
                uint elemEdgeNodeCnt = tetFE.EdgeCount;
                uint elemNodeCnt = tetFE.NodeCount;
                int[] edgeIds = new int[elemEdgeNodeCnt];
                bool[] isReverses = new bool[elemEdgeNodeCnt];
                int[] edgeNodes = new int[elemEdgeNodeCnt];
                for (int iENode = 0; iENode < elemEdgeNodeCnt; iENode++)
                {
                    int[] coIds = tetFE.EdgeCoordIdss[iENode];
                    bool isReverse;
                    int edgeId = World.GetEdgeIdFromCoords(quantityId, coIds[0], coIds[1], out isReverse);
                    int edgeNodeId = World.Edge2EdgeNode(quantityId, edgeId);

                    edgeIds[iENode] = edgeId;
                    isReverses[iENode] = isReverse;
                    edgeNodes[iENode] = edgeNodeId;
                }

                // 辺方向成分の値
                double[] et = new double[elemEdgeNodeCnt];
                for (int iENode = 0; iENode < elemEdgeNodeCnt; iENode++)
                {
                    int edgeNodeId = edgeNodes[iENode];
                    if (edgeNodeId == -1)
                    {
                        continue;
                    }
                    double sgn = isReverses[iENode] ? -1.0 : 1.0;
                    et[iENode] = sgn * E[edgeNodeId];
                }

                // 節点のL
                double[][] nodeL = null;
                if (elemNodeCnt == 4)
                {
                    nodeL = new double[4][]
                    {
                        new double[4] { 1.0, 0.0, 0.0, 0.0 },
                        new double[4] { 0.0, 1.0, 0.0, 0.0 },
                        new double[4] { 0.0, 0.0, 1.0, 0.0 },
                        new double[4] { 0.0, 0.0, 0.0, 1.0 }
                    };
                }
                else if (elemNodeCnt == 10)
                {
                    nodeL = new double[10][]
                    {
                        new double[4] { 1.0, 0.0, 0.0, 0.0 },
                        new double[4] { 0.0, 1.0, 0.0, 0.0 },
                        new double[4] { 0.0, 0.0, 1.0, 0.0 },
                        new double[4] { 0.0, 0.0, 0.0, 1.0 },
                        new double[4] { 0.5, 0.5, 0.0, 0.0 },
                        new double[4] { 0.0, 0.5, 0.5, 0.0 },
                        new double[4] { 0.5, 0.0, 0.5, 0.0 },
                        new double[4] { 0.5, 0.0, 0.0, 0.5 },
                        new double[4] { 0.0, 0.5, 0.0, 0.5 },
                        new double[4] { 0.0, 0.0, 0.5, 0.5 }
                    };
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
                System.Diagnostics.Debug.Assert(nodeL.Length == elemNodeCnt);
                // 節点におけるEx,Ey,Ezを求める
                double[][] exyz = new double[elemNodeCnt][];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    double[] L = nodeL[iNode];
                    double[][] N = tetFE.CalcEdgeN(L);

                    exyz[iNode] = new double[dof];
                    for (int kENode = 0; kENode < elemEdgeNodeCnt; kENode++)
                    {
                        for (int idim = 0; idim < dof; idim++)
                        {
                            exyz[iNode][idim] += N[kENode][idim] * et[kENode];
                        }
                    }
                }

                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    //int nodeId = nodes[iNode];
                    int coId = tetFE.NodeCoordIds[iNode];
                    for (int idim = 0; idim < dof; idim++)
                    {
                        coordExyz[coId * dof + idim] += exyz[iNode][idim];
                    }
                    coordValueCnt[coId]++;
                }
            }

            for (int coId = 0; coId < coCnt; coId++)
            {
                for (int idim = 0; idim < dof; idim++)
                {
                    coordExyz[coId * dof + idim] /= coordValueCnt[coId];
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

            if (RefPortCount > 0)
            {
                CalcSParameterAtRefPorts(freqDomainAmpsInc, out freqs, out freqDomainAmpss, out Sss);
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
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
                int timeCnt = RefTimeEsss[refIndex].Count;
                int edgeNodeCntB = RefTimeEsss[refIndex][0].Length;
                double[][] datasRef = new double[edgeNodeCntB][];
                for (int edgeNodeIdB = 0; edgeNodeIdB < edgeNodeCntB; edgeNodeIdB++)
                {
                    datasRef[edgeNodeIdB] = new double[timeCnt];
                    for (int timeIndex = 0; timeIndex < timeCnt; timeIndex++)
                    {
                        datasRef[edgeNodeIdB][timeIndex] = RefTimeEsss[refIndex][timeIndex][edgeNodeIdB];
                    }
                }
                datasRefs[refIndex] = datasRef;
            }

            int dataCnt = datasRefs[0][0].Length;
            double dt = TimeStep;

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
                int edgeNodeCntB = datass.Length;
                System.Diagnostics.Debug.Assert(datass[0].Length == dataCnt);
                System.Numerics.Complex[][] freqDomainDatass = new System.Numerics.Complex[dataCnt][];
                for (int freqIndex = 0; freqIndex < dataCnt; freqIndex++)
                {
                    freqDomainDatass[freqIndex] = new System.Numerics.Complex[edgeNodeCntB];
                }
                for (int edgeNodeIdB = 0; edgeNodeIdB < edgeNodeCntB; edgeNodeIdB++)
                {
                    double[] datas = datass[edgeNodeIdB];
                    double[] _freqs;
                    System.Numerics.Complex[] _freqDomainDatas;
                    IvyFEM.FFT.Functions.DoFFT(times, datas, out _freqs, out _freqDomainDatas);
                    System.Diagnostics.Debug.Assert(_freqs.Length == dataCnt);
                    if (refIndex == 0 && edgeNodeIdB == 0)
                    {
                        freqs = _freqs;
                    }
                    for (int freqIndex = 0; freqIndex < dataCnt; freqIndex++)
                    {
                        System.Numerics.Complex data = _freqDomainDatas[freqIndex];
                        freqDomainDatass[freqIndex][edgeNodeIdB] = data;
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
                    // モード
                    System.Numerics.Complex[] betas;
                    System.Numerics.Complex[][] etEVecs;
                    System.Numerics.Complex[][] ezEVecs;
                    EMWaveguide2DEigenFEMForPort eigenFEM;
                    CalcEigen(
                        portId, freq,
                        out betas, out etEVecs, out ezEVecs, out eigenFEM);

                    // 振幅分布
                    System.Numerics.Complex[] freqEt = freqDomainDatass[freqIndex];
                    // モード振幅の算出
                    int iMode = 0; // 基本モード
                    System.Numerics.Complex b = eigenFEM.CalcModeAmp(
                        omega, iMode, betas, etEVecs, ezEVecs, freqEt);
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
