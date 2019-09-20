using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class EMWaveguide2DHPlaneHigherOrderABCTDFEM : FEM
    {
        public uint QuantityId { get; private set; } = 0;
        /// <summary>
        /// 吸収境界条件の次数
        /// </summary>
        public int ABCOrder { get; set; } = 5;
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
        /// 観測点の頂点ID
        /// </summary>
        public IList<uint> RefVIds { get; private set; } = new List<uint>();

        /// <summary>
        /// [A]
        /// </summary>
        private IvyFEM.Lapack.DoubleMatrix A = null;
        // 逆行列を計算するのでIvyFEM.Lapack.DoubleMatrixを使用
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
        private IList<double> SrcBetaXs = null;
        /// <summary>
        /// 境界の界のモード分布(ポート単位)
        /// </summary>
        private IList<double[]> SrcProfiles = null;

        /// <summary>
        /// 吸収境界波の速度リスト(ポート単位)
        /// </summary>
        private IList<double[]> Velos = null;
        /// <summary>
        /// 吸収境界 b0_abc
        /// </summary>
        private IList<double> B0Abcs = null;
        /// <summary>
        /// 吸収境界 a_abc
        /// </summary>
        private IList<double[]> AAbcs = null;
        /// <summary>
        /// 吸収境界 b_abc
        /// </summary>
        private IList<double[]> BAbcs = null;

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
        /// 観測点1の電界(時間変化リスト)
        /// </summary>
        public IList<double> TimeEzsPort1 { get; private set; } = null;
        /// <summary>
        /// 観測点2の電界(時間変化リスト)
        /// </summary>
        public IList<double> TimeEzsPort2 { get; private set; }  = null;

        public EMWaveguide2DHPlaneHigherOrderABCTDFEM
            (FEWorld world)
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
                TimeEzsPort1 = new List<double>();
                TimeEzsPort2 = new List<double>();
            }

            // 時刻の取得
            if (TimeIndex == 0)
            {
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
            //t = System.Environment.TickCount;
            {
                ////----------------------------------
                //double[] X;
                //Solver.DoubleSolve(out X, A, B);
                //EzPz = X;
                ////----------------------------------

                // 逆行列を用いる
                EzPz = A * B;
            }
            //System.Diagnostics.Debug.WriteLine("Solve t = " + (System.Environment.TickCount - t));

            // ABC分を除いた電界を取得
            Ez = new double[nodeCnt];
            for (int nodeId = 0; nodeId < nodeCnt; nodeId++)
            {
                Ez[nodeId] = EzPz[nodeId];
            }

            // 観測点
            int portCnt = (int)World.GetPortCount(QuantityId) - 1; // 励振分を引く
            for (int portId = 0; portId < portCnt; portId++)
            {
                uint vId = RefVIds[portId];
                IList<int> coIds = World.GetCoordIdsFromCadId(QuantityId, vId, CadElementType.Vertex);
                System.Diagnostics.Debug.Assert(coIds.Count == 1);
                int coId = coIds[0];
                int nodeId = World.Coord2Node(QuantityId, coId);
                double fValue = Ez[nodeId];
                if (portId == 0)
                {
                    TimeEzsPort1.Add(fValue);
                }
                else if (portId == 1)
                {
                    TimeEzsPort2.Add(fValue);
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
            }

            if ((TimeIndex + 1) % 50 == 0)
            {
                System.Diagnostics.Debug.WriteLine("timeIndex: {0}", (TimeIndex + 1));
            }
        }

        private void CalcA()
        {
            System.Diagnostics.Debug.Assert(TimeIndex == 0);
            A = null;
            Qbs = new List<IvyFEM.Lapack.DoubleMatrix>();
            Rbs = new List<IvyFEM.Lapack.DoubleMatrix>();
            Tbs = new List<IvyFEM.Lapack.DoubleMatrix>();
            SrcBetaXs = new List<double>();
            SrcProfiles = new List<double[]>();
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

            int portCnt = (int)World.GetPortCount(QuantityId) - 1; // ポート3を除く(励振源)

            // 境界の界に導波路開口条件を追加
            for (int portId = 0; portId < (portCnt + 1); portId++)
            {
                IvyFEM.Lapack.DoubleMatrix ryy1D = null;
                IvyFEM.Lapack.DoubleMatrix txx1D = null;
                IvyFEM.Lapack.DoubleMatrix uzz1D = null;
                System.Numerics.Complex[] betas = null;
                System.Numerics.Complex[][] eVecs = null;
                var eigen1DFEM = new EMWaveguide1DEigenFEM(World, QuantityId, (uint)portId);
                eigen1DFEM.Frequency = srcFreq;
                eigen1DFEM.Solve();
                ryy1D = eigen1DFEM.Ryy;
                txx1D = eigen1DFEM.Txx;
                uzz1D = eigen1DFEM.Uzz;
                betas = eigen1DFEM.Betas;
                eVecs = eigen1DFEM.EzEVecs;
                int nodeCntB = ryy1D.RowLength;
                Qbs.Add(ryy1D);
                Rbs.Add(txx1D);
                Tbs.Add(uzz1D);

                // 基本モード
                uint iMode = 0;
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
            }

            /////////////////////////////////////////////////////////

            //------------------------------------------------------
            // 電界
            //-----------------------------------------------------
            int nodeCnt = (int)World.GetNodeCount(QuantityId);
            int nodeCntPlusABC = nodeCnt;
            for (int portId = 0; portId < portCnt; portId++)
            {
                int nodeCntB = (int)World.GetPortNodeCount(QuantityId, (uint)portId);
                System.Diagnostics.Debug.Assert(Qbs[portId].RowLength == nodeCntB);
                nodeCntPlusABC += nodeCntB * (ABCOrder - 1);
            }
            EzPz = new double[nodeCntPlusABC];
            EzPzPrev = new double[nodeCntPlusABC];
            EzPzPrev2 = new double[nodeCntPlusABC];

            //------------------------------------------------------
            // 全体係数行列の作成
            //------------------------------------------------------
            A = new IvyFEM.Lapack.DoubleMatrix(nodeCntPlusABC, nodeCntPlusABC);
            double dt = TimeDelta;
            for (int rowNodeId = 0; rowNodeId < nodeCnt; rowNodeId++)
            {
                for (int colNodeId = 0; colNodeId < nodeCnt; colNodeId++)
                {
                    A[rowNodeId, colNodeId] =
                        (1.0 / (dt * dt)) * M[rowNodeId, colNodeId] +
                        NewmarkBeta * K[rowNodeId, colNodeId];
                }
            }

            // 吸収境界パラメータ
            Velos = new List<double[]>();
            B0Abcs = new List<double>();
            AAbcs = new List<double[]>();
            BAbcs = new List<double[]>();

            for (int portId = 0; portId < portCnt; portId++)
            {
                double srcBetaX = SrcBetaXs[portId];
                double vpx = srcOmega / srcBetaX;
                double[] velo = new double[ABCOrder];
                for (int order = 0; order < ABCOrder; order++)
                {
                    if (GaussianType == GaussianType.Normal)
                    {
                        velo[order] = Constants.C0; // ガウシアンパルスの場合
                    }
                    else if (GaussianType == GaussianType.SinModulation)
                    {
                        velo[order] = vpx; // 正弦波変調ガウシアンパルスの場合
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                }
                Velos.Add(velo);
            }
            for (int portId = 0; portId < portCnt; portId++)
            {
                double[] velo = Velos[portId];
                double B0Abc = 1.0 / velo[0];
                double[] AAbc = new double[ABCOrder - 1]; // 一様媒質の場合
                double[] BAbc = new double[ABCOrder - 1];
                for (int order = 0; order < (ABCOrder - 1); order++)
                {
                    AAbc[order] = 1.0 / (velo[order] * velo[order]) - 1.0 / (Constants.C0 * Constants.C0);
                    BAbc[order] = 1.0 / velo[order] + 1.0 / velo[order + 1];
                }
                B0Abcs.Add(B0Abc);
                AAbcs.Add(AAbc);
                BAbcs.Add(BAbc);
            }

            // 吸収境界
            for (int portId = 0; portId < portCnt; portId++)
            {
                double B0Abc = B0Abcs[portId];

                var Qb = Qbs[portId];
                int nodeCntB = Qb.RowLength;
                for (int rowNodeIdB = 0; rowNodeIdB < nodeCntB; rowNodeIdB++)
                {
                    int rowCoId = World.PortNode2Coord(QuantityId, (uint)portId, rowNodeIdB);
                    int rowNodeId = World.Coord2Node(QuantityId, rowCoId);

                    for (int colNodeIdB = 0; colNodeIdB < nodeCntB; colNodeIdB++)
                    {
                        int colCoId = World.PortNode2Coord(QuantityId, (uint)portId, colNodeIdB);
                        int colNodeId = World.Coord2Node(QuantityId, colCoId);

                        // C
                        A[rowNodeId, colNodeId] +=
                            (B0Abc / (2.0 * dt)) * Qb[rowNodeIdB, colNodeIdB];
                    }
                }
            }

            // Φ1の開始位置
            int[] Pz1NodeIds = new int[portCnt];
            for (int portId = 0; portId < portCnt; portId++)
            {
                if (portId == 0)
                {
                    Pz1NodeIds[portId] = nodeCnt;
                }
                else
                {
                    int nodeCntBPrevPort = (int)World.GetPortNodeCount(QuantityId, (uint)(portId - 1));
                    Pz1NodeIds[portId] =
                        Pz1NodeIds[portId - 1] + nodeCntBPrevPort * (ABCOrder - 1);
                }
            }

            // Ez - Φ1
            if (ABCOrder > 1)
            {
                for (int portId = 0; portId < portCnt; portId++)
                {
                    var Qb = Qbs[portId];
                    int nodeCntB = Qb.RowLength;

                    for (int rowNodeIdB = 0; rowNodeIdB < nodeCntB; rowNodeIdB++)
                    {
                        int rowCoId = World.PortNode2Coord(QuantityId, (uint)portId, rowNodeIdB);
                        int rowNodeId = World.Coord2Node(QuantityId, rowCoId);

                        for (int colNodeIdB = 0; colNodeIdB < nodeCntB; colNodeIdB++)
                        {
                            // Φ1
                            int colNodeId1 = colNodeIdB + Pz1NodeIds[portId];
                            // G
                            A[rowNodeId, colNodeId1] = -1.0 * NewmarkBeta * Qb[rowNodeIdB, colNodeIdB];
                        }
                    }
                }
            }

            // Φ1 ～ Φ(ABC_order - 1)
            for (int portId = 0; portId < portCnt; portId++)
            {
                var Qb = Qbs[portId];
                var Rb = Rbs[portId];
                var Tb = Tbs[portId];
                double[] velo = Velos[portId];
                //double B0Abc = B0Abcs[portId];
                double[] AAbc = AAbcs[portId]; // 一様媒質のとき
                double[] BAbc = BAbcs[portId];
                int nodeCntB = Qb.RowLength;

                // Φ(order)に関する式
                for (int order = 0; order < (ABCOrder - 1); order++)
                {
                    for (int rowNodeIdB = 0; rowNodeIdB < nodeCntB; rowNodeIdB++)
                    {
                        // order
                        int rowNodeId1 = rowNodeIdB + nodeCntB * order + Pz1NodeIds[portId];

                        for (int colNodeIdB = 0; colNodeIdB < nodeCntB; colNodeIdB++)
                        {
                            int colNodeId1 = 0;
                            // Φ(order)
                            colNodeId1 = colNodeIdB + nodeCntB * order + Pz1NodeIds[portId];
                            // Cj
                            A[rowNodeId1, colNodeId1] =
                                (BAbc[order] / (2.0 * dt)) * Qb[rowNodeIdB, colNodeIdB];

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
                                colNodeId1 = colNodeIdB + nodeCntB * (order - 1) + Pz1NodeIds[portId];
                            }
                            // 一様媒質の場合
                            A[rowNodeId1, colNodeId1] =
                                // Pj
                                -(AAbc[order] / (dt * dt)) * Qb[rowNodeIdB, colNodeIdB] +
                                // Qj
                                // Note: Rb + λQb (λ = 0)
                                NewmarkBeta * Rb[rowNodeIdB, colNodeIdB];
                            /*
                            // 媒質定数がy方向に変化する場合
                            A[rowNodeId1, colNodeId1] =
                                -1.0 * (1.0 / (velo[order] * velo[order] * dt * dt)) * Qb[rowNodeIdB, colNodeIdB]
                                + 1.0 * (1.0 / (Constants.C0 * Constants.C0 * dt * dt)) * Tb[rowNodeIdB, colNodeIdB]
                                + NewmarkBeta * Rb[rowNodeIdB, colNodeIdB];
                            */

                            // Φ(order + 1)
                            colNodeId1 = 0;
                            if (order == (ABCOrder - 2))
                            {
                                // なし
                            }
                            else
                            {
                                colNodeId1 = colNodeIdB + nodeCntB * (order + 1) + Pz1NodeIds[portId];
                                // Rj
                                A[rowNodeId1, colNodeId1] =
                                    -1.0 * NewmarkBeta * Qb[rowNodeIdB, colNodeIdB];
                            }
                        }
                    }
                }
            }

            // 逆行列を計算
            System.Diagnostics.Debug.WriteLine("calc [A]-1");
            A = IvyFEM.Lapack.DoubleMatrix.Inverse(A);
            System.Diagnostics.Debug.WriteLine("calc [A]-1 done");

            System.Diagnostics.Debug.WriteLine("ABCOrder:{0}", ABCOrder);
            System.Diagnostics.Debug.WriteLine("dt:{0}", TimeDelta);
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

                        // clapack形式の行列格納方法で格納
                        K[rowNodeId, colNodeId] += kValue;
                        M[rowNodeId, colNodeId] += mValue;
                    }
                }
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
            int portCnt = (int)World.GetPortCount(QuantityId) - 1;  // 励振源分引く
            // 吸収境界
            // Φ1の開始位置
            int[] Pz1NodeIds = new int[portCnt];
            for (int portId = 0; portId < portCnt; portId++)
            {
                if (portId == 0)
                {
                    Pz1NodeIds[portId] = nodeCnt;
                }
                else
                {
                    int nodeCntBPrevPort = (int)World.GetPortNodeCount(QuantityId, (uint)(portId - 1));
                    Pz1NodeIds[portId] = Pz1NodeIds[portId - 1] + nodeCntBPrevPort * (ABCOrder - 1);
                }
            }            
            
            // Ez - Φ1
            for (int portId = 0; portId < portCnt; portId++)
            {
                var Qb = Qbs[portId];
                int nodeCntB = Qb.RowLength;
                double B0Abc = B0Abcs[portId];

                double[] workEzPrev2 = new double[nodeCntB];
                double[] workPzOrder2Prev = new double[nodeCntB];
                double[] workPzOrder2Prev2 = new double[nodeCntB];
                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                {
                    int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                    int nodeId = World.Coord2Node(QuantityId, coId);
                    workEzPrev2[nodeIdB] = EzPzPrev2[nodeId];
                }
                
                if (ABCOrder > 1)
                {
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        // Φ1
                        int nodeId = nodeIdB + Pz1NodeIds[portId];
                        workPzOrder2Prev[nodeIdB] = EzPzPrev[nodeId];
                        workPzOrder2Prev2[nodeIdB] = EzPzPrev2[nodeId];
                    }
                }

                double[] vecQb = new double[nodeCntB];
                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                {
                    // C
                    vecQb[nodeIdB] = (B0Abc / (2.0 * dt)) * workEzPrev2[nodeIdB];
                    if (ABCOrder > 1)
                    {
                        vecQb[nodeIdB] +=
                            // G
                            (1.0 - 2.0 * NewmarkBeta) * workPzOrder2Prev[nodeIdB] +
                            NewmarkBeta * workPzOrder2Prev2[nodeIdB];
                    }
                }
                vecQb = Qb * vecQb;

                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                {
                    int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                    int nodeId = World.Coord2Node(QuantityId, coId);
                    B[nodeId] += vecQb[nodeIdB];
                }
            }

            // Φ1 ～ Φ(ABC_order - 1)
            for (int portId = 0; portId < portCnt; portId++)
            {
                var Qb = Qbs[portId];
                var Rb = Rbs[portId];
                var Tb = Tbs[portId];
                double[] velo = Velos[portId];
                //double B0Abc = B0Abcs[portId];
                double[] AAbc = AAbcs[portId]; // 一様媒質の場合
                double[] BAbc = BAbcs[portId];
                int nodeCntB = Qb.RowLength;

                // Φorderに関する式
                for (int order = 0; order < (ABCOrder - 1); order++)
                {
                    double[] workPzOrder0Prev2 = new double[nodeCntB];
                    double[] workPzOrder1Prev = new double[nodeCntB];
                    double[] workPzOrder1Prev2 = new double[nodeCntB];
                    double[] workPzOrder2Prev = new double[nodeCntB];
                    double[] workPzOrder2Prev2 = new double[nodeCntB];
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int nodeId = 0;
                        // Φorder
                        nodeId = nodeIdB + nodeCntB * (order) + Pz1NodeIds[portId];
                        workPzOrder0Prev2[nodeIdB] = EzPzPrev2[nodeId];

                        // Φorder-1
                        nodeId = 0;
                        if (order == 0)
                        {
                            // Φ0 (Ez)
                            int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                            nodeId = World.Coord2Node(QuantityId, coId);
                        }
                        else
                        {
                            nodeId = nodeIdB + nodeCntB * (order - 1) + Pz1NodeIds[portId];
                        }
                        workPzOrder1Prev[nodeIdB] = EzPzPrev[nodeId];
                        workPzOrder1Prev2[nodeIdB] = EzPzPrev2[nodeId];

                        // Φ(order + 1)
                        if (order == (ABCOrder - 2))
                        {
                            // なし
                            workPzOrder2Prev[nodeIdB] = 0.0;
                            workPzOrder2Prev2[nodeIdB] = 0.0;
                        }
                        else
                        {
                            nodeId = nodeIdB + nodeCntB * (order + 1) + Pz1NodeIds[portId];
                            workPzOrder2Prev[nodeIdB] = EzPzPrev[nodeId];
                            workPzOrder2Prev2[nodeIdB] = EzPzPrev2[nodeId];
                        }
                    }
                    double[] vecQb = new double[nodeCntB];
                    double[] vecRb = new double[nodeCntB];
                    double[] vecTb = new double[nodeCntB];
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        // 一様媒質の場合
                        vecQb[nodeIdB] =
                            // Cj
                            (BAbc[order] / (2.0 * dt)) * workPzOrder0Prev2[nodeIdB] +
                            // Pj
                            (AAbc[order] / (dt * dt)) * (
                            -2.0 * workPzOrder1Prev[nodeIdB] + workPzOrder1Prev2[nodeIdB]) +
                            // Rj
                            ((1.0 - 2.0 * NewmarkBeta) * workPzOrder2Prev[nodeIdB] +
                            NewmarkBeta * workPzOrder2Prev2[nodeIdB]);
                        /*
                        // 媒質定数がy方向に変化する場合 
                        vecQb[nodeIdB] =
                            (BAbc[order] / (2.0 * dt)) * workPzOrder0Prev2[nodeIdB] +
                            (1.0 / (velo[order] * velo[order] * dt * dt)) * 
                            (-2.0 * workPzOrder1Prev[nodeIdB] + workPzOrder1Prev2[nodeIdB]) +
                            ((1.0 - 2.0 * NewmarkBeta) * workPzOrder2Prev[nodeIdB] +
                            NewmarkBeta * workPzOrder2Prev2[nodeIdB]);
                        */

                        /*
                        // 一様媒質の場合
                        vecTb[nodeIdB] = 0;
                        */
                        /*
                        // 媒質定数がy方向に変化する場合
                        vecTb[nodeIdB] =
                            (-1.0 / (Constants.C0 * Constants.C0 * dt * dt)) *
                            (-2.0 * workPzOrder1Prev[nodeIdB] + workPzOrder1Prev2[nodeIdB]);
                        */

                        // 共用
                        vecRb[nodeIdB] =
                            // Qj
                            -1.0 * (
                            (1.0 - 2.0 * NewmarkBeta) * workPzOrder1Prev[nodeIdB] +
                            NewmarkBeta * workPzOrder1Prev2[nodeIdB]);
                    }
                    vecQb = Qb * vecQb;
                    vecRb = Rb * vecRb;
                    // Note: Rb + λQb (λ = 0)

                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int nodeId = nodeIdB + nodeCntB * (order) + Pz1NodeIds[portId];
                        B[nodeId] = vecQb[nodeIdB] + vecTb[nodeIdB] + vecRb[nodeIdB];
                    }
                }
            }

            //--------------------------------------------------------------
            // 励振源
            //--------------------------------------------------------------
            {
                int portId = portCnt; // ポートリストの最後の要素が励振境界
                var Qb = Qbs[portId];
                int nodeCntB = Qb.RowLength;
                double srcBetaX = SrcBetaXs[portId];
                double[] srcProfile = SrcProfiles[portId];

                double srcF1 = 0.0;
                double srcF2 = 0.0;
                double srcF0 = 0.0;

                int n = TimeIndex;
                // 周波数
                double srcFreq = SrcFrequency;
                // 角周波数
                double srcOmega = 2.0 * Math.PI * srcFreq;

                // ガウシアンパルス
                if ((n * dt) <= (2.0 * GaussianT0 + dt))
                {
                    if (GaussianType == GaussianType.Normal)
                    {
                        // ガウシアンパルス
                        // Note: veloはC0
                        srcF1 = Math.Exp(-1.0 * ((n + 1) * dt - GaussianT0) * ((n + 1) * dt - GaussianT0) /
                            (2.0 * GaussianTp * GaussianTp));
                        srcF2 = Math.Exp(-1.0 * ((n - 1) * dt - GaussianT0) * ((n - 1) * dt - GaussianT0) /
                            (2.0 * GaussianTp * GaussianTp));
                        srcF0 = Math.Exp(-1.0 * ((n) * dt - GaussianT0) * ((n) * dt - GaussianT0) /
                            (2.0 * GaussianTp * GaussianTp));
                    }
                    else if (GaussianType == GaussianType.SinModulation)
                    {
                        // 正弦波変調ガウシアンパルス
                        // Note: veloはvpx
                        // 正弦波変調ガウシアンパルス
                        srcF1 = Math.Cos(srcOmega * ((n + 1) * dt - GaussianT0)) *
                            Math.Exp(-1.0 * ((n + 1) * dt - GaussianT0) * ((n + 1) * dt - GaussianT0) /
                            (2.0 * GaussianTp * GaussianTp));
                        srcF2 = Math.Cos(srcOmega * ((n - 1) * dt - GaussianT0)) *
                            Math.Exp(-1.0 * ((n - 1) * dt - GaussianT0) * ((n - 1) * dt - GaussianT0) /
                            (2.0 * GaussianTp * GaussianTp));
                        srcF0 = Math.Cos(srcOmega * ((n) * dt - GaussianT0)) *
                            Math.Exp(-1.0 * ((n) * dt - GaussianT0) * ((n) * dt - GaussianT0) /
                            (2.0 * GaussianTp * GaussianTp));
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                }
                /*
                // 正弦波
                srcF1 = Math.Sin(srcOmega * (n + 1) * dt);
                srcF2 = Math.Sin(srcOmega * (n - 1) * dt);
                srcF0 = Math.Sin(srcOmega * (n) * dt);
                */
                //System.Diagnostics.Debug.WriteLine("Input: srcF2 = {0}", srcF2);

                // 境界積分
                // dF/dt (F = Ez)
                double[] srcFt = new double[nodeCntB];
                System.Diagnostics.Debug.Assert(srcProfile.Length == nodeCntB);
                for (int i = 0; i < nodeCntB; i++)
                {
                    double normalizeFactor = -1.0;
                    srcFt[i] = normalizeFactor * srcProfile[i] * (srcF1 - srcF2) / (2.0 * dt);
                }
                double vpx = srcOmega / srcBetaX;
                double[] vecQb = new double[nodeCntB];
                for (int i = 0; i < nodeCntB; i++)
                {
                    vecQb[i] = (-2.0 / vpx) * srcFt[i];
                }
                vecQb = Qb * vecQb;
                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                {
                    int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                    int nodeId = World.Coord2Node(QuantityId, coId);
                    B[nodeId] += vecQb[nodeIdB];
                }

                /*
                // 領域積分
                double[] srcFt2 = new double[nodeCnt];
                double[] srcF = new double[nodeCnt];
                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                {
                    int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                    int nodeId = World.Coord2Node(QuantityId, coId);
                    double normalizeFactor = -srcOmega * Constants.Mu0 / srcBetaX;
                    srcF[nodeId] = normalizeFactor * srcProfile[nodeIdB] * srcF0;
                    srcFt2[nodeId] = 
                        normalizeFactor * srcProfile[nodeIdB] *
                        (srcF1 - 2.0 * srcF0 + srcF2) / (dt * dt);
                }
                double[] vecM = M * srcFt2;
                double[] vecK = K * srcF;
                for (int i = 0; i < nodeCnt; i++)
                {
                    B[i] += vecM[i] + vecK[i];
                }
                */
            }
        }

        /// <summary>
        /// Sパラメータの計算
        /// </summary>
        /// <returns></returns>
        public void CalcSParameter(
            IList<double> TimeEzsPort1Inc,
            out IList<double> freqs, out IList<System.Numerics.Complex[]> Sss)
        {
            int dataCnt = TimeEzsPort1Inc.Count;
            double dt = TimeDelta;
            // 周波数のリスト
            freqs = null;

            double[] times = new double[dataCnt];
            for (int i = 0; i < dataCnt; i++)
            {
                times[i] = i * dt;
            }

            // 散乱パラメータを計算する
            double[] datasPort1Inc = TimeEzsPort1Inc.ToArray();
            double[] datasPort1 = TimeEzsPort1.ToArray();
            double[] datasPort2 = TimeEzsPort2.ToArray();
            System.Numerics.Complex[] S11s = null;
            System.Numerics.Complex[] S21s = null;
            {
                double[] _freqs;
                _CalcSParameter(
                    times,
                    datasPort1Inc,
                    datasPort1,
                    datasPort2,
                    out _freqs,
                    out S11s,
                    out S21s
                    );
                freqs = _freqs.ToList();
            }

            // ポート数
            //  励振面の分を引く
            int portCnt = (int)World.GetPortCount(QuantityId) - 1;
            // 散乱パラメータ(各周波数に対するS11とS21)のリスト
            System.Diagnostics.Debug.Assert(freqs.Count == dataCnt);
            Sss = new List<System.Numerics.Complex[]>();
            for (int i = 0; i < dataCnt; i++)
            {
                System.Numerics.Complex[] Ss = new System.Numerics.Complex[portCnt];
                for (int portId = 0; portId < portCnt; portId++)
                {

                    System.Numerics.Complex value = 0;
                    if (portId == 0)
                    {
                        value = S11s[i];
                    }
                    else if (portId == 1)
                    {
                        value = S21s[i];
                    }
                    Ss[portId] = value;
                }
                // 格納
                Sss.Add(Ss);
            }
        }

        private bool _CalcSParameter(
            double[] times,
            double[] datasPort1Inc,
            double[] datasPort1,
            double[] datasPort2,
            out double[] freqs,
            out System.Numerics.Complex[] S11s,
            out System.Numerics.Complex[] S21s
            )
        {
            int dataCnt = times.Length;

            // 入射波
            double[] freqsInc = null;
            System.Numerics.Complex[] freqDomainDatasInc = null;
            IvyFEM.FFT.Functions.DoFFT(times, datasPort1Inc, out freqsInc, out freqDomainDatasInc);

            // 反射波
            double[] workDatasR = new double[dataCnt];
            for (int i = 0; i < dataCnt; i++)
            {
                workDatasR[i] = datasPort1[i] - datasPort1Inc[i];
            }
            double[] freqsR = null;
            System.Numerics.Complex[] freqDomainDatasR = null;
            IvyFEM.FFT.Functions.DoFFT(times, workDatasR, out freqsR, out freqDomainDatasR);

            // 透過波
            double[] freqsT = null;
            System.Numerics.Complex[] freqDomainDatasT = null;
            IvyFEM.FFT.Functions.DoFFT(times, datasPort2, out freqsT, out freqDomainDatasT);

            freqs = freqsInc;
            System.Diagnostics.Debug.Assert(freqs.Length == dataCnt);
            S11s = new System.Numerics.Complex[dataCnt];
            S21s = new System.Numerics.Complex[dataCnt];
            for (int i = 0; i < dataCnt; i++)
            {
                System.Numerics.Complex S11 = freqDomainDatasR[i] / freqDomainDatasInc[i];
                System.Numerics.Complex S21 = freqDomainDatasT[i] / freqDomainDatasInc[i];
                S11s[i] = S11;
                S21s[i] = S21;
            }

            return true;
        }

    }
}
