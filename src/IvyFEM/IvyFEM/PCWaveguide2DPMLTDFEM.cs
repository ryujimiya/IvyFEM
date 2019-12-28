using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class PCWaveguide2DPMLTDFEM : FEM
    {
        //---------------------------------
        // ※ポートの設定順序
        // 参照(観測)面、励振源の順
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
        /// [A]
        /// </summary>
        private IvyFEM.Linear.DoubleSparseMatrix A = null;
        /// <summary>
        /// {b}
        /// </summary>
        private double[] B = null;
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
        /// 電界(1つ前)
        /// </summary>
        private double[] EzPrev = null;
        /// <summary>
        /// 電界(2つ前)
        /// </summary>
        private double[] EzPrev2 = null;
        /// <summary>
        /// ψx(1つ前→現在値に変化する) Px[要素番号][要素内節点番号]
        /// </summary>
        private double[][] Pxs = null;
        /// <summary>
        /// ψy(1つ前→現在値に変化する) Py[要素番号][要素内節点番号]
        /// </summary>
        private double[][] Pys = null;
        /// <summary>
        /// vx(1つ前→現在値に変化する) Vx[要素番号][要素内節点番号]
        /// </summary>
        private double[][] Vxs = null;
        /// <summary>
        /// vy(1つ前→現在値に変化する) Vy[要素番号][要素内節点番号]
        /// </summary>
        private double[][] Vys = null;

        /// <summary>
        /// 観測点の電界(時間変化リスト)
        /// </summary>
        public IList<IList<double[]>> RefTimeEzsss { get; private set; } = new List<IList<double[]>>();

        public PCWaveguide2DPMLTDFEM(FEWorld world)
        {
            World = world;
        }

        public override void Solve()
        {
            int t;
            int refPortCnt = (int)World.GetPortCount(QuantityId) - 1; // 励振面を引く

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
                else if (refPortCnt > 0)
                {
                    cnt = refPortCnt; 
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
            EzPrev.CopyTo(EzPrev2, 0);
            Ez.CopyTo(EzPrev, 0);
            Ez = null;

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
                double[] X;
                Solver.DoubleSolve(out X, A, B);
                Ez = X;
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
            else if (refPortCnt > 0)
            {
                for (int refIndex = 0; refIndex < refPortCnt; refIndex++)
                {
                    int portId = refIndex;
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

            IvyFEM.Linear.DoubleSparseMatrix _A;
            _A = null;
            Qbs = new List<IvyFEM.Lapack.DoubleMatrix>();
            Rbs = new List<IvyFEM.Lapack.DoubleMatrix>();
            Tbs = new List<IvyFEM.Lapack.DoubleMatrix>();
            SrcBetaXs = new List<double>();
            SrcProfiles = new List<System.Numerics.Complex[]>();
            SrcModifyProfiles = new List<System.Numerics.Complex[]>();
            Ez = null;
            EzPrev = null;
            EzPrev2 = null;

            //------------------------------------------------------
            // 電界
            //-----------------------------------------------------
            int nodeCnt = (int)World.GetNodeCount(QuantityId);
            Ez = new double[nodeCnt];
            EzPrev = new double[nodeCnt];
            EzPrev2 = new double[nodeCnt];
            int feCnt = World.GetTriangleFEIds(QuantityId).Count;
            Pxs = new double[feCnt][];
            Pys = new double[feCnt][];
            Vxs = new double[feCnt][];
            Vys = new double[feCnt][];

            //------------------------------------------------------
            // 全体係数行列の作成
            //------------------------------------------------------
            _A = new IvyFEM.Linear.DoubleSparseMatrix(nodeCnt, nodeCnt);
            CalcA(_A);

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

            int refPortCnt = (int)World.GetPortCount(QuantityId) - 1; //励振源を除く

            for (int portId = 0; portId < (refPortCnt + 1); portId++)
            {
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
                System.Diagnostics.Debug.Assert(fVec.Length  == nodeCntB);
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
            }

            A = _A;
        }

        private void CalcA(IvyFEM.Linear.DoubleSparseMatrix _A)
        {
            double dt = TimeDelta;

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
                System.Diagnostics.Debug.Assert(
                    ma0 is DielectricMaterial ||
                    ma0 is DielectricPMLMaterial);
                DielectricMaterial ma = null;
                DielectricPMLMaterial maPML = null;
                double epxx = 0;
                double epyy = 0;
                double epzz = 0;
                double muxx = 0;
                double muyy = 0;
                double muzz = 0;
                double rotAngle = 0.0;
                OpenTK.Vector2d rotOrigin = new OpenTK.Vector2d();
                if (ma0 is DielectricMaterial)
                {
                    ma = ma0 as DielectricMaterial;
                    epxx = ma.Epxx;
                    epyy = ma.Epyy;
                    epzz = ma.Epzz;
                    muxx = ma.Muxx;
                    muyy = ma.Muyy;
                    muzz = ma.Muzz;
                }
                else if (ma0 is DielectricPMLMaterial)
                {
                    maPML = ma0 as DielectricPMLMaterial;
                    epxx = maPML.Epxx;
                    epyy = maPML.Epyy;
                    epzz = maPML.Epzz;
                    muxx = maPML.Muxx;
                    muyy = maPML.Muyy;
                    muzz = maPML.Muzz;
                    rotOrigin = maPML.RotOriginPoint;
                    rotAngle = maPML.RotAngle;
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
                // 回転移動
                World.RotAngle = rotAngle;
                World.RotOrigin = new double[] { rotOrigin.X, rotOrigin.Y };
                double maPxx = 0;
                double maPyy = 0;
                double maQzz = 0;
                if (IsTMMode)
                {
                    // TMモード
                    maPxx = 1.0 / epxx;
                    maPyy = 1.0 / epyy;
                    maQzz = muzz;
                }
                else
                {
                    // TEモード
                    maPxx = 1.0 / muxx;
                    maPyy = 1.0 / muyy;
                    maQzz = epzz;
                }

                // 重心を求める
                OpenTK.Vector2d cPt;
                {
                    OpenTK.Vector2d[] vertexPts = new OpenTK.Vector2d[triFE.VertexCount];
                    for (int iVertex = 0; iVertex < triFE.VertexCount; iVertex++)
                    {
                        int coId = triFE.VertexCoordIds[iVertex];
                        double[] coord = World.GetCoord(QuantityId, coId);
                        vertexPts[iVertex] = new OpenTK.Vector2d(coord[0], coord[1]);
                    }
                    cPt = (vertexPts[0] + vertexPts[1] + vertexPts[2]) / 3.0;
                }

                // PML
                bool isXDirection = false;
                bool isYDirection = false;
                // X方向PML
                double sigmaX = 0.0;
                double c1PX = 0.0;
                double c2PX = 0.0;
                double c1VX = 0.0;
                // Y方向PML
                double sigmaY = 0.0;
                double c1PY = 0.0;
                double c2PY = 0.0;
                double c1VY = 0.0;
                if (maPML != null)
                {
                    isXDirection = maPML.IsXDirection();
                    isYDirection = maPML.IsYDirection();
                    if (isXDirection && !isYDirection)
                    {
                        maPML.CalcSigmaXForTD(cPt, dt, out sigmaX, out c1PX, out c2PX, out c1VX);
                    }
                    else if (isYDirection && !isXDirection)
                    {
                        maPML.CalcSigmaYForTD(cPt, dt, out sigmaY, out c1PY, out c2PY, out c1VY);
                    }
                    else if (isXDirection && isYDirection)
                    {
                        maPML.CalcSigmaXForTD(cPt, dt, out sigmaX, out c1PX, out c2PX, out c1VX);
                        maPML.CalcSigmaYForTD(cPt, dt, out sigmaY, out c1PY, out c2PY, out c1VY);
                    }
                    else
                    {
                        // 方向がない?
                        System.Diagnostics.Debug.Assert(false);
                    }
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
                        double kxValue = maPyy * sNxNx[row, col];
                        double kyValue = maPxx * sNyNy[row, col];
                        // 要素質量行列
                        double mValue = Constants.Ep0 * Constants.Mu0 * maQzz * sNN[row, col];

                        double aValue = 0;
                        // 共通
                        aValue += (1.0 / (dt * dt)) * mValue + NewmarkBeta * (kxValue + kyValue);
                        
                        if (maPML != null)
                        {
                            if (isXDirection && !isYDirection)
                            {
                                // X方向PML
                                aValue += (1.0 / (2.0 * dt)) * (sigmaX / (Constants.Ep0 * epxx)) * mValue;
                            }
                            else if (isYDirection && !isXDirection)
                            {
                                // Y方向PML
                                aValue += (1.0 / (2.0 * dt)) * (sigmaY / (Constants.Ep0 * epyy)) * mValue;
                            }
                            else if (isXDirection && isYDirection)
                            {
                                // XY方向PML
                                double sigmaEpX = sigmaX / epxx;
                                double sigmaEpY = sigmaY / epyy;
                                aValue +=
                                    (1.0 / (2.0 * dt)) * ((sigmaEpX + sigmaEpY) / Constants.Ep0) * mValue +
                                    NewmarkBeta * (sigmaEpX * sigmaEpY / (Constants.Ep0 * Constants.Ep0)) * mValue;
                            }
                            else
                            {
                                // 方向がない?
                                System.Diagnostics.Debug.Assert(false);
                            }

                        }

                        _A[rowNodeId, colNodeId] += aValue;
                    }
                }

                // 回転移動
                // 後片付け
                World.RotAngle = 0.0;
                World.RotOrigin = null;
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
            int refPortCnt = (int)World.GetPortCount(QuantityId) - 1; // 励振源を除く
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
            //--------------------------------------------------------------
            // 残差の計算
            //--------------------------------------------------------------
            B = new double[nodeCnt];
            CalcB(B);

            int refPortCnt = (int)World.GetPortCount(QuantityId) - 1;  // 励振源分引く

            //--------------------------------------------------------------
            // 励振源
            //--------------------------------------------------------------
            {
                int portId = refPortCnt; // ポートリストの最後の要素が励振境界
                var Qb = Qbs[portId];
                int nodeCntB = Qb.RowLength;
                double srcBetaX = SrcBetaXs[portId];
                double vpx = srcOmega / srcBetaX;
                //System.Numerics.Complex[] srcProfile = SrcProfiles[portId];
                System.Numerics.Complex[] srcModifyProfile = SrcModifyProfiles[portId];

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

                {
                    // dEzinc/dt
                    double[] srcUt = new double[nodeCntB];
                    for (int i = 0; i < nodeCntB; i++)
                    {
                        srcUt[i] = (2.0 / vpx) * (srcModifyProfile[i].Real * (srcU2 - srcU1) / (2.0 * dt) +
                            srcModifyProfile[i].Imaginary * (1.0 / srcOmega) * (1.0 / (dt * dt)) * 
                            (srcU2 - 2.0 * srcU0 + srcU1));
                    }
                    double[] vecQb = Qb * srcUt;
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);
                        B[nodeId] += vecQb[nodeIdB];
                    }
                }
            }
        }

        private void CalcB(double[] _B)
        {
            double dt = TimeDelta;

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
                System.Diagnostics.Debug.Assert(
                    ma0 is DielectricMaterial ||
                    ma0 is DielectricPMLMaterial);
                DielectricMaterial ma = null;
                DielectricPMLMaterial maPML = null;
                double epxx = 0;
                double epyy = 0;
                double epzz = 0;
                double muxx = 0;
                double muyy = 0;
                double muzz = 0;
                double rotAngle = 0.0;
                OpenTK.Vector2d rotOrigin = new OpenTK.Vector2d();
                if (ma0 is DielectricMaterial)
                {
                    ma = ma0 as DielectricMaterial;
                    epxx = ma.Epxx;
                    epyy = ma.Epyy;
                    epzz = ma.Epzz;
                    muxx = ma.Muxx;
                    muyy = ma.Muyy;
                    muzz = ma.Muzz;
                }
                else if (ma0 is DielectricPMLMaterial)
                {
                    maPML = ma0 as DielectricPMLMaterial;
                    epxx = maPML.Epxx;
                    epyy = maPML.Epyy;
                    epzz = maPML.Epzz;
                    muxx = maPML.Muxx;
                    muyy = maPML.Muyy;
                    muzz = maPML.Muzz;
                    rotOrigin = maPML.RotOriginPoint;
                    rotAngle = maPML.RotAngle;
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
                // 回転移動
                World.RotAngle = rotAngle;
                World.RotOrigin = new double[] { rotOrigin.X, rotOrigin.Y };
                double maPxx = 0;
                double maPyy = 0;
                double maQzz = 0;
                if (IsTMMode)
                {
                    // TMモード
                    maPxx = 1.0 / epxx;
                    maPyy = 1.0 / epyy;
                    maQzz = muzz;
                }
                else
                {
                    // TEモード
                    maPxx = 1.0 / muxx;
                    maPyy = 1.0 / muyy;
                    maQzz = epzz;
                }

                // 重心を求める
                OpenTK.Vector2d cPt;
                {
                    OpenTK.Vector2d[] vertexPts = new OpenTK.Vector2d[triFE.VertexCount];
                    for (int iVertex = 0; iVertex < triFE.VertexCount; iVertex++)
                    {
                        int coId = triFE.VertexCoordIds[iVertex];
                        double[] coord = World.GetCoord(QuantityId, coId);
                        vertexPts[iVertex] = new OpenTK.Vector2d(coord[0], coord[1]);
                    }
                    cPt = (vertexPts[0] + vertexPts[1] + vertexPts[2]) / 3.0;
                }

                // PML
                bool isXDirection = false;
                bool isYDirection = false;
                // X方向PML
                double sigmaX = 0.0;
                double c1PX = 0.0;
                double c2PX = 0.0;
                double c1VX = 0.0;
                // Y方向PML
                double sigmaY = 0.0;
                double c1PY = 0.0;
                double c2PY = 0.0;
                double c1VY = 0.0;
                // ψx,y
                double[] px = null;
                double[] py = null;
                // vx,y
                double[] vx = null;
                double[] vy = null;
                if (maPML != null)
                {
                    isXDirection = maPML.IsXDirection();
                    isYDirection = maPML.IsYDirection();
                    if (isXDirection && !isYDirection)
                    {
                        maPML.CalcSigmaXForTD(cPt, dt, out sigmaX, out c1PX, out c2PX, out c1VX);
                    }
                    else if (isYDirection && !isXDirection)
                    {
                        maPML.CalcSigmaYForTD(cPt, dt, out sigmaY, out c1PY, out c2PY, out c1VY);
                    }
                    else if (isXDirection && isYDirection)
                    {
                        maPML.CalcSigmaXForTD(cPt, dt, out sigmaX, out c1PX, out c2PX, out c1VX);
                        maPML.CalcSigmaYForTD(cPt, dt, out sigmaY, out c1PY, out c2PY, out c1VY);
                    }
                    else
                    {
                        // 方向がない?
                        System.Diagnostics.Debug.Assert(false);
                    }

                    // ψx,y
                    px = Pxs[feId - 1];
                    if (px == null)
                    {
                        px = new double[elemNodeCnt];
                        Pxs[feId - 1] = px;
                    }
                    py = Pys[feId - 1];
                    if (py == null)
                    {
                        py = new double[elemNodeCnt];
                        Pys[feId - 1] = py;
                    }
                    // vx,y
                    vx = Vxs[feId - 1];
                    if (vx == null)
                    {
                        vx = new double[elemNodeCnt];
                        Vxs[feId - 1] = vx;
                    }
                    vy = Vys[feId - 1];
                    if (vy == null)
                    {
                        vy = new double[elemNodeCnt];
                        Vys[feId - 1] = vy;
                    }

                    // ψx,y,vx,yの更新
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int nodeId = nodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }
                        px[iNode] = c2PX * px[iNode] + c1PX * EzPrev[nodeId];
                        py[iNode] = c2PY * py[iNode] + c1PY * EzPrev[nodeId];
                        vx[iNode] = vx[iNode] + c1VX * EzPrev[nodeId];
                        vy[iNode] = vy[iNode] + c1VY * EzPrev[nodeId];
                    }
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
                        double kxValue = maPyy * sNxNx[row, col];
                        double kyValue = maPxx * sNyNy[row, col];
                        // 要素質量行列
                        double mValue = Constants.Ep0 * Constants.Mu0 * maQzz * sNN[row, col];

                        double bValue = 0;
                        // 共通
                        bValue += 
                            ((2.0 / (dt * dt)) * EzPrev[colNodeId] -
                            (1.0 / (dt * dt)) * EzPrev2[colNodeId]) * mValue;
                        bValue += 
                            (-(1.0 - 2.0 * NewmarkBeta) * EzPrev[colNodeId] -
                            NewmarkBeta * EzPrev2[colNodeId]) * 
                            (kxValue + kyValue);

                        if (maPML != null)
                        {
                            if (isXDirection && !isYDirection)
                            {
                                // X方向PML
                                // Mσ
                                bValue +=
                                    (1.0 / (2.0 * dt)) * (sigmaX / (Constants.Ep0 * epxx)) *
                                    EzPrev2[colNodeId] * mValue;

                                // ψx
                                bValue += kxValue * px[col];

                                // vx
                                bValue += -1.0 * kyValue * vx[col];

                            }
                            else if (isYDirection && !isXDirection)
                            {
                                // Y方向PML
                                // Mσ
                                bValue +=
                                    (1.0 / (2.0 * dt)) * (sigmaY / (Constants.Ep0 * epyy)) *
                                    EzPrev2[colNodeId] * mValue;

                                // ψy
                                bValue += kyValue * py[col];

                                // vy
                                bValue += -1.0 * kxValue * vy[col];
                            }
                            else if (isXDirection && isYDirection)
                            {
                                // XY方向PML

                                double sigmaEpX = sigmaX / epxx;
                                double sigmaEpY = sigmaY / epyy;
                                // Mσ1、Mσ2
                                bValue +=
                                    (1.0 / (2.0 * dt)) * ((sigmaEpX + sigmaEpY) / Constants.Ep0) *
                                    EzPrev2[colNodeId] * mValue +
                                    -1.0 * ((1.0 - 2.0 * NewmarkBeta) * EzPrev[colNodeId] +
                                    NewmarkBeta * EzPrev2[colNodeId]) *
                                    (sigmaEpX * sigmaEpY / (Constants.Ep0 * Constants.Ep0)) * mValue;

                                // ψx
                                bValue += ((sigmaEpX - sigmaEpY) / sigmaEpX) * kxValue * px[col];

                                // ψy
                                bValue += ((sigmaEpY - sigmaEpX) / sigmaEpY) * kyValue * py[col];
                            }
                            else
                            {
                                // 方向がない?
                                System.Diagnostics.Debug.Assert(false);
                            }
                        }

                        _B[rowNodeId] += bValue;
                    }
                }

                // 回転移動
                // 後片付け
                World.RotAngle = 0.0;
                World.RotOrigin = null;
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

            int refPortCnt = (int)World.GetPortCount(QuantityId) - 1; // 励振面を引く

            if (RefVIds.Count > 0)
            {
                CalcSParameterAtRefPoints(freqDomainAmpsInc, out freqs, out freqDomainAmpss, out Sss);
            }
            else if (refPortCnt > 0)
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
            int refPortCnt = (int)World.GetPortCount(QuantityId) - 1; // 励振面を引く
            System.Diagnostics.Debug.Assert(refPortCnt > 0);
            double[][][] datasRefs = new double[refPortCnt][][];
            for (int refIndex = 0; refIndex < refPortCnt; refIndex++)
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

            for (int refIndex = 0; refIndex < refPortCnt; refIndex++)
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

                    int portId = refIndex;
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
            for (int refIndex = 0; refIndex < refPortCnt; refIndex++)
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
