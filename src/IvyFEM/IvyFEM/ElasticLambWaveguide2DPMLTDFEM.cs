using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class ElasticLambWaveguide2DPMLTDFEM : FEM
    {
        //---------------------------------
        // ※ポートの設定順序
        // 参照(観測)面、励振源の順
        //---------------------------------

        public uint QuantityId { get; private set; } = 0;

        public double TimeStep { get; set; } = 0.0;
        public double NewmarkBeta { get; private set; } = 1.0 / 4.0;
        public double NewmarkGamma { get; private set; } = 1.0 / 2.0;
        public uint UValueId { get; private set; } = 0;
        public uint PrevUValueId { get; private set; } = 0;

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
        /// 観測点の頂点ID
        /// </summary>
        public IList<uint> RefVIds { get; set; } = new List<uint>();

        /// <summary>
        /// [A]
        /// </summary>
        private IvyFEM.Linear.DoubleSparseMatrix A = null;
        private IvyFEM.Linear.DoubleSparseMatrix A0 = null; // 固定境界適用前
        /// <summary>
        /// {b}
        /// </summary>
        private double[] B = null;
        
        /// <summary>
        /// 境界行列リスト(ポート単位) sNN
        /// </summary>
        private IList<IvyFEM.Lapack.ComplexMatrix> SNNs = null;
        /// <summary>
        /// 境界行列リスト(ポート単位) sNNy
        /// </summary>
        private IList<IvyFEM.Lapack.ComplexMatrix> SNNys = null;

        /// <summary>
        /// 境界の界の伝搬定数(ポート単位)
        /// </summary>
        public IList<System.Numerics.Complex> SrcBetaXs { get; private set; } = null;
        /// <summary>
        /// 境界の界のモード分布(ポート単位)
        /// </summary>
        public IList<System.Numerics.Complex[]> SrcHFs { get; private set; } = null;
        public IList<System.Numerics.Complex[]> SrcHGs { get; private set; } = null;
        public ElasticLambWaveguide1DEigenFEM[] EigenFEMs { get; private set; }

        /// <summary>
        /// 変位
        /// </summary>
        public double[] U { get; private set; } = null;
        public double[] CoordU { get; private set; }
        public double[] CoordSigma { get; private set; }

        /// <summary>
        /// ψdx(1つ前→現在値に変化する) Pdxs[要素番号][要素内節点番号 * 2 + d]
        /// </summary>
        private double[][] Pdxs = null;
        private double[][] Pdys = null;
        /// <summary>
        /// wdx(1つ前→現在値に変化する) Wdxs[要素番号][要素内節点番号 * 2 + d]
        /// </summary>
        private double[][] Wdxs = null;
        private double[][] Wdys = null;

        /// <summary>
        /// 観測点の変位(時間変化リスト)
        /// </summary>
        public IList<IList<double[]>> RefTimeUsss { get; private set; } = new List<IList<double[]>>();
        public IList<IList<double[]>> RefTimeSigmasss { get; private set; } = new List<IList<double[]>>();

        public ElasticLambWaveguide2DPMLTDFEM(
            FEWorld world,
            double timeStep,
            double newmarkBeta, double newmarkGamma,
            uint uValueId, uint prevUValueId)
        {
            World = world;
            TimeStep = timeStep;
            NewmarkBeta = newmarkBeta;
            NewmarkGamma = newmarkGamma;
            UValueId = uValueId;
            PrevUValueId = prevUValueId;
        }

        public void UpdateFieldValuesTimeDomain()
        {
            UpdateFieldValuesNewmarkBetaTimeDomain(
                U, UValueId, PrevUValueId,
                TimeStep,
                NewmarkBeta, NewmarkGamma);
        }

        public override void Solve()
        {
            int t;
            int uDof = 2;
            int sDof = 2;
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
                RefTimeUsss.Clear();
                RefTimeSigmasss.Clear();
                for (int i = 0; i < cnt; i++)
                {
                    var refTimeUss = new List<double[]>();
                    var refTimeSigmass = new List<double[]>();
                    RefTimeUsss.Add(refTimeUss);
                    RefTimeSigmasss.Add(refTimeSigmass);
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

            //--------------------------------------------------------------
            // 残差
            //--------------------------------------------------------------
            //t = System.Environment.TickCount;
            CalcB();
            //System.Diagnostics.Debug.WriteLine("CalcB t = " + (System.Environment.TickCount - t));

            //--------------------------------------------------------------
            // 固定境界条件
            //--------------------------------------------------------------
            A = new IvyFEM.Linear.DoubleSparseMatrix(A0); // 固定境界適用前のAをコピーする
            DoubleSetFixedCadsCondtion(A, B);

            //------------------------------------------------------------------
            // Ux,Uyを求める
            //------------------------------------------------------------------
            {
                double[] X;
                Solver.DoubleSolve(out X, A, B);
                U = X;
            }

            double[] nodeUonly = U;
            // uを抽出(座標ベース)
            double[] coordUonly;
            // σを算出(座標ベース)
            double[] coordSigma;
            CalcCoordUCoordSigma(World,nodeUonly, out coordUonly, out coordSigma);
            CoordU = coordUonly;
            CoordSigma = coordSigma;

            if (RefVIds.Count > 0)
            {
                System.Diagnostics.Debug.Assert(false);
                throw new NotImplementedException();
            }
            else if (refPortCnt > 0)
            {
                for (int refIndex = 0; refIndex < refPortCnt; refIndex++)
                {
                    int portId = refIndex;
                    int nodeCntB = (int)World.GetPortNodeCount(QuantityId, (uint)portId);
                    double[] uValues = new double[nodeCntB * uDof];
                    double[] sigmaValues = new double[nodeCntB * sDof];
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            //uValues[nodeIdB * uDof + iDof] = U[nodeId * uDof + iDof];
                            uValues[nodeIdB * uDof + iDof] = CoordU[coId * uDof + iDof];
                            sigmaValues[nodeIdB * sDof + iDof] = CoordSigma[coId * sDof + iDof];
                        }
                    }
                    var refTimeUss = RefTimeUsss[refIndex];
                    var refTimeSigmass = RefTimeSigmasss[refIndex];
                    refTimeUss.Add(uValues);
                    refTimeSigmass.Add(sigmaValues);
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

            int uDof = 2;
            int sDof = 2;

            //------------------------------------------------------
            // 全体係数行列の作成
            //------------------------------------------------------
            int nodeCnt = (int)World.GetNodeCount(QuantityId);
            A = new IvyFEM.Linear.DoubleSparseMatrix(nodeCnt * uDof, nodeCnt * uDof);
            _CalcA();

            //--------------------------------------------------------------
            // 固定境界条件を適用するとAが変更されるので
            //--------------------------------------------------------------
            A0 = new IvyFEM.Linear.DoubleSparseMatrix(A); // 退避

            // バンド幅縮小を座標を用いて行う
            SetupCoordsForBandMatrix();

            //////////////////////////////////////////////////////
            int portCnt = (int)World.GetPortCount(QuantityId) - 1; // 励振源を除く
            // 周波数
            double srcFreq = SrcFrequency;
            SNNs = new List<IvyFEM.Lapack.ComplexMatrix>();
            SNNys = new List<IvyFEM.Lapack.ComplexMatrix>();
            SrcBetaXs = new List<System.Numerics.Complex>();
            SrcHFs = new List<System.Numerics.Complex[]>();
            SrcHGs = new List<System.Numerics.Complex[]>();

            EigenFEMs = new ElasticLambWaveguide1DEigenFEM[(portCnt + 1)];
            var portConditions = World.GetPortConditions(QuantityId);
            for (uint portId = 0; portId < (portCnt + 1); portId++)
            {
                uint nodeCntB = World.GetPortNodeCount(QuantityId, portId);
                double normalX = 1.0;
                {
                    var portConditon = portConditions[(int)portId];
                    normalX = portConditon.GetDoubleAdditionalParameters()[0];
                    System.Diagnostics.Debug.Assert(Math.Abs(Math.Abs(normalX) - 1.0) < Constants.PrecisionLowerLimit);
                }

                //------------------------------------------------------
                // モード分布計算
                //------------------------------------------------------
                var eigenFEM = new ElasticLambWaveguide1DEigenFEM(World, QuantityId, portId);
                EigenFEMs[portId] = eigenFEM;

                eigenFEM.Frequency = srcFreq;
                eigenFEM.NormalX = normalX;
                eigenFEM.Solve();
                System.Numerics.Complex[] betas = eigenFEM.Betas;
                System.Numerics.Complex[][] hUEVecs = eigenFEM.HUEVecs;
                System.Numerics.Complex[][] hSigmaEVecs = eigenFEM.HSigmaEVecs;

                System.Diagnostics.Debug.WriteLine("port = {0} mode Count = {1}", portId, betas.Length);
                // 基本モード
                int iMode = 0;
                System.Numerics.Complex beta = betas[iMode];
                SrcBetaXs.Add(beta);

                var srcHF = new System.Numerics.Complex[nodeCntB * uDof];
                var srcHG = new System.Numerics.Complex[nodeCntB * sDof];
                SrcHFs.Add(srcHF);
                SrcHGs.Add(srcHG);

                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                {
                    srcHF[nodeIdB * uDof] = hUEVecs[iMode][nodeIdB * uDof];
                    srcHF[nodeIdB * uDof + 1] = hUEVecs[iMode][nodeIdB * uDof + 1];
                    srcHG[nodeIdB * sDof] = hSigmaEVecs[iMode][nodeIdB * sDof];
                    srcHG[nodeIdB * sDof + 1] = hSigmaEVecs[iMode][nodeIdB * sDof + 1];
                }

                //------------------------------------------------
                IvyFEM.Lapack.ComplexMatrix sNN;
                IvyFEM.Lapack.ComplexMatrix sNNy;
                eigenFEM.CalcBoundaryMatrixForPML(out sNN, out sNNy);
                SNNs.Add(sNN);
                SNNys.Add(sNNy);
            }

            int feCnt = World.GetTriangleFEIds(QuantityId).Count;
            Pdxs = new double[feCnt][];
            Pdys = new double[feCnt][];
            Wdxs = new double[feCnt][];
            Wdys = new double[feCnt][];
        }

        // バンド幅縮小を座標を用いて行う
        private void SetupCoordsForBandMatrix()
        {
            int uDof = 2;
            int nodeCnt = (int)World.GetNodeCount(QuantityId);
            if (Solver is IvyFEM.Linear.LapackEquationSolver)
            {
                var solver = Solver as IvyFEM.Linear.LapackEquationSolver;
                if (solver.IsOrderingToBandMatrix)
                {
                    int[] coIds = new int[nodeCnt * uDof];
                    for (int nodeId = 0; nodeId < nodeCnt; nodeId++)
                    {
                        int coId = World.Node2Coord(QuantityId, nodeId);
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            coIds[nodeId * uDof + iDof] = coId;
                        }
                    }
                    double[][] coords = new double[coIds.Length][];
                    for (int rowIndex = 0; rowIndex < coIds.Length; rowIndex++)
                    {
                        int coId = coIds[rowIndex];
                        double[] coord = World.GetCoord(QuantityId, coId);
                        coords[rowIndex] = coord;
                    }

                    solver.CoordsForBandMatrix = coords;
                }
            }
        }

        private void _CalcA()
        {
            int uDof = 2;
            double dt = TimeStep;
            double beta = NewmarkBeta;
            double gamma = NewmarkGamma;

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
                    ma0 is LinearElasticMaterial ||
                    ma0 is LinearElasticPMLMaterial);
                LinearElasticMaterial ma = null;
                LinearElasticPMLMaterial maPML = null;
                double rho = 0;
                double lambda = 0;
                double mu = 0;
                double rotAngle = 0.0;
                OpenTK.Vector2d rotOrigin = new OpenTK.Vector2d();
                if (ma0 is LinearElasticMaterial)
                {
                    ma = ma0 as LinearElasticMaterial;
                    rho = ma.MassDensity;
                    lambda = ma.LameLambda;
                    mu = ma.LameMu;
                }
                else if (ma0 is LinearElasticPMLMaterial)
                {
                    maPML = ma0 as LinearElasticPMLMaterial;
                    rho = maPML.MassDensity;
                    lambda = maPML.LameLambda;
                    mu = maPML.LameMu;
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
                double AX = 1.0;
                double DX = 0.0;
                double c1PX = 0.0;
                double c2PX = 0.0;
                double c1WX = 0.0;
                // Y方向PML
                double AY = 1.0;
                double DY = 0.0;
                double c1PY = 0.0;
                double c2PY = 0.0;
                double c1WY = 0.0;
                if (maPML != null)
                {
                    isXDirection = maPML.IsXDirection();
                    isYDirection = maPML.IsYDirection();
                    if (isXDirection && !isYDirection)
                    {
                        maPML.CalcScalingDampingFactorXForTD(cPt, dt,
                            out AX, out DX, out c1PX, out c2PX, out c1WX);
                    }
                    else if (isYDirection && !isXDirection)
                    {
                        maPML.CalcScalingDampingFactorYForTD(cPt, dt,
                            out AY, out DY, out c1PY, out c2PY, out c1WY);
                    }
                    else if (isXDirection && isYDirection)
                    {
                        maPML.CalcScalingDampingFactorXForTD(cPt, dt,
                            out AX, out DX, out c1PX, out c2PX, out c1WX);
                        maPML.CalcScalingDampingFactorYForTD(cPt, dt,
                            out AY, out DY, out c1PY, out c2PY, out c1WY);
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
                double[,] sNxNy = sNuNv[0, 1];
                double[,] sNyNx = sNuNv[1, 0];
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

                        // K1
                        double kxValue1 = (lambda + 2.0 * mu) * (1.0 / AX) * sNxNx[row, col];
                        double kxyValue1 = lambda * sNxNy[row, col] + mu * sNyNx[row, col];
                        double kyValue1 = lambda * AX * sNyNy[row, col];
                        // M11
                        double mValue11 = rho * AX * sNN[row, col];
                        // M12
                        double mValue12 = rho * AX * DX * sNN[row, col];
                        // K2
                        double kyValue2 = (lambda + 2.0 * mu) * AX * sNyNy[row, col];
                        double kxyValue2 = lambda * sNyNx[row, col] + mu * sNxNy[row, col];
                        double kxValue2 = mu * (1.0 / AX) * sNxNx[row, col];
                        // M21
                        double mValue21 = rho * AX * sNN[row, col];
                        // M22
                        double mValue22 = rho * AX * DX * sNN[row, col];

                        double axxValue = 0.0;
                        double axyValue = 0.0;
                        double ayxValue = 0.0;
                        double ayyValue = 0.0;
                        //!!!!!!!!!!!!!
                        // X方向PML
                        System.Diagnostics.Debug.Assert(maPML == null ||
                            (maPML != null && isXDirection && !isYDirection));
                        //
                        axxValue +=
                            (1.0 / (beta * dt * dt)) * mValue11 +
                            (gamma / (beta * dt)) * mValue12 +
                            kxValue1 + kyValue1;
                        axyValue += kxyValue1;
                        //
                        ayxValue += kxyValue2;
                        ayyValue +=
                            (1.0 / (beta * dt * dt)) * mValue21 +
                            (gamma / (beta * dt)) * mValue22 +
                            kyValue2 + kxValue2;

                        A[rowNodeId * uDof, colNodeId * uDof] += axxValue;
                        A[rowNodeId * uDof, colNodeId * uDof + 1] += axyValue;
                        A[rowNodeId * uDof + 1, colNodeId * uDof] += ayxValue;
                        A[rowNodeId * uDof + 1, colNodeId * uDof + 1] += ayyValue;
                    }
                }

                // 回転移動
                // 後片付け
                World.RotAngle = 0.0;
                World.RotOrigin = null;
            }
        }

        private void CalcB()
        {
            int uDof = 2;
            int sDof = 2;
            double dt = TimeStep;

            // 周波数
            double srcFreq = SrcFrequency;
            // 角周波数
            double srcOmega = 2.0 * Math.PI * srcFreq;

            //--------------------------------------------------------------
            // 残差の計算
            //--------------------------------------------------------------
            B = new double[A.RowLength];
            // 残差
            _CalcB();

            int refPortCnt = (int)World.GetPortCount(QuantityId) - 1;  // 励振源分引く
            //--------------------------------------------------------------
            // 励振源
            //--------------------------------------------------------------
            {
                int portId = refPortCnt; // ポートリストの最後の要素が励振境界
                var sNN = SNNs[portId];
                var sNNReal = sNN.Real();
                int nodeCntB = sNN.RowLength;
                System.Numerics.Complex srcBetaX = SrcBetaXs[portId];
                System.Numerics.Complex[] srcHF = SrcHFs[portId];
                System.Numerics.Complex[] srcHG = SrcHGs[portId];
                double rho = 0;
                double lambda = 0;
                double mu = 0;
                {
                    var feId = World.GetPortLineFEIds(QuantityId, (uint)portId)[0];
                    var lineFE = World.GetLineFE(QuantityId, feId);
                    Material ma0 = World.GetMaterial(lineFE.MaterialId);
                    System.Diagnostics.Debug.Assert(ma0 is LinearElasticMaterial);
                    OpenTK.Vector2d rotOrigin = new OpenTK.Vector2d();
                    LinearElasticMaterial ma = ma0 as LinearElasticMaterial;
                    rho = ma.MassDensity;
                    lambda = ma.LameLambda;
                    mu = ma.LameMu;
                }
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
                    double[] srcSigmaXX = new double[nodeCntB];
                    double[] srcSigmaYX = new double[nodeCntB];
                    System.Diagnostics.Debug.Assert(srcHG.Length == (nodeCntB * uDof));
                    for (int i = 0; i < nodeCntB; i++)
                    {
                        System.Numerics.Complex gXX =
                            -1.0 * System.Numerics.Complex.ImaginaryOne * srcHG[i * sDof];
                        System.Numerics.Complex gYX = srcHG[i * sDof + 1];
                        //----------
                        // 位相をσxx基準にする
                        gXX *= System.Numerics.Complex.ImaginaryOne;
                        gYX *= System.Numerics.Complex.ImaginaryOne;
                        //----------
                        double amp = GaussianAmp;

                        srcSigmaXX[i] = gXX.Real * amp * srcU0;
                        srcSigmaYX[i] = gYX.Real * amp * srcU0;
                    }
                    double[] vecX = sNNReal * srcSigmaXX;
                    double[] vecY = sNNReal * srcSigmaYX;
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);
                        B[nodeId * uDof] += vecX[nodeIdB];
                        B[nodeId * uDof + 1] += vecY[nodeIdB];
                    }
                }
                /*
                // viscous wave
                {
                    // dux/dt
                    double[] srcUxt = new double[nodeCntB];
                    double[] srcUyt = new double[nodeCntB];
                    System.Diagnostics.Debug.Assert(srcHF.Length == (nodeCntB * uDof));
                    for (int i = 0; i < nodeCntB; i++)
                    {
                        System.Numerics.Complex ux = srcHF[i * uDof];
                        System.Numerics.Complex uy =
                            System.Numerics.Complex.ImaginaryOne * srcHF[i * uDof + 1];
                        double amp = GaussianAmp;
                        //----------
                        // 位相をux基準にする
                        //ux *= 1.0;
                        //uy *= 1.0;
                        //----------

                        srcUxt[i] = ux.Real * amp * (srcU2 - srcU1) / (2.0 * dt);
                        srcUyt[i] = uy.Real * amp * (srcU2 - srcU1) / (2.0 * dt);
                    }
                    double vp = Math.Sqrt((lambda + 2.0 * mu) / rho);
                    double vs = Math.Sqrt(mu / rho);
                    double[] workX = IvyFEM.Lapack.Functions.dscal(srcUxt, (-1.0 * rho * vp));
                    double[] workY = IvyFEM.Lapack.Functions.dscal(srcUyt, (-1.0 * rho * vs));
                    double[] vecX = sNNReal * workX;
                    double[] vecY = sNNReal * workY;
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);
                        B[nodeId * uDof] += vecX[nodeIdB];
                        B[nodeId * uDof + 1] += vecY[nodeIdB];
                    }
                }
                */
            }
        }

        private void _CalcB()
        {
            int uDof = 2;
            double dt = TimeStep;
            double beta = NewmarkBeta;
            double gamma = NewmarkGamma;
            var uFV = World.GetFieldValue(UValueId);

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
                    ma0 is LinearElasticMaterial ||
                    ma0 is LinearElasticPMLMaterial);
                LinearElasticMaterial ma = null;
                LinearElasticPMLMaterial maPML = null;
                double rho = 0;
                double lambda = 0;
                double mu = 0;
                double rotAngle = 0.0;
                OpenTK.Vector2d rotOrigin = new OpenTK.Vector2d();
                if (ma0 is LinearElasticMaterial)
                {
                    ma = ma0 as LinearElasticMaterial;
                    rho = ma.MassDensity;
                    lambda = ma.LameLambda;
                    mu = ma.LameMu;
                }
                else if (ma0 is LinearElasticPMLMaterial)
                {
                    maPML = ma0 as LinearElasticPMLMaterial;
                    rho = maPML.MassDensity;
                    lambda = maPML.LameLambda;
                    mu = maPML.LameMu;
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
                double AX = 1.0;
                double DX = 0.0;
                double c1PX = 0.0;
                double c2PX = 0.0;
                double c1WX = 0.0;
                // Y方向PML
                double AY = 1.0;
                double DY = 0.0;
                double c1PY = 0.0;
                double c2PY = 0.0;
                double c1WY = 0.0;
                // ψdx,ψdy
                double[] pdx = null;
                double[] pdy = null;
                // wdx,wdy
                double[] wdx = null;
                double[] wdy = null;
                if (maPML != null)
                {
                    isXDirection = maPML.IsXDirection();
                    isYDirection = maPML.IsYDirection();
                    if (isXDirection && !isYDirection)
                    {
                        maPML.CalcScalingDampingFactorXForTD(cPt, dt,
                            out AX, out DX, out c1PX, out c2PX, out c1WX);
                    }
                    else if (isYDirection && !isXDirection)
                    {
                        maPML.CalcScalingDampingFactorYForTD(cPt, dt,
                            out AY, out DY, out c1PY, out c2PY, out c1WY);
                    }
                    else if (isXDirection && isYDirection)
                    {
                        maPML.CalcScalingDampingFactorXForTD(cPt, dt,
                            out AX, out DX, out c1PX, out c2PX, out c1WX);
                        maPML.CalcScalingDampingFactorYForTD(cPt, dt,
                            out AY, out DY, out c1PY, out c2PY, out c1WY);
                    }
                    else
                    {
                        // 方向がない?
                        System.Diagnostics.Debug.Assert(false);
                    }

                    // ψdx,ψdy
                    pdx = Pdxs[feId - 1];
                    if (pdx == null)
                    {
                        pdx = new double[elemNodeCnt * uDof];
                        Pdxs[feId - 1] = pdx;
                    }
                    pdy = Pdys[feId - 1];
                    if (pdy == null)
                    {
                        pdy = new double[elemNodeCnt * uDof];
                        Pdys[feId - 1] = pdy;
                    }
                    // wdx,wdy
                    wdx = Wdxs[feId - 1];
                    if (wdx == null)
                    {
                        wdx = new double[elemNodeCnt * uDof];
                        Wdxs[feId - 1] = wdx;
                    }
                    wdy = Wdys[feId - 1];
                    if (wdy == null)
                    {
                        wdy = new double[elemNodeCnt * uDof];
                        Wdys[feId - 1] = wdy;
                    }

                    // ψdx,ψdy,wdx,wdyの更新
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int nodeId = nodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }
                        int iCoId = triFE.NodeCoordIds[iNode];
                        double[] u = uFV.GetDoubleValue(iCoId, FieldDerivativeType.Value);
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            pdx[iNode * uDof + iDof] =
                                c2PX * pdx[iNode * uDof + iDof] + c1PX * u[iDof];
                            pdy[iNode * uDof + iDof] =
                                c2PY * pdy[iNode * uDof + iDof] + c1PY * u[iDof];
                            wdx[iNode * uDof + iDof] =
                                wdx[iNode * uDof + iDof] + c1WX * u[iDof];
                            wdy[iNode * uDof + iDof] =
                                wdy[iNode * uDof + iDof] + c1WY * u[iDof];
                        }
                    }
                }

                double[,] sNN = triFE.CalcSNN();
                double[,][,] sNuNv = triFE.CalcSNuNv();
                double[,] sNxNx = sNuNv[0, 0];
                double[,] sNxNy = sNuNv[0, 1];
                double[,] sNyNx = sNuNv[1, 0];
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

                        int colCoId = triFE.NodeCoordIds[col];
                        double[] u = uFV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                        double[] velU = uFV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                        double[] accU = uFV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);
                        System.Diagnostics.Debug.Assert(u.Length == uDof);
                        System.Diagnostics.Debug.Assert(velU.Length == uDof);
                        System.Diagnostics.Debug.Assert(accU.Length == uDof);

                        // K1
                        double kxValue1 = (lambda + 2.0 * mu) * (1.0 / AX) * sNxNx[row, col];
                        double kxyValue1 = lambda * sNxNy[row, col] + mu * sNyNx[row, col];
                        double kyValue1 = lambda * AX * sNyNy[row, col];
                        // M11
                        double mValue11 = rho * AX * sNN[row, col];
                        // M12
                        double mValue12 = rho * AX * DX * sNN[row, col];
                        // K2
                        double kyValue2 = (lambda + 2.0 * mu) * AX * sNyNy[row, col];
                        double kxyValue2 = lambda * sNyNx[row, col] + mu * sNxNy[row, col];
                        double kxValue2 = mu * (1.0 / AX) * sNxNx[row, col];
                        // M21
                        double mValue21 = rho * AX * sNN[row, col];
                        // M22
                        double mValue22 = rho * AX * DX * sNN[row, col];

                        double bxValue = 0.0;
                        double byValue = 0.0;
                        //!!!!!!!!!!!!!
                        // X方向PML
                        System.Diagnostics.Debug.Assert(maPML == null ||
                            (maPML != null && isXDirection && !isYDirection));
                        //
                        bxValue +=
                            mValue11 * (
                            (1.0 / (beta * dt * dt)) * u[0] +
                            (1.0 / (beta * dt)) * velU[0] +
                            ((1.0 / (2.0 * beta)) - 1.0) * accU[0]) +
                            mValue12 * (
                            (gamma / (beta * dt)) * u[0] +
                            ((gamma / beta) - 1.0) * velU[0] +
                            ((gamma / (2.0 * beta)) - 1.0) * dt * accU[0]
                            );
                        //
                        byValue +=
                            mValue21 * (
                            (1.0 / (beta * dt * dt)) * u[1] +
                            (1.0 / (beta * dt)) * velU[1] +
                            ((1.0 / (2.0 * beta)) - 1.0) * accU[1]) +
                            mValue22 * (
                            (gamma / (beta * dt)) * u[1] +
                            ((gamma / beta) - 1.0) * velU[1] +
                            ((gamma / (2.0 * beta)) - 1.0) * dt * accU[1]
                            );

                        if (maPML != null)
                        {
                            if (isXDirection && !isYDirection)
                            {
                                // X方向PML

                                // ψxx
                                bxValue += kxValue1 * pdx[col * uDof];

                                // wxx
                                bxValue += -1.0 * kyValue1 * wdx[col * uDof];

                                // ψyx
                                byValue += kxValue2 * pdx[col * uDof + 1];

                                // wyx
                                byValue += -1.0 * kyValue2 * wdx[col * uDof + 1];

                            }
                            else if (isYDirection && !isXDirection)
                            {
                                // Y方向PML

                                System.Diagnostics.Debug.Assert(false);
                                throw new NotImplementedException();
                            }
                            else if (isXDirection && isYDirection)
                            {
                                // XY方向PML

                                System.Diagnostics.Debug.Assert(false);
                                throw new NotImplementedException();
                            }
                            else
                            {
                                // 方向がない?
                                System.Diagnostics.Debug.Assert(false);
                            }
                        }

                        B[rowNodeId * uDof] += bxValue;
                        B[rowNodeId * uDof + 1] += byValue;
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
                System.Diagnostics.Debug.Assert(false);
                throw new NotImplementedException();
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

        private void CalcSParameterAtRefPorts(
            System.Numerics.Complex[] freqDomainAmpsInc,
            out double[] freqs,
            out IList<System.Numerics.Complex[]> freqDomainAmpss,
            out IList<System.Numerics.Complex[]> Sss)
        {
            int uDof = 2;
            int sDof = 2;
            int refPortCnt = (int)World.GetPortCount(QuantityId) - 1; // 励振面を引く
            System.Diagnostics.Debug.Assert(refPortCnt > 0);
            double[][][] datasRefUxs = new double[refPortCnt][][];
            double[][][] datasRefUys = new double[refPortCnt][][];
            double[][][] datasRefSigmaXXs = new double[refPortCnt][][];
            double[][][] datasRefSigmaYXs = new double[refPortCnt][][];
            for (int refIndex = 0; refIndex < refPortCnt; refIndex++)
            {
                uint portId = (uint)refIndex;
                int timeCnt = RefTimeUsss[refIndex].Count;
                int nodeCntB = (int)World.GetPortNodeCount(QuantityId, portId);
                System.Diagnostics.Debug.Assert(nodeCntB == RefTimeUsss[refIndex][0].Length / uDof);
                double[][] datasRefUx = new double[nodeCntB][];
                double[][] datasRefUy = new double[nodeCntB][];
                double[][] datasRefSigmaXX = new double[nodeCntB][];
                double[][] datasRefSigmaYX = new double[nodeCntB][];
                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                {
                    datasRefUx[nodeIdB] = new double[timeCnt];
                    datasRefUy[nodeIdB] = new double[timeCnt];
                    datasRefSigmaXX[nodeIdB] = new double[timeCnt];
                    datasRefSigmaYX[nodeIdB] = new double[timeCnt];
                    for (int timeIndex = 0; timeIndex < timeCnt; timeIndex++)
                    {
                        datasRefUx[nodeIdB][timeIndex] = RefTimeUsss[refIndex][timeIndex][nodeIdB * uDof];
                        datasRefUy[nodeIdB][timeIndex] = RefTimeUsss[refIndex][timeIndex][nodeIdB * uDof + 1];
                        datasRefSigmaXX[nodeIdB][timeIndex] = RefTimeSigmasss[refIndex][timeIndex][nodeIdB * sDof];
                        datasRefSigmaYX[nodeIdB][timeIndex] = RefTimeSigmasss[refIndex][timeIndex][nodeIdB * sDof + 1];
                    }
                }
                datasRefUxs[refIndex] = datasRefUx;
                datasRefUys[refIndex] = datasRefUy;
                datasRefSigmaXXs[refIndex] = datasRefSigmaXX;
                datasRefSigmaYXs[refIndex] = datasRefSigmaYX;
            }

            int dataCnt = datasRefUxs[0][0].Length;
            System.Diagnostics.Debug.Assert(datasRefUys[0][0].Length == dataCnt);
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

            double[][][][] datasRefUSigmass = { datasRefUxs, datasRefUys, datasRefSigmaXXs, datasRefSigmaYXs };
            for (int refIndex = 0; refIndex < refPortCnt; refIndex++)
            {
                int nodeCntB = datasRefUxs[refIndex].Length;
                System.Numerics.Complex[][] freqDomainUxDatass = new System.Numerics.Complex[dataCnt][];
                System.Numerics.Complex[][] freqDomainUyDatass = new System.Numerics.Complex[dataCnt][];
                System.Numerics.Complex[][] freqDomainSigmaXXDatass = new System.Numerics.Complex[dataCnt][];
                System.Numerics.Complex[][] freqDomainSigmaYXDatass = new System.Numerics.Complex[dataCnt][];
                for (int freqIndex = 0; freqIndex < dataCnt; freqIndex++)
                {
                    freqDomainUxDatass[freqIndex] = new System.Numerics.Complex[nodeCntB];
                    freqDomainUyDatass[freqIndex] = new System.Numerics.Complex[nodeCntB];
                    freqDomainSigmaXXDatass[freqIndex] = new System.Numerics.Complex[nodeCntB];
                    freqDomainSigmaYXDatass[freqIndex] = new System.Numerics.Complex[nodeCntB];
                }
                System.Numerics.Complex[][][] freqDomainDatasss = {
                    freqDomainUxDatass, freqDomainUyDatass, freqDomainSigmaXXDatass, freqDomainSigmaYXDatass
                };
                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                {
                    for (int dataId = 0; dataId < datasRefUSigmass.Length; dataId++)
                    {
                        double[] datas = datasRefUSigmass[dataId][refIndex][nodeIdB];
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
                            freqDomainDatasss[dataId][freqIndex][nodeIdB] = data;
                        }
                    }
                }

                // モード振幅の計算
                System.Numerics.Complex[] freqDomainAmps = new System.Numerics.Complex[dataCnt];
                var portConditions = World.GetPortConditions(QuantityId);
                for (int freqIndex = 0; freqIndex < dataCnt; freqIndex++)
                {
                    // 周波数
                    double freq = freqs[freqIndex];
                    // 角周波数
                    double omega = 2.0 * Math.PI * freq;

                    // 0周波数は計算できない
                    if (Math.Abs(freq) < IvyFEM.Constants.PrecisionLowerLimit)
                    {
                        freqDomainAmps[freqIndex] = 0;
                        continue;
                    }
                    if (freq < StartFrequencyForSMatrix || freq > EndFrequencyForSMatrix)
                    {
                        freqDomainAmps[freqIndex] = 0;
                        continue;
                    }

                    uint portId = (uint)refIndex;
                    uint _nodeCntB = World.GetPortNodeCount(QuantityId, portId);
                    System.Diagnostics.Debug.Assert(nodeCntB == _nodeCntB);
                    double normalX = 1.0;
                    {
                        var portConditon = portConditions[(int)portId];
                        normalX = portConditon.GetDoubleAdditionalParameters()[0];
                        System.Diagnostics.Debug.Assert(Math.Abs(Math.Abs(normalX) - 1.0) < Constants.PrecisionLowerLimit);
                    }

                    //------------------------------------------------------
                    // モード分布計算
                    //------------------------------------------------------
                    var eigenFEM = new ElasticLambWaveguide1DEigenFEM(World, QuantityId, portId);

                    eigenFEM.Frequency = freq;
                    eigenFEM.NormalX = normalX;
                    eigenFEM.Solve();
                    System.Numerics.Complex[] betas = eigenFEM.Betas;
                    System.Numerics.Complex[][] hUEVecs = eigenFEM.HUEVecs;
                    System.Numerics.Complex[][] hSigmaEVecs = eigenFEM.HSigmaEVecs;

                    //System.Diagnostics.Debug.WriteLine("port = {0} mode Count = {1}", portId, betas.Length);
                    // 基本モード
                    int iMode = 0;
                    System.Numerics.Complex beta = betas[iMode];

                    var srcHF = new System.Numerics.Complex[nodeCntB * uDof];
                    var srcHG = new System.Numerics.Complex[nodeCntB * sDof];
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        srcHF[nodeIdB * uDof] = hUEVecs[iMode][nodeIdB * uDof];
                        srcHF[nodeIdB * uDof + 1] = hUEVecs[iMode][nodeIdB * uDof + 1];
                        srcHG[nodeIdB * sDof] = hSigmaEVecs[iMode][nodeIdB * sDof];
                        srcHG[nodeIdB * sDof + 1] = hSigmaEVecs[iMode][nodeIdB * sDof + 1];
                    }

                    // 振幅分布
                    System.Numerics.Complex[] freqhU = new System.Numerics.Complex[nodeCntB * uDof];
                    System.Numerics.Complex[] freqhSigma = new System.Numerics.Complex[nodeCntB * uDof];
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        freqhU[nodeIdB * uDof] = freqDomainUxDatass[freqIndex][nodeIdB];
                        freqhU[nodeIdB * uDof + 1] = 
                            -1.0 * System.Numerics.Complex.ImaginaryOne * freqDomainUyDatass[freqIndex][nodeIdB];
                        freqhSigma[nodeIdB * sDof] = 
                            System.Numerics.Complex.ImaginaryOne * freqDomainSigmaXXDatass[freqIndex][nodeIdB];
                        freqhSigma[nodeIdB * sDof + 1] = freqDomainSigmaYXDatass[freqIndex][nodeIdB];
                    }
                    // モード振幅の算出
                    System.Numerics.Complex b = eigenFEM.CalcModeAmp(omega, beta, srcHF, srcHG, freqhU, freqhSigma);
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

        public static void CalcCoordUCoordSigma(
            FEWorld world, double[] Uonly,
            out double[] coordUonly,
            out double[] coordSigma)
        {
            uint quantityId = 0;

            int uDof = 2;
            int sDof = 2;
            int coCnt = (int)world.GetCoordCount(quantityId);
            int nodeCnt = (int)world.GetNodeCount(quantityId);

            //-----------------------------------------------
            // u
            coordUonly = new double[coCnt * uDof];
            for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
            {
                int coId = world.Node2Coord(quantityId, iNodeId);
                for (int iDof = 0; iDof < uDof; iDof++)
                {
                    coordUonly[coId * uDof + iDof] = Uonly[iNodeId * uDof + iDof];
                }
            }

            //-----------------------------------------------
            // σの算出
            var coSigmaXXValues = new Dictionary<int, IList<double>>();
            var coSigmaYXValues = new Dictionary<int, IList<double>>();

            var feIds = world.GetTriangleFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = world.GetTriangleFE(quantityId, feId);
                uint elemNodeCnt = triFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = triFE.NodeCoordIds[iNode];
                    int nodeId = world.Coord2Node(quantityId, coId);
                    nodes[iNode] = nodeId;
                }

                Material ma0 = world.GetMaterial(triFE.MaterialId);
                System.Diagnostics.Debug.Assert(
                    ma0 is LinearElasticMaterial ||
                    ma0 is LinearElasticPMLMaterial);
                LinearElasticMaterial ma = null;
                LinearElasticPMLMaterial maPML = null;
                double rho = 0;
                double lambda = 0;
                double mu = 0;
                if (ma0 is LinearElasticMaterial)
                {
                    ma = ma0 as LinearElasticMaterial;
                    rho = ma.MassDensity;
                    lambda = ma.LameLambda;
                    mu = ma.LameMu;
                }
                else if (ma0 is LinearElasticPMLMaterial)
                {
                    maPML = ma0 as LinearElasticPMLMaterial;
                    rho = maPML.MassDensity;
                    lambda = maPML.LameLambda;
                    mu = maPML.LameMu;
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }

                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = triFE.NodeCoordIds[iNode];
                    double[] L = triFE.GetNodeL(iNode);
                    double[][] Nu = triFE.CalcNu(L);
                    double[] Nx = Nu[0];
                    double[] Ny = Nu[1];
                    double sigmaXXValue = 0.0;
                    double sigmaYXValue = 0.0;
                    for (int kNode = 0; kNode < elemNodeCnt; kNode++)
                    {
                        int kNodeId = nodes[kNode];
                        double uxValue = (kNodeId != -1) ? Uonly[kNodeId * uDof] : 0.0;
                        double uyValue = (kNodeId != -1) ? Uonly[kNodeId * uDof + 1] : 0.0;
                        sigmaXXValue +=
                            (lambda + 2.0 * mu) * Nx[kNode] * uxValue + lambda * Ny[kNode] * uyValue;
                        sigmaYXValue +=
                            mu * Nx[kNode] * uyValue + mu * Ny[kNode] * uxValue;
                    }

                    {
                        IList<double> sigmaXXValues = null;
                        if (coSigmaXXValues.ContainsKey(coId))
                        {
                            sigmaXXValues = coSigmaXXValues[coId];
                        }
                        else
                        {
                            sigmaXXValues = new List<double>();
                            coSigmaXXValues.Add(coId, sigmaXXValues);
                        }
                        sigmaXXValues.Add(sigmaXXValue);
                    }
                    {
                        IList<double> sigmaYXValues = null;
                        if (coSigmaYXValues.ContainsKey(coId))
                        {
                            sigmaYXValues = coSigmaYXValues[coId];
                        }
                        else
                        {
                            sigmaYXValues = new List<double>();
                            coSigmaYXValues.Add(coId, sigmaYXValues);
                        }
                        sigmaYXValues.Add(sigmaYXValue);
                    }
                }
            }

            coordSigma = new double[coCnt * sDof];
            for (int coId = 0; coId < coCnt; coId++)
            {
                IList<double>[] sigmaValues = {
                    coSigmaXXValues.ContainsKey(coId) ? coSigmaXXValues[coId] : new List<double>(),
                    coSigmaYXValues.ContainsKey(coId) ? coSigmaYXValues[coId] : new List<double>()
                };
                for (int iDof = 0; iDof < sDof; iDof++)
                {
                    var values = sigmaValues[iDof];
                    double sum = 0.0;
                    foreach (double value in values)
                    {
                        sum += value;
                    }
                    coordSigma[coId * sDof + iDof] = values.Count > 0 ? sum / values.Count : 0.0;
                }
            }
        }
    }
}
