using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class EMWaveguide2DHPlanePMLFEM : FEM
    {
        //---------------------------------
        // ※ポートの設定順序
        // 参照(観測)面、励振源の順
        //---------------------------------

        public uint QuantityId { get; private set; } = 0;
        /// <summary>
        /// 周波数
        /// </summary>
        public double Frequency { get; set; } = 0.0;
        /// <summary>
        /// 1D固有値問題で減衰定数を用いる？
        /// </summary>
        public IList<bool> IsEigen1DUseDecayParameters { get; set; } = new List<bool>(); 
        /// <summary>
        /// 1D固有値問題のクラッド比誘電率
        /// </summary>
        public IList<double> Eigen1DCladdingEps { get; set; } = new List<double>();

        /// <summary>
        /// TEモードで実装した式をTMモードに流用するため
        ///   TEモードの場合は μ0
        ///   TMモードの場合は ε0
        /// </summary>
        public double ReplacedMu0 { get; set; } = Constants.Mu0;

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
        ///  電界
        /// </summary>
        public System.Numerics.Complex[] Ez { get; private set; } = null;
        /// <summary>
        /// Sパラメータ
        /// </summary>
        public System.Numerics.Complex[][] S { get; private set; }
        /// <summary>
        /// 固有値問題 EMWaveguide1DEigenFEM or EMWaveguide1DOepnEigenFEM
        /// </summary>
        public EMWaveguide1DEigenBaseFEM[] EigenBaseFEM { get; private set; }

        public EMWaveguide2DHPlanePMLFEM(FEWorld world)
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

            //--------------------------------------------------------------
            // 残差
            //--------------------------------------------------------------
            CalcB();

            //------------------------------------------------------------------
            // Ezを求める
            //------------------------------------------------------------------
            System.Numerics.Complex[] X;
            Solver.ComplexSolve(out X, A, B);
            Ez = X;

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

            //------------------------------------------------------
            // 剛性行列、質量行列を作成
            //------------------------------------------------------
            int nodeCnt = (int)World.GetNodeCount(QuantityId);
            A = new IvyFEM.Linear.ComplexSparseMatrix(nodeCnt, nodeCnt);

            CalcA(k0, A);

            //------------------------------------------------------
            // モード分布計算
            //------------------------------------------------------
            int refPortCnt = (int)World.GetPortCount(QuantityId) - 1; // 励振源を除く
            System.Diagnostics.Debug.Assert(
                IsEigen1DUseDecayParameters.Count == 0 ||
                IsEigen1DUseDecayParameters.Count == (refPortCnt + 1));
            EigenBaseFEM = new EMWaveguide1DEigenBaseFEM[(refPortCnt + 1)];
            for (int portId = 0; portId < (refPortCnt + 1); portId++)
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

                EigenBaseFEM[portId] = eigen1DBaseFEM;
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

            // バンド幅縮小を座標を用いて行う
            if (Solver is IvyFEM.Linear.LapackEquationSolver)
            {
                var solver = Solver as IvyFEM.Linear.LapackEquationSolver;
                if (solver.IsOrderingToBandMatrix)
                {
                    int[] coIds = new int[nodeCnt];
                    for (int nodeId = 0; nodeId < nodeCnt; nodeId++)
                    {
                        int coId = World.Node2Coord(QuantityId, nodeId);
                        coIds[nodeId] = coId;
                    }
                    double[][] coords = new double[nodeCnt][];
                    for (int nodeId = 0; nodeId < nodeCnt; nodeId++)
                    {
                        int coId = coIds[nodeId];
                        double[] coord = World.GetCoord(QuantityId, coId);
                        coords[nodeId] = coord;
                    }

                    solver.CoordsForBandMatrix = coords;
                }
            }
        }

        private void CalcA(double k0, IvyFEM.Linear.ComplexSparseMatrix _A)
        {
            // 角周波数
            double omega = k0 * Constants.C0;

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
                double epzz = 0;
                double muxx = 0;
                double muyy = 0;
                if (ma0 is DielectricMaterial)
                {
                    ma = ma0 as DielectricMaterial;
                    epzz = ma.Epzz;
                    muxx = ma.Muxx;
                    muyy = ma.Muyy;
                }
                else if (ma0 is DielectricPMLMaterial)
                {
                    maPML = ma0 as DielectricPMLMaterial;
                    epzz = maPML.Epzz;
                    muxx = maPML.Muxx;
                    muyy = maPML.Muyy;
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
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
                System.Numerics.Complex sx = 1.0;
                System.Numerics.Complex sy = 1.0;
                if (maPML != null)
                {
                    bool isXDirection = false;
                    bool isYDirection = false;
                    // X方向PML
                    double sigmaX = 0.0;
                    double epxx = 0.0;
                    // Y方向PML
                    double sigmaY = 0.0;
                    double epyy = 0.0;

                    isXDirection = maPML.IsXDirection();
                    isYDirection = maPML.IsYDirection();
                    if (isXDirection && !isYDirection)
                    {
                        sigmaX = maPML.CalcSigmaX(cPt);
                    }
                    else if (isYDirection && !isXDirection)
                    {
                        sigmaY = maPML.CalcSigmaY(cPt);
                    }
                    else if (isXDirection && isYDirection)
                    {
                        sigmaX = maPML.CalcSigmaX(cPt);
                        sigmaY = maPML.CalcSigmaY(cPt);
                    }
                    else
                    {
                        // 方向がない?
                        System.Diagnostics.Debug.Assert(false);
                    }

                    if (maPML.IsTMMode)
                    {
                        // TMモードのときMuにEpが格納されている
                        epxx = maPML.Muxx;
                        epyy = maPML.Muyy;
                    }
                    else
                    {
                        epxx = maPML.Epxx;
                        epyy = maPML.Epyy;
                    }

                    sx = 1.0 + sigmaX / (System.Numerics.Complex.ImaginaryOne * omega * Constants.Ep0 * epxx);
                    sy = 1.0 + sigmaY / (System.Numerics.Complex.ImaginaryOne * omega * Constants.Ep0 * epyy);
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

                        System.Numerics.Complex a = 
                            (1.0 / muxx) * (sx / sy) * sNyNy[row, col] +
                            (1.0 / muyy) * (sy / sx) * sNxNx[row, col] -
                            k0 * k0 * epzz * sx * sy * sNN[row, col];
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
            int refPortCnt = (int)World.GetPortCount(QuantityId) - 1; // 励振源を除く
            alpha = 0;
            if (IsEigen1DUseDecayParameters.Count == (refPortCnt + 1) &&
                IsEigen1DUseDecayParameters[portId])
            {
                // 減衰定数を考慮した固有値問題
                var eigen1DFEM = new EMWaveguide1DOpenEigenFEM(World, QuantityId, (uint)portId);
                eigen1DBaseFEM = eigen1DFEM;
                eigen1DFEM.CladdingEp = Eigen1DCladdingEps[portId];
                eigen1DFEM.ReplacedMu0 = ReplacedMu0;
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
                eigen1DFEM.ReplacedMu0 = ReplacedMu0;
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
            int nodeCnt = A.RowLength;
            B = new System.Numerics.Complex[nodeCnt];

            int refPortCnt = (int)World.GetPortCount(QuantityId) - 1;  // 励振源分引く

            //--------------------------------------------------------------
            // 励振源
            //--------------------------------------------------------------
            {
                int portId = refPortCnt; // ポートリストの最後の要素が励振境界
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
                        //work[nodeIdB] = (-2.0 * System.Numerics.Complex.ImaginaryOne * betaX) * srcEVec[nodeIdB];
                        // 符号が逆
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
            int refPortCnt = (int)World.GetPortCount(QuantityId) - 1; // 励振源を除く
            int incidentPortId = World.GetIncidentPortId(QuantityId);
            int excitationPortId = refPortCnt;

            // 励振面から入射参照面までの距離と位相差の計算
            var portConditions = World.GetPortConditions(QuantityId);
            PortCondition[] tagtPortConditions = { portConditions[excitationPortId], portConditions[incidentPortId] };
            IList<OpenTK.Vector2d[]> portSEPts = new List<OpenTK.Vector2d[]>();
            foreach (PortCondition portCondition in tagtPortConditions)
            {
                OpenTK.Vector2d[] sePt = new OpenTK.Vector2d[2]; 
                IList<uint> eIds = portCondition.EIds;
                uint eId1 = eIds[0];
                Edge2D e1 = World.Mesh.Cad2D.GetEdge(eId1);
                sePt[0] = e1.GetVertexCoord(true);
                uint eId2 = eIds[eIds.Count - 1];
                Edge2D e2 = World.Mesh.Cad2D.GetEdge(eId2);
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
                double distanceX = Math.Abs(IvyFEM.CadUtils.TriHeight(v1, v2, v3));

                System.Numerics.Complex betaX = SrcBetaXs[excitationPortId];
                // 入射面（ポート1)における振幅を計算
                a = 1.0 * System.Numerics.Complex.Exp(-1.0 * System.Numerics.Complex.ImaginaryOne * betaX * distanceX);
            }

            // Sマトリクスの計算
            var S = new System.Numerics.Complex[refPortCnt][];
            
            for (int refIndex = 0; refIndex < refPortCnt; refIndex++)
            {
                int portId = refIndex;
                int nodeCntB = (int)World.GetPortNodeCount(QuantityId, (uint)portId);
                System.Numerics.Complex[] portEz = GetPortEz((uint)portId, Ez);
                int incidentModeId = -1;
                if (incidentPortId == portId)
                {
                    incidentModeId = (int)World.GetIncidentModeId(QuantityId);
                    // 現状0固定
                    System.Diagnostics.Debug.Assert(incidentModeId == 0);
                }

                EMWaveguide1DEigenBaseFEM eigenBaseFEM = EigenBaseFEM[portId];
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
