using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class PCWaveguide2DPBCFEM : FEM
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

        public IList<PCWaveguidePortInfo> WgPortInfos { get; set; } = null;
        /// <summary>
        /// TMモード？
        /// </summary>
        public bool IsTMMode { get; set; } = false;
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
        public IList<System.Numerics.Complex[]> SrcFEVecs { get; private set; } = null;
        /// <summary>
        /// 境界の界（周期構造用に修正)のモード分布(ポート単位)
        /// </summary>
        public IList<System.Numerics.Complex[]> SrcModifyFEVecs { get; private set; } = null;

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
        public PCWaveguide2DEigenFEM[] EigenFEMs { get; private set; }

        public PCWaveguide2DPBCFEM(FEWorld world)
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
            SrcBetaXs = new List<double>();
            SrcFEVecs = new List<System.Numerics.Complex[]>();
            SrcModifyFEVecs = new List<System.Numerics.Complex[]>();

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
            int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1; // 参照面と励振源を除く

            EigenFEMs = new PCWaveguide2DEigenFEM[(portCnt + RefPortCount + 1)];
            for (int portId = 0; portId < (portCnt + RefPortCount + 1); portId++)
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
                    portId, freq,
                    out ryy1D, out txx1D, out uzz1D, out betas, out eVecs, out fxEVecs, out eigenFEM);

                EigenFEMs[portId] = eigenFEM;
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
                SrcFEVecs.Add(fVec);
                SrcModifyFEVecs.Add(fVecModify);
            }

            ////////////////////////////////////////////////////////////////////////////////////////
            // 周期境界 境界1
            for (int portId = 0; portId < portCnt; portId++)
            {
                var wgPortInfo = WgPortInfos[portId];
                double srcBetaX = SrcBetaXs[portId];
                double periodicDistance = wgPortInfo.PeriodicDistanceX;
                bool isPortBc2Reverse = wgPortInfo.IsPortBc2Reverse;
                var bcNodes1 = wgPortInfo.BcNodess[0]; // 境界1
                var bcNodes2 = wgPortInfo.BcNodess[1]; // 境界2
                int nodeCntB = bcNodes1.Count;
                for (int iB1 = 0; iB1 < nodeCntB; iB1++)
                {
                    int rowPortNodeId = bcNodes1[iB1];
                    int rowCoId = World.PortNode2Coord(QuantityId, (uint)portId, rowPortNodeId);
                    int rowNodeId = World.Coord2Node(QuantityId, rowCoId);
                    for (int colNodeId = 0; colNodeId < nodeCnt; colNodeId++)
                    {
                        A[rowNodeId, colNodeId] = 0.0;
                    }
                    // 境界1
                    {
                        int colNodeId = rowNodeId;
                        A[rowNodeId, colNodeId] = 1.0;
                    }
                    // 境界2
                    {
                        int jB1 = iB1;
                        int jB2 = isPortBc2Reverse ? (nodeCntB - 1 - jB1) : jB1;
                        int colPortNodeId = bcNodes2[jB2];
                        int colCoId = World.PortNode2Coord(QuantityId, (uint)portId, colPortNodeId);
                        int colNodeId = World.Coord2Node(QuantityId, colCoId);
                        A[rowNodeId, colNodeId] = -1.0 * System.Numerics.Complex.Exp(
                            -1.0 * System.Numerics.Complex.ImaginaryOne * srcBetaX * periodicDistance);
                    }
                }
            }
        }

        private void CalcA(double k0,
            IvyFEM.Linear.ComplexSparseMatrix A)
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
                        double a = maPxx * sNyNy[row, col] + maPyy * sNxNx[row, col] -
                            k0 * k0 * maQzz * sNN[row, col];

                        A[rowNodeId, colNodeId] += (System.Numerics.Complex)a;
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

            int nodeCnt = A.RowLength;
            B = new System.Numerics.Complex[nodeCnt];

            int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1; // 参照面と励振源を除く

            //--------------------------------------------------------------
            // 励振源
            //--------------------------------------------------------------
            {
                int portId = portCnt + RefPortCount; // ポートリストの最後の要素が励振境界
                var Qb = Qbs[portId];
                int nodeCntB = Qb.RowLength;
                System.Numerics.Complex srcBetaX = SrcBetaXs[portId];
                System.Numerics.Complex[] srcModifyEVec = SrcModifyFEVecs[portId];

                {
                    System.Numerics.Complex betaX = srcBetaX;
                    System.Diagnostics.Debug.Assert(srcModifyEVec.Length == nodeCntB);
                    System.Diagnostics.Debug.Assert(betaX.Real >= 0);
                    System.Numerics.Complex[] work = new System.Numerics.Complex[nodeCntB];
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        work[nodeIdB] = (2.0 * System.Numerics.Complex.ImaginaryOne * betaX) * srcModifyEVec[nodeIdB];
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
                System.Diagnostics.Debug.Assert(portCondition.IsPeriodic);
                OpenTK.Vector2d[] sePt = new OpenTK.Vector2d[2];
                IList<uint> eIds = portCondition.BcEIdsForPeriodic1;
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

                // FIXME: 三角形格子TMモードだと入射が逆になっているようだ なぜ？
                if (IsTMMode)
                {
                    a *= -1.0;
                }
            }

            // Sマトリクスの計算
            var S = new System.Numerics.Complex[RefPortCount][];

            for (int refIndex = 0; refIndex < RefPortCount; refIndex++)
            {
                int portId = portCnt + refIndex;
                System.Numerics.Complex[] portEzB1 = GetPortEzB1((uint)portId, Ez);
                int incidentModeId = -1;
                if (incidentPortId == portId)
                {
                    incidentModeId = (int)World.GetIncidentModeId(QuantityId);
                    // 現状0固定
                    System.Diagnostics.Debug.Assert(incidentModeId == 0);
                }

                PCWaveguide2DEigenFEM eigenFEM = EigenFEMs[portId];
                System.Numerics.Complex[] betas = eigenFEM.Betas;
                System.Numerics.Complex[][] ezEVecs = eigenFEM.BcEVecs;
                System.Numerics.Complex[][] ezFxEVecs = eigenFEM.BcFxEVecs;
                int modeCnt = betas.Length;
                System.Numerics.Complex[] S1 = new System.Numerics.Complex[modeCnt];
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    System.Numerics.Complex b = eigenFEM.CalcModeAmp(omega, iMode, betas, ezEVecs, ezFxEVecs, portEzB1);
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

        private System.Numerics.Complex[] GetPortEzB1(uint portId, System.Numerics.Complex[] Ez)
        {
            var wgPortInfo = WgPortInfos[(int)portId];
            var bcNodes1 = wgPortInfo.BcNodess[0];
            int nodeCntB = bcNodes1.Count;
            System.Numerics.Complex[] portEzB1 = new System.Numerics.Complex[nodeCntB];
            for (int row = 0; row < nodeCntB; row++)
            {
                int portNodeId = bcNodes1[row];
                int coId = World.PortNode2Coord(QuantityId, portId, portNodeId);
                int nodeId = World.Coord2Node(QuantityId, coId);
                portEzB1[row] = Ez[nodeId];
            }
            return portEzB1;
        }
    }
}
