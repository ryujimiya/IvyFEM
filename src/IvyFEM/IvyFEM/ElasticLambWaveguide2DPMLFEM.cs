using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Animation;

namespace IvyFEM
{
    public class ElasticLambWaveguide2DPMLFEM : FEM
    {
        //---------------------------------
        // ※ポートの設定順序
        // ポート境界(=PMLの開始位置)、励振源の順
        //---------------------------------

        public uint QuantityId { get; private set; } = 0;

        /// <summary>
        /// [A]
        /// </summary>
        private IvyFEM.Linear.ComplexSparseMatrix A = null;
        /// <summary>
        /// {b}
        /// </summary>
        private System.Numerics.Complex[] B = null;

        /// <summary>
        /// 境界の界の伝搬定数(ポート単位)
        /// </summary>
        public IList<System.Numerics.Complex> SrcBetaXs { get; private set; } = null;
        /// <summary>
        /// 境界の界のモード分布(ポート単位)
        /// </summary>
        public IList<System.Numerics.Complex[]> SrcHFs { get; private set; } = null;
        public IList<System.Numerics.Complex[]> SrcHGs { get; private set; } = null;
        /// <summary>
        /// 境界行列リスト(ポート単位) sNN
        /// </summary>
        private IList<IvyFEM.Lapack.ComplexMatrix> SNNs = null;
        /// <summary>
        /// 境界行列リスト(ポート単位) sNNy
        /// </summary>
        private IList<IvyFEM.Lapack.ComplexMatrix> SNNys = null;

        // Solve
        // input
        public double Frequency { get; set; }
        // output
        public System.Numerics.Complex[] U { get; private set; }
        public System.Numerics.Complex[] CoordU { get; private set; }
        public System.Numerics.Complex[] CoordSigma { get; private set; }
        public System.Numerics.Complex[][] S { get; private set; }
        public ElasticLambWaveguide1DEigenFEM[] EigenFEMs { get; private set; }

        public ElasticLambWaveguide2DPMLFEM(FEWorld world)
        {
            World = world;
        }

        public override void Solve()
        {
            // 周波数
            double freq = Frequency;
            // 角周波数
            double omega = 2.0 * Math.PI * freq;

            //--------------------------------------------------------------
            // 全体行列
            //--------------------------------------------------------------
            CalcA(omega);

            //--------------------------------------------------------------
            // 残差
            //--------------------------------------------------------------
            CalcB(omega);

            //--------------------------------------------------------------
            // 固定境界条件
            //--------------------------------------------------------------
            ComplexSetFixedCadsCondtion(A, B);

            //------------------------------------------------------------------
            // uを求める
            //------------------------------------------------------------------
            System.Numerics.Complex[] X;
            Solver.ComplexSolve(out X, A, B);
            U = X;

            // uを抽出(座標ベース)
            System.Numerics.Complex[] nodeUonly = U;
            System.Numerics.Complex[] coordUonly;
            // σを算出(座標ベース)
            System.Numerics.Complex[] coordSigma;
            CalcCoordUCoordSigma(World, Frequency, nodeUonly, out coordUonly, out coordSigma);
            CoordU = coordUonly;
            CoordSigma = coordSigma;

            //------------------------------------------------------------------
            // Sマトリクスを求める
            //------------------------------------------------------------------
            S = CalcS(omega);
        }

        private void CalcA(double omega)
        {
            int uDof = 2;
            int sDof = 2;

            int nodeCnt = (int)World.GetNodeCount(QuantityId);
            A = new IvyFEM.Linear.ComplexSparseMatrix(nodeCnt * uDof, nodeCnt * uDof);

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
                System.Numerics.Complex sx = 1.0;
                System.Numerics.Complex sy = 1.0;
                if (maPML != null)
                {
                    bool isXDirection = false;
                    bool isYDirection = false;
                    // X方向PML
                    double AX = 1.0;
                    double DX = 0.0;
                    // Y方向PML
                    double AY = 1.0;
                    double DY = 0.0;

                    isXDirection = maPML.IsXDirection();
                    isYDirection = maPML.IsYDirection();
                    if (isXDirection && !isYDirection)
                    {
                        maPML.CalcScalingDampingFactorX(cPt, out AX, out DX);
                    }
                    else if (isYDirection && !isXDirection)
                    {
                        maPML.CalcScalingDampingFactorY(cPt, out AY, out DY);
                    }
                    else if (isXDirection && isYDirection)
                    {
                        maPML.CalcScalingDampingFactorX(cPt, out AX, out DX);
                        maPML.CalcScalingDampingFactorY(cPt, out AY, out DY);
                    }
                    else
                    {
                        // 方向がない?
                        System.Diagnostics.Debug.Assert(false);
                    }

                    //---------------------------------
                    sx = AX * (
                        1.0 - System.Numerics.Complex.ImaginaryOne * (DX / omega));
                    sy = AY * (
                        1.0 - System.Numerics.Complex.ImaginaryOne * (DY / omega));
                    //---------------------------------
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
                        //
                        //-----------------------------------------------------
                        System.Numerics.Complex kxxVal =
                            (lambda + 2.0 * mu) * (sy / sx) * sNxNx[row, col] +
                            mu * (sx / sy) * sNyNy[row, col];
                        System.Numerics.Complex kxyVal =
                            lambda * sNxNy[row, col] +
                            mu * sNyNx[row, col];
                        System.Numerics.Complex kyxVal =
                            lambda * sNyNx[row, col] +
                            mu * sNxNy[row, col];
                        System.Numerics.Complex kyyVal =
                            (lambda + 2.0 * mu) * (sx / sy) * sNyNy[row, col] +
                            mu * (sy / sx) * sNxNx[row, col];
                        System.Numerics.Complex mxxVal = sx * sy * rho * sNN[row, col];
                        System.Numerics.Complex myyVal = sx * sy * rho * sNN[row, col];
                        //-----------------------------------------------------
                        //
                        System.Numerics.Complex axxVal = kxxVal - omega * omega * mxxVal;
                        System.Numerics.Complex axyVal = kxyVal;
                        System.Numerics.Complex ayxVal = kyxVal;
                        System.Numerics.Complex ayyVal = kyyVal - omega * omega * myyVal;

                        A[rowNodeId * uDof, colNodeId * uDof] += axxVal;
                        A[rowNodeId * uDof, colNodeId * uDof + 1] += axyVal;
                        A[rowNodeId * uDof + 1, colNodeId * uDof] += ayxVal;
                        A[rowNodeId * uDof + 1, colNodeId * uDof + 1] += ayyVal;
                    }
                }

                // 回転移動
                // 後片付け
                World.RotAngle = 0.0;
                World.RotOrigin = null;
            }

            //////////////////////////////////////////////////////
            int portCnt = (int)World.GetPortCount(QuantityId) - 1; // 励振源を除く
            SNNs = new List<IvyFEM.Lapack.ComplexMatrix>();
            SNNys = new List<IvyFEM.Lapack.ComplexMatrix>();
            SrcBetaXs = new List<System.Numerics.Complex>();
            SrcHFs = new List<System.Numerics.Complex[]>();
            SrcHGs = new List<System.Numerics.Complex[]>();

            EigenFEMs = new ElasticLambWaveguide1DEigenFEM[(portCnt + 1)];
            var portConditions = World.GetPortConditions(QuantityId);
            for (uint portId = 0; portId < (portCnt + 1); portId++)
            {
                uint portNodeCnt = World.GetPortNodeCount(QuantityId, portId);
                double normalX = 1.0;
                {
                    var portConditon = portConditions[(int)portId];
                    normalX = portConditon.GetComplexAdditionalParameters()[0].Real;
                    System.Diagnostics.Debug.Assert(Math.Abs(Math.Abs(normalX) - 1.0) < Constants.PrecisionLowerLimit);
                }

                //------------------------------------------------------
                // モード分布計算
                //------------------------------------------------------
                var eigenFEM = new ElasticLambWaveguide1DEigenFEM(World, QuantityId, portId);
                EigenFEMs[portId] = eigenFEM;

                eigenFEM.Frequency = Frequency;
                eigenFEM.NormalX = normalX;
                eigenFEM.Solve();
                System.Numerics.Complex[] betas = eigenFEM.Betas;
                System.Numerics.Complex[][] hUEVecs = eigenFEM.HUEVecs;
                System.Numerics.Complex[][] hSigmaEVecs = eigenFEM.HSigmaEVecs;

                System.Diagnostics.Debug.WriteLine("port = {0} mode Count = {1}", portId, betas.Length);
                // 基本モード
                int iMode = 0;
                System.Diagnostics.Debug.Assert(World.GetIncidentModeId(QuantityId) == 0);//現状基本モード固定
                System.Numerics.Complex beta = betas[iMode];
                SrcBetaXs.Add(beta);

                var srcHF = new System.Numerics.Complex[portNodeCnt * uDof];
                var srcHG = new System.Numerics.Complex[portNodeCnt * sDof];
                SrcHFs.Add(srcHF);
                SrcHGs.Add(srcHG);

                for (int nodeIdB = 0; nodeIdB < portNodeCnt; nodeIdB++)
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
        }

        private void CalcB(double omega)
        {
            int uDof = 2;
            int sDof = 2;

            int nodeCnt = A.RowLength / uDof;
            B = new System.Numerics.Complex[nodeCnt * uDof];

            int portCnt = (int)World.GetPortCount(QuantityId) - 1;  // 参照面、励振源分引く

            //--------------------------------------------------------------
            // 励振源
            //--------------------------------------------------------------
            {
                int portId = portCnt; // ポートリストの最後の要素が励振境界
                // 一様媒質とする
                uint maId0;
                {
                    IList<uint> lineFEIds = World.GetPortLineFEIds(QuantityId, (uint)portId);
                    uint lineFEId0 = lineFEIds[0];
                    LineFE lineFE = World.GetLineFE(QuantityId, lineFEId0);
                    maId0 = lineFE.MaterialId;
                }
                var ma0 = World.GetMaterial(maId0);
                System.Diagnostics.Debug.Assert(ma0 is LinearElasticMaterial);
                var ma = ma0 as LinearElasticMaterial;
                double rho = ma.MassDensity;
                double lambda = ma.LameLambda;
                double mu = ma.LameMu;

                var sNN = SNNs[portId];
                var sNNy = SNNys[portId];
                int nodeCntB = sNN.RowLength;
                System.Numerics.Complex srcBetaX = SrcBetaXs[portId];
                System.Numerics.Complex[] srcHF = SrcHFs[portId];
                System.Numerics.Complex[] srcHG = SrcHGs[portId];

                {
                    System.Numerics.Complex betaX = srcBetaX;
                    System.Diagnostics.Debug.Assert(betaX.Real >= 0);
                    System.Diagnostics.Debug.Assert(srcHF.Length == nodeCntB * uDof);
                    System.Diagnostics.Debug.Assert(srcHG.Length == nodeCntB * sDof);

                    System.Numerics.Complex[] workFx = new System.Numerics.Complex[nodeCntB];
                    System.Numerics.Complex[] workFy = new System.Numerics.Complex[nodeCntB];
                    //===============================================================
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        //OK?
                        //// σxxS- == -σxxS+ etc
                        ////--------------------------------------------------------
                        //System.Numerics.Complex hGXX0Value = srcHG[nodeIdB * sDof];
                        //System.Numerics.Complex GYX0Value = srcHG[nodeIdB * sDof + 1];
                        //
                        //workFx[nodeIdB] =
                        //    -2.0 * (-1.0 * System.Numerics.Complex.ImaginaryOne) * hGXX0Value;
                        //workFy[nodeIdB] =
                        //    -2.0 * GYX0Value;
                        ////--------------------------------------------------------

                        //--------------------------------------------------------
                        System.Numerics.Complex hGXX0Value = srcHG[nodeIdB * sDof];
                        System.Numerics.Complex GYX0Value = srcHG[nodeIdB * sDof + 1];

                        workFx[nodeIdB] =
                            -2.0 * (-1.0 * System.Numerics.Complex.ImaginaryOne) * hGXX0Value;
                        workFy[nodeIdB] =
                            -2.0 * GYX0Value;
                        //--------------------------------------------------------
                    }
                    System.Numerics.Complex[] vecIx = sNN * workFx;
                    System.Numerics.Complex[] vecIy = sNN * workFy;
                    //===============================================================

                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);

                        B[nodeId * uDof] += vecIx[nodeIdB];
                        B[nodeId * uDof + 1] += vecIy[nodeIdB];
                    }
                }
            }
        }

        private System.Numerics.Complex[][] CalcS(double omega)
        {
            int portCnt = (int)World.GetPortCount(QuantityId) - 1; // 励振源を除く
            int incidentPortId = World.GetIncidentPortId(QuantityId);
            int excitationPortId = portCnt;

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
            var S = new System.Numerics.Complex[portCnt][];
            
            for (int refIndex = 0; refIndex < portCnt; refIndex++)
            {
                int portId = refIndex; //!!!!!!!!!この実装では参照面はポート境界と一致
                int nodeCntB = (int)World.GetPortNodeCount(QuantityId, (uint)portId);
                System.Numerics.Complex[] portHU;
                System.Numerics.Complex[] portHSigma;
                GetPortField((uint)portId, U, CoordSigma, out portHU, out portHSigma);
                int incidentModeId = -1;
                if (incidentPortId == portId)
                {
                    incidentModeId = World.GetIncidentModeId(QuantityId);
                    // 現状0固定
                    System.Diagnostics.Debug.Assert(incidentModeId == 0);
                }

                ElasticLambWaveguide1DEigenFEM eigenFEM = EigenFEMs[portId];
                double normalX = eigenFEM.NormalX;
                System.Numerics.Complex[] betas = eigenFEM.Betas;
                System.Numerics.Complex[][] hUEVecs = eigenFEM.HUEVecs;
                System.Numerics.Complex[][] hSigmaEVecs = eigenFEM.HSigmaEVecs;
                int modeCnt = betas.Length;
                System.Numerics.Complex[] S1 = new System.Numerics.Complex[modeCnt];
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    System.Numerics.Complex beta = betas[iMode];
                    System.Numerics.Complex[] hUEVec = hUEVecs[iMode];
                    System.Numerics.Complex[] hSigmaEVec = hSigmaEVecs[iMode];
                    System.Numerics.Complex b = eigenFEM.CalcModeAmp(
                        omega, beta, hUEVec, hSigmaEVec, portHU, portHSigma);

                    //!!!!!!!!!!!!!!!
                    //b *= normalX;
                    //!!!!!!!!!!!!!!

                    b *= (1.0 / a);
                    if (incidentModeId == iMode)
                    {
                        //入射面がPML境界のとき
                        //-------------------
                        b += -1.0;
                        //-------------------
                    }
                    S1[iMode] = b;
                }
                S[refIndex] = S1;
            }
            return S;
        }

        private void GetPortField(
            uint portId,
            System.Numerics.Complex[] nodeUonly,
            System.Numerics.Complex[] coordSigma,
            out System.Numerics.Complex[] portHU,
            out System.Numerics.Complex[] portHSigma)
        {
            int uDof = 2;
            int uyOffset = 1;
            int sDof = 2;
            int syxOffset = 1;
            int bcNodeCnt = (int)World.GetPortNodeCount(QuantityId, portId);
            portHU = new System.Numerics.Complex[bcNodeCnt * uDof];
            portHSigma = new System.Numerics.Complex[bcNodeCnt * sDof];
            for (int bcNodeId = 0; bcNodeId < bcNodeCnt; bcNodeId++)
            {
                int coId = World.PortNode2Coord(QuantityId, portId, bcNodeId);
                int nodeId = World.Coord2Node(QuantityId, coId);

                portHU[bcNodeId * uDof] = nodeUonly[nodeId * uDof];
                portHU[bcNodeId * uDof + uyOffset] =
                    -1.0 * System.Numerics.Complex.ImaginaryOne * nodeUonly[nodeId * uDof + uyOffset];
                portHSigma[bcNodeId * sDof] =
                    System.Numerics.Complex.ImaginaryOne * coordSigma[coId * sDof];
                portHSigma[bcNodeId * sDof + syxOffset] = coordSigma[coId * sDof + syxOffset];
            }
        }

        public static void CalcCoordUCoordSigma(
            FEWorld world, double frequency, System.Numerics.Complex[] Uonly,
            out System.Numerics.Complex[] coordUonly,
            out System.Numerics.Complex[] coordSigma)
        {
            uint quantityId = 0;
            // 角周波数
            double omega = 2.0 * Math.PI * frequency;

            int uDof = 2;
            int sDof = 2;
            int coCnt = (int)world.GetCoordCount(quantityId);
            int nodeCnt = (int)world.GetNodeCount(quantityId);

            //-----------------------------------------------
            // u
            coordUonly = new System.Numerics.Complex[coCnt * uDof];
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
            var coSigmaXXValues = new Dictionary<int, IList<System.Numerics.Complex>>();
            var coSigmaYXValues = new Dictionary<int, IList<System.Numerics.Complex>>();

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
                world.RotAngle = rotAngle;
                world.RotOrigin = new double[] { rotOrigin.X, rotOrigin.Y };

                // 重心を求める
                OpenTK.Vector2d cPt;
                {
                    OpenTK.Vector2d[] vertexPts = new OpenTK.Vector2d[triFE.VertexCount];
                    for (int iVertex = 0; iVertex < triFE.VertexCount; iVertex++)
                    {
                        int coId = triFE.VertexCoordIds[iVertex];
                        double[] coord = world.GetCoord(quantityId, coId);
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
                    double AX = 1.0;
                    double DX = 0.0;
                    // Y方向PML
                    double AY = 1.0;
                    double DY = 0.0;

                    isXDirection = maPML.IsXDirection();
                    isYDirection = maPML.IsYDirection();
                    if (isXDirection && !isYDirection)
                    {
                        maPML.CalcScalingDampingFactorX(cPt, out AX, out DX);
                    }
                    else if (isYDirection && !isXDirection)
                    {
                        maPML.CalcScalingDampingFactorY(cPt, out AY, out DY);
                    }
                    else if (isXDirection && isYDirection)
                    {
                        maPML.CalcScalingDampingFactorX(cPt, out AX, out DX);
                        maPML.CalcScalingDampingFactorY(cPt, out AY, out DY);
                    }
                    else
                    {
                        // 方向がない?
                        System.Diagnostics.Debug.Assert(false);
                    }

                    //---------------------------------
                    sx = AX * (
                        1.0 - System.Numerics.Complex.ImaginaryOne * (DX / omega));
                    sy = AY * (
                        1.0 - System.Numerics.Complex.ImaginaryOne * (DY / omega));
                    //---------------------------------
                }

                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = triFE.NodeCoordIds[iNode];
                    double[] L = triFE.GetNodeL(iNode);
                    double[][] Nu = triFE.CalcNu(L);
                    double[] Nx = Nu[0];
                    double[] Ny = Nu[1];
                    System.Numerics.Complex sigmaXXValue = 0.0;
                    System.Numerics.Complex sigmaYXValue = 0.0;
                    for (int kNode = 0; kNode < elemNodeCnt; kNode++)
                    {
                        int kNodeId = nodes[kNode];
                        System.Numerics.Complex uxValue = (kNodeId != -1) ? Uonly[kNodeId * uDof] : 0.0;
                        System.Numerics.Complex uyValue = (kNodeId != -1) ? Uonly[kNodeId * uDof + 1] : 0.0;
                        sigmaXXValue +=
                            (lambda + 2.0 * mu) * (1.0 / sx) * Nx[kNode] * uxValue +
                            lambda * (1.0 / sy) * Ny[kNode] * uyValue;
                        sigmaYXValue +=
                            mu * (1.0 / sx) * Nx[kNode] * uyValue +
                            mu * (1.0 / sy) * Ny[kNode] * uxValue;
                    }

                    {
                        IList<System.Numerics.Complex> sigmaXXValues = null;
                        if (coSigmaXXValues.ContainsKey(coId))
                        {
                            sigmaXXValues = coSigmaXXValues[coId];
                        }
                        else
                        {
                            sigmaXXValues = new List<System.Numerics.Complex>();
                            coSigmaXXValues.Add(coId, sigmaXXValues);
                        }
                        sigmaXXValues.Add(sigmaXXValue);
                    }
                    {
                        IList<System.Numerics.Complex> sigmaYXValues = null;
                        if (coSigmaYXValues.ContainsKey(coId))
                        {
                            sigmaYXValues = coSigmaYXValues[coId];
                        }
                        else
                        {
                            sigmaYXValues = new List<System.Numerics.Complex>();
                            coSigmaYXValues.Add(coId, sigmaYXValues);
                        }
                        sigmaYXValues.Add(sigmaYXValue);
                    }
                }
                // 回転移動
                // 後片付け
                world.RotAngle = 0.0;
                world.RotOrigin = null;
            }

            coordSigma = new System.Numerics.Complex[coCnt * sDof];
            for (int coId = 0; coId < coCnt; coId++)
            {
                IList<System.Numerics.Complex>[] sigmaValues = {
                    coSigmaXXValues.ContainsKey(coId) ? coSigmaXXValues[coId] : new List<System.Numerics.Complex>(),
                    coSigmaYXValues.ContainsKey(coId) ? coSigmaYXValues[coId] : new List<System.Numerics.Complex>()
                };
                for (int iDof = 0; iDof < sDof; iDof++)
                {
                    var values = sigmaValues[iDof];
                    System.Numerics.Complex sum = 0.0;
                    foreach (System.Numerics.Complex value in values)
                    {
                        sum += value;
                    }
                    coordSigma[coId * sDof + iDof] = values.Count > 0 ? sum / values.Count : 0.0;
                }
            }

            //////////////////////////////////////////////////////////////////////////////
            //----------------------------------------------------------------------------
            // ポート境界だけ別計算
            // PML領域と内部領域の境界だが、PML領域からの寄与を無視する
            int portCnt = (int)world.GetPortCount(quantityId) - 1; // 励振源を引く
            // 励振源も含めて
            for (uint portId = 0; portId < (portCnt + 1); portId++)
            {
                int nodeCntB = (int)world.GetPortNodeCount(quantityId, portId);
                var sigmaXXValuesB = new Dictionary<int, IList<System.Numerics.Complex>>();
                var sigmaYXValuesB = new Dictionary<int, IList<System.Numerics.Complex>>();

                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                {
                    int coIdB = world.PortNode2Coord(quantityId, portId, nodeIdB);
                    IList<uint> feIdsB = world.GetTriangleFEIdsFromCoord(quantityId, coIdB);
                    foreach (uint feIdB in feIdsB)
                    {
                        TriangleFE triFEB = world.GetTriangleFE(quantityId, feIdB);
                        uint elemNodeCnt = triFEB.NodeCount;
                        int[] nodes = new int[elemNodeCnt];
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            int coId = triFEB.NodeCoordIds[iNode];
                            int nodeId = world.Coord2Node(quantityId, coId);
                            nodes[iNode] = nodeId;
                        }

                        Material ma0 = world.GetMaterial(triFEB.MaterialId);
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
                            // skip
                            continue;
                        }
                        else
                        {
                            System.Diagnostics.Debug.Assert(false);
                        }

                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            int coId = triFEB.NodeCoordIds[iNode];
                            if (coId != coIdB)
                            {
                                continue;
                            }
                            double[] L = triFEB.GetNodeL(iNode);
                            double[][] Nu = triFEB.CalcNu(L);
                            double[] Nx = Nu[0];
                            double[] Ny = Nu[1];
                            System.Numerics.Complex sigmaXXValue = 0.0;
                            System.Numerics.Complex sigmaYXValue = 0.0;
                            for (int kNode = 0; kNode < elemNodeCnt; kNode++)
                            {
                                int kNodeId = nodes[kNode];
                                System.Numerics.Complex uxValue = (kNodeId != -1) ? Uonly[kNodeId * uDof] : 0.0;
                                System.Numerics.Complex uyValue = (kNodeId != -1) ? Uonly[kNodeId * uDof + 1] : 0.0;
                                sigmaXXValue +=
                                    (lambda + 2.0 * mu) * Nx[kNode] * uxValue +
                                    lambda * Ny[kNode] * uyValue;
                                sigmaYXValue +=
                                    mu * Nx[kNode] * uyValue +
                                    mu * Ny[kNode] * uxValue;
                            }

                            // ポート節点番号ベース
                            {
                                IList<System.Numerics.Complex> worksigmaXXValuesB = null;
                                if (sigmaXXValuesB.ContainsKey(nodeIdB))
                                {
                                    worksigmaXXValuesB = sigmaXXValuesB[nodeIdB];
                                }
                                else
                                {
                                    worksigmaXXValuesB = new List<System.Numerics.Complex>();
                                    sigmaXXValuesB.Add(nodeIdB, worksigmaXXValuesB);
                                }
                                worksigmaXXValuesB.Add(sigmaXXValue);
                            }
                            {
                                IList<System.Numerics.Complex> worksigmaYXValuesB = null;
                                if (sigmaYXValuesB.ContainsKey(nodeIdB))
                                {
                                    worksigmaYXValuesB = sigmaYXValuesB[nodeIdB];
                                }
                                else
                                {
                                    worksigmaYXValuesB = new List<System.Numerics.Complex>();
                                    sigmaYXValuesB.Add(nodeIdB, worksigmaYXValuesB);
                                }
                                worksigmaYXValuesB.Add(sigmaYXValue);
                            }
                        }
                    }
                }

                System.Numerics.Complex[] sigmaB = new System.Numerics.Complex[nodeCntB * sDof];
                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                {
                    IList<System.Numerics.Complex>[] sigmaValuesB = {
                        sigmaXXValuesB.ContainsKey(nodeIdB) ? sigmaXXValuesB[nodeIdB] : new List<System.Numerics.Complex>(),
                        sigmaYXValuesB.ContainsKey(nodeIdB) ? sigmaYXValuesB[nodeIdB] : new List<System.Numerics.Complex>()
                    };
                    for (int iDof = 0; iDof < sDof; iDof++)
                    {
                        var values = sigmaValuesB[iDof];
                        System.Numerics.Complex sum = 0.0;
                        foreach (System.Numerics.Complex value in values)
                        {
                            sum += value;
                        }
                        sigmaB[nodeIdB * sDof + iDof] = values.Count > 0 ? sum / values.Count : 0.0;
                    }
                }
                // 領域全体の配列に上書きする
                for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                {
                    int coId = world.PortNode2Coord(quantityId, portId, nodeIdB);
                    for (int iDof = 0; iDof < sDof; iDof++)
                    {
                        coordSigma[coId * sDof + iDof] = sigmaB[nodeIdB * sDof + iDof]; 
                    }
                }
            }
            //----------------------------------------------------------------------------------------------
            ////////////////////////////////////////////////////////////////////////////////////////////////
        }
    }
}
