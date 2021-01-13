using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class ElasticSHWaveguide2DFirstOrderABCFEM : FEM
    {
        //---------------------------------
        // ※ポートの設定順序
        // ABC境界、参照(観測)面、励振源の順
        //---------------------------------

        public uint QuantityId { get; private set; } = 0;

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

        // Solve
        // input
        public double Frequency { get; set; }
        // output
        public System.Numerics.Complex[] U { get; private set; }
        public System.Numerics.Complex[] CoordU { get; private set; }
        public System.Numerics.Complex[] CoordSigmaZX { get; private set; }
        public System.Numerics.Complex[] CoordSigmaZY { get; private set; }
        public System.Numerics.Complex[][] S { get; private set; }
        public ElasticSHWaveguide1DEigenFEM[] EigenFEMs { get; private set; }

        public ElasticSHWaveguide2DFirstOrderABCFEM(FEWorld world)
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
            System.Numerics.Complex[] coordSigmaZX;
            System.Numerics.Complex[] coordSigmaZY;
            CalcCoordUCoordSigma(World, Frequency, nodeUonly, out coordUonly, out coordSigmaZX, out coordSigmaZY);
            CoordU = coordUonly;
            CoordSigmaZX = coordSigmaZX;
            CoordSigmaZY = coordSigmaZY;

            //------------------------------------------------------------------
            // Sマトリクスを求める
            //------------------------------------------------------------------
            S = CalcS(omega);
        }

        private void CalcA(double omega)
        {
            int nodeCnt = (int)World.GetNodeCount(QuantityId);
            A = new IvyFEM.Linear.ComplexSparseMatrix(nodeCnt, nodeCnt);

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
                System.Diagnostics.Debug.Assert(ma0 is LinearElasticMaterial);
                var ma = ma0 as LinearElasticMaterial;
                double rho = ma.MassDensity;
                double lambda = ma.LameLambda;
                double mu = ma.LameMu;

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
                        //
                        //-----------------------------------------------------
                        System.Numerics.Complex kzzVal = mu * sNxNx[row, col] +
                            mu * sNyNy[row, col];
                        System.Numerics.Complex mzzVal = rho * sNN[row, col];
                        //-----------------------------------------------------
                        //
                        System.Numerics.Complex azzVal = kzzVal - omega * omega * mzzVal;

                        A[rowNodeId, colNodeId] += azzVal;
                    }
                }
            }

            //////////////////////////////////////////////////////
            int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1; // 参照面と励振源を除く
            var Bzzs = new List<IvyFEM.Lapack.ComplexMatrix>();
            SNNs = new List<IvyFEM.Lapack.ComplexMatrix>();
            SrcBetaXs = new List<System.Numerics.Complex>();
            SrcHFs = new List<System.Numerics.Complex[]>();
            SrcHGs = new List<System.Numerics.Complex[]>();

            EigenFEMs = new ElasticSHWaveguide1DEigenFEM[(portCnt + RefPortCount + 1)];
            var portConditions = World.GetPortConditions(QuantityId);
            for (uint portId = 0; portId < (portCnt + RefPortCount + 1); portId++)
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
                var eigenFEM = new ElasticSHWaveguide1DEigenFEM(World, QuantityId, portId);
                EigenFEMs[portId] = eigenFEM;

                eigenFEM.Frequency = Frequency;
                eigenFEM.NormalX = normalX;
                eigenFEM.Solve();
                System.Numerics.Complex[] betas = eigenFEM.Betas;
                System.Numerics.Complex[][] uEVecs = eigenFEM.UEVecs;
                System.Numerics.Complex[][] hSigmaZXEVecs = eigenFEM.HSigmaZXEVecs;

                System.Diagnostics.Debug.WriteLine("port = {0} mode Count = {1}", portId, betas.Length);
                // 基本モード
                int iMode = 0;
                System.Diagnostics.Debug.Assert(World.GetIncidentModeId(QuantityId) == 0);//現状基本モード固定
                System.Numerics.Complex beta = betas[iMode];
                SrcBetaXs.Add(beta);

                var srcF = new System.Numerics.Complex[portNodeCnt];
                var srcHG = new System.Numerics.Complex[portNodeCnt];
                SrcHFs.Add(srcF);
                SrcHGs.Add(srcHG);

                for (int nodeIdB = 0; nodeIdB < portNodeCnt; nodeIdB++)
                {
                    srcF[nodeIdB] = uEVecs[iMode][nodeIdB];
                    srcHG[nodeIdB] = hSigmaZXEVecs[iMode][nodeIdB];
                }

                //----------------------------------------------------------------------
                // Bzz
                //----------------------------------------------------------------------
                System.Numerics.Complex beta0 = betas[0]; // 基本モードと仮定する
                IvyFEM.Lapack.ComplexMatrix Bzz;
                IvyFEM.Lapack.ComplexMatrix sNN;
                eigenFEM.CalcBoundaryMatrixForFirstOrderABC(
                    omega, beta0, out Bzz, out sNN);
                Bzzs.Add(Bzz);
                SNNs.Add(sNN);
            }

            ////////////////////////////////////////////////////////////////
            // 吸収境界条件
            // ABC境界(参照面、励振面を除くポート)
            for (uint portId = 0; portId < portCnt; portId++)
            {
                uint portNodeCnt = World.GetPortNodeCount(QuantityId, portId);
                for (int iNodeIdB = 0; iNodeIdB < portNodeCnt; iNodeIdB++)
                {
                    int iCoId = World.PortNode2Coord(QuantityId, portId, iNodeIdB);
                    int iNodeId = World.Coord2Node(QuantityId, iCoId);

                    for (int jNodeIdB = 0; jNodeIdB < portNodeCnt; jNodeIdB++)
                    {
                        int jCoId = World.PortNode2Coord(QuantityId, portId, jNodeIdB);
                        int jNodeId = World.Coord2Node(QuantityId, jCoId);

                        // zz
                        A[iNodeId, jNodeId] +=
                            -1.0 * Bzzs[(int)portId][iNodeIdB, jNodeIdB];
                    }
                }
            }
        }

        private void CalcB(double omega)
        {
            int nodeCnt = A.RowLength;
            B = new System.Numerics.Complex[nodeCnt];

            int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1; // 参照面と励振源を除く

            //--------------------------------------------------------------
            // 励振源
            //--------------------------------------------------------------
            {
                int portId = portCnt + RefPortCount; // ポートリストの最後の要素が励振境界
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
                int nodeCntB = sNN.RowLength;
                System.Numerics.Complex srcBetaX = SrcBetaXs[portId];
                System.Numerics.Complex[] srcHF = SrcHFs[portId];
                System.Numerics.Complex[] srcHG = SrcHGs[portId];

                {
                    System.Numerics.Complex betaX = srcBetaX;
                    System.Diagnostics.Debug.Assert(betaX.Real >= 0);
                    System.Diagnostics.Debug.Assert(srcHF.Length == nodeCntB);
                    System.Diagnostics.Debug.Assert(srcHG.Length == nodeCntB);

                    System.Numerics.Complex[] workFz = new System.Numerics.Complex[nodeCntB];
                    //===============================================================
                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        //--------------------------------------------------------
                        System.Numerics.Complex uz0Value = srcHF[nodeIdB];

                        workFz[nodeIdB] =
                            -2.0 * (-1.0 * System.Numerics.Complex.ImaginaryOne * betaX * mu) * (omega) * uz0Value;
                        //--------------------------------------------------------
                    }
                    System.Numerics.Complex[] vecIz = sNN * workFz;
                    //===============================================================

                    for (int nodeIdB = 0; nodeIdB < nodeCntB; nodeIdB++)
                    {
                        int coId = World.PortNode2Coord(QuantityId, (uint)portId, nodeIdB);
                        int nodeId = World.Coord2Node(QuantityId, coId);

                        B[nodeId] = vecIz[nodeIdB];
                    }
                }
            }
        }

        private System.Numerics.Complex[][] CalcS(double omega)
        {
            int portCnt = (int)World.GetPortCount(QuantityId) - RefPortCount - 1; // 参照面と励振源を除く
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
            
            for (int refIndex = 0; refIndex < RefPortCount; refIndex++)
            {
                int portId = refIndex + portCnt;
                int nodeCntB = (int)World.GetPortNodeCount(QuantityId, (uint)portId);
                System.Numerics.Complex[] portHU;
                System.Numerics.Complex[] portHSigma;
                GetPortField((uint)portId, U, CoordSigmaZX, out portHU, out portHSigma);
                int incidentModeId = -1;
                if (incidentPortId == portId)
                {
                    incidentModeId = World.GetIncidentModeId(QuantityId);
                    // 現状0固定
                    System.Diagnostics.Debug.Assert(incidentModeId == 0);
                }

                ElasticSHWaveguide1DEigenFEM eigenFEM = EigenFEMs[portId];
                double normalX = eigenFEM.NormalX;
                System.Numerics.Complex[] betas = eigenFEM.Betas;
                System.Numerics.Complex[][] uEVecs = eigenFEM.UEVecs;
                System.Numerics.Complex[][] hSigmaZXEVecs = eigenFEM.HSigmaZXEVecs;
                int modeCnt = betas.Length;
                System.Numerics.Complex[] S1 = new System.Numerics.Complex[modeCnt];
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    System.Numerics.Complex beta = betas[iMode];
                    System.Numerics.Complex[] uEVec = uEVecs[iMode];
                    System.Numerics.Complex[] hSigmaZXEVec = hSigmaZXEVecs[iMode];
                    System.Numerics.Complex b = eigenFEM.CalcModeAmp(
                        omega, beta, uEVec, hSigmaZXEVec, portHU, portHSigma);

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
            System.Numerics.Complex[] coordSigmaZX,
            out System.Numerics.Complex[] portU,
            out System.Numerics.Complex[] portHSigmaZX)
        {
            int bcNodeCnt = (int)World.GetPortNodeCount(QuantityId, portId);
            portU = new System.Numerics.Complex[bcNodeCnt];
            portHSigmaZX = new System.Numerics.Complex[bcNodeCnt];
            for (int bcNodeId = 0; bcNodeId < bcNodeCnt; bcNodeId++)
            {
                int coId = World.PortNode2Coord(QuantityId, portId, bcNodeId);
                int nodeId = World.Coord2Node(QuantityId, coId);

                portU[bcNodeId] = nodeUonly[nodeId];
                portHSigmaZX[bcNodeId] =
                    System.Numerics.Complex.ImaginaryOne * coordSigmaZX[coId];
            }
        }

        public static void CalcCoordUCoordSigma(
            FEWorld world, double frequency, System.Numerics.Complex[] Uonly,
            out System.Numerics.Complex[] coordUonly,
            out System.Numerics.Complex[] coordSigmaZX,
            out System.Numerics.Complex[] coordSigmaZY)
        {
            uint quantityId = 0;
            // 角周波数
            double omega = 2.0 * Math.PI * frequency;

            int coCnt = (int)world.GetCoordCount(quantityId);
            int nodeCnt = (int)world.GetNodeCount(quantityId);

            //-----------------------------------------------
            // u
            coordUonly = new System.Numerics.Complex[coCnt];
            for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
            {
                int coId = world.Node2Coord(quantityId, iNodeId);
                coordUonly[coId] = Uonly[iNodeId];
            }

            //-----------------------------------------------
            // σの算出
            var coSigmaZXValues = new Dictionary<int, IList<System.Numerics.Complex>>();
            var coSigmaZYValues = new Dictionary<int, IList<System.Numerics.Complex>>();

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
                System.Diagnostics.Debug.Assert(ma0 is LinearElasticMaterial);
                LinearElasticMaterial ma = ma0 as LinearElasticMaterial;
                double rho = ma.MassDensity;
                double lambda = ma.LameLambda;
                double mu = ma.LameMu;

                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = triFE.NodeCoordIds[iNode];
                    double[] L = triFE.GetNodeL(iNode);
                    double[][] Nu = triFE.CalcNu(L);
                    double[] Nx = Nu[0];
                    double[] Ny = Nu[1];
                    System.Numerics.Complex sigmaZXValue = 0.0;
                    System.Numerics.Complex sigmaZYValue = 0.0;
                    for (int kNode = 0; kNode < elemNodeCnt; kNode++)
                    {
                        int kNodeId = nodes[kNode];
                        System.Numerics.Complex uzValue = (kNodeId != -1) ? Uonly[kNodeId] : 0.0;
                        sigmaZXValue +=
                            mu * Nx[kNode] * uzValue;
                        sigmaZYValue +=
                            mu * Ny[kNode] * uzValue;
                    }

                    {
                        IList<System.Numerics.Complex> sigmaZXValues = null;
                        if (coSigmaZXValues.ContainsKey(coId))
                        {
                            sigmaZXValues = coSigmaZXValues[coId];
                        }
                        else
                        {
                            sigmaZXValues = new List<System.Numerics.Complex>();
                            coSigmaZXValues.Add(coId, sigmaZXValues);
                        }
                        sigmaZXValues.Add(sigmaZXValue);
                    }
                    {
                        IList<System.Numerics.Complex> sigmaZYValues = null;
                        if (coSigmaZYValues.ContainsKey(coId))
                        {
                            sigmaZYValues = coSigmaZYValues[coId];
                        }
                        else
                        {
                            sigmaZYValues = new List<System.Numerics.Complex>();
                            coSigmaZYValues.Add(coId, sigmaZYValues);
                        }
                        sigmaZYValues.Add(sigmaZYValue);
                    }
                }
            }

            coordSigmaZX = new System.Numerics.Complex[coCnt];
            coordSigmaZY = new System.Numerics.Complex[coCnt];
            for (int coId = 0; coId < coCnt; coId++)
            {
                IList<System.Numerics.Complex>[] sigmaValues = {
                    coSigmaZXValues.ContainsKey(coId) ? coSigmaZXValues[coId] : new List<System.Numerics.Complex>(),
                    coSigmaZYValues.ContainsKey(coId) ? coSigmaZYValues[coId] : new List<System.Numerics.Complex>()
                };
                System.Numerics.Complex[] averageValues = { 0.0, 0.0 };
                for (int iDof = 0; iDof < sigmaValues.Length; iDof++)
                {
                    var values = sigmaValues[iDof];
                    System.Numerics.Complex sum = 0.0;
                    foreach (System.Numerics.Complex value in values)
                    {
                        sum += value;
                    }
                    averageValues[iDof] = values.Count > 0 ? sum / values.Count : 0.0;
                }
                coordSigmaZX[coId] = averageValues[0];
                coordSigmaZY[coId] = averageValues[1];
            }

        }
    }
}
