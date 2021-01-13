using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class ElasticSHWaveguide2DFEM : FEM
    {
        public uint UQuantityId { get; private set; } = 0;
        public uint SigmaQuantityId { get; private set; } = 1;

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

        public ElasticSHWaveguide2DFEM(FEWorld world)
        {
            World = world;
        }

        public override void Solve()
        {
            int t;
            U = null;
            S = null;
            EigenFEMs = null;

            // 角周波数
            double omega = 2.0 * Math.PI * Frequency;

            t = System.Environment.TickCount;
            int nodeCnt = 0;
            int quantityCnt = World.GetQuantityCount();
            for (uint quantityId = 0; quantityId < quantityCnt; quantityId++)
            {
                int quantityDof = (int)World.GetDof(quantityId);
                int quantityNodeCnt = (int)World.GetNodeCount(quantityId);
                nodeCnt += quantityNodeCnt * quantityDof;
            }

            var A = new IvyFEM.Linear.ComplexSparseMatrix(nodeCnt, nodeCnt);
            var B = new System.Numerics.Complex[nodeCnt];
            CalcA(omega, A);
            System.Diagnostics.Debug.WriteLine("CalcA t = " + (System.Environment.TickCount - t));

            t = System.Environment.TickCount;
            uint portCnt = World.GetPortCount(UQuantityId);
            ElasticSHWaveguide1DEigenFEM[] eigenFEMs;
            SetBoundaryCondition(omega, A, B, portCnt, out eigenFEMs);
            System.Diagnostics.Debug.WriteLine("SetBoundaryCondition t = " + (System.Environment.TickCount - t));
            EigenFEMs = eigenFEMs;

            //--------------------------------------------------------------
            // 固定境界条件
            //--------------------------------------------------------------
            ComplexSetFixedCadsCondtion(A, B);

            t = System.Environment.TickCount;
            //----------------------------------
            System.Numerics.Complex[] X;
            Solver.ComplexSolve(out X, A, B);
            U = X;
            //----------------------------------
            System.Diagnostics.Debug.WriteLine("Solve t = " + (System.Environment.TickCount - t));

            t = System.Environment.TickCount;
            // Uだけを抽出
            int uNodeCnt = (int)World.GetNodeCount(UQuantityId);
            System.Numerics.Complex[] nodeUonly = new System.Numerics.Complex[uNodeCnt];
            Array.Copy(U, 0, nodeUonly, 0, nodeUonly.Length);

            // uを抽出(座標ベース)
            System.Numerics.Complex[] coordUonly;
            // σを算出(座標ベース)
            System.Numerics.Complex[] coordSigmaZX;
            System.Numerics.Complex[] coordSigmaZY;
            CalcCoordUCoordSigma(World, Frequency, nodeUonly, out coordUonly, out coordSigmaZX, out coordSigmaZY);
            CoordU = coordUonly;
            CoordSigmaZX = coordSigmaZX;
            CoordSigmaZY = coordSigmaZY;

            // Sパラメータを算出
            S = CalcS(omega, nodeUonly, CoordSigmaZX, portCnt, eigenFEMs);
            System.Diagnostics.Debug.WriteLine("CalcS t = " + (System.Environment.TickCount - t));
        }

        private void CalcA(double omega, IvyFEM.Linear.ComplexSparseMatrix A)
        {
            IList<uint> feIds = World.GetTriangleFEIds(UQuantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = World.GetTriangleFE(UQuantityId, feId);
                uint elemNodeCnt = triFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = triFE.NodeCoordIds[iNode];
                    int nodeId = World.Coord2Node(UQuantityId, coId);
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
                        double kzzVal = mu * sNxNx[row, col] + mu * sNyNy[row, col];
                        double mzzVal = rho * sNN[row, col];
                        //
                        double azzVal = kzzVal - omega * omega * mzzVal;

                        A[rowNodeId, colNodeId] += azzVal;
                    }
                }
            }
        }

        private void SetBoundaryCondition(double omega,
            IvyFEM.Linear.ComplexSparseMatrix A, System.Numerics.Complex[] B,
            uint portCnt, out ElasticSHWaveguide1DEigenFEM[] eigenFEMs)
        {
            int uNodeCnt = (int)World.GetNodeCount(UQuantityId);
            int sOffset = uNodeCnt;

            eigenFEMs = new ElasticSHWaveguide1DEigenFEM[portCnt];
            var portConditions = World.GetPortConditions(UQuantityId);
            for (uint portId = 0; portId < portCnt; portId++)
            {
                uint portNodeCnt = World.GetPortNodeCount(UQuantityId, portId);
                double normalX = 1.0;
                {
                    var portConditon = portConditions[(int)portId];
                    normalX = portConditon.GetComplexAdditionalParameters()[0].Real;
                    System.Diagnostics.Debug.Assert(Math.Abs(Math.Abs(normalX) - 1.0) < Constants.PrecisionLowerLimit);
                }

                var eigenFEM = new ElasticSHWaveguide1DEigenFEM(World, UQuantityId, portId);
                eigenFEMs[portId] = eigenFEM;

                eigenFEM.Frequency = Frequency;
                eigenFEM.NormalX = normalX;
                eigenFEM.Solve();
                System.Numerics.Complex[] betas = eigenFEM.Betas;
                System.Numerics.Complex[][] uEVecs = eigenFEM.UEVecs;
                System.Numerics.Complex[][] hSigmaEVecs = eigenFEM.HSigmaZXEVecs;

                IvyFEM.Lapack.ComplexMatrix BMat;
                IvyFEM.Lapack.ComplexMatrix Fzz;
                eigenFEM.CalcBoundaryMatrix(
                    omega, betas, uEVecs, hSigmaEVecs,
                    out BMat, out Fzz);

                // u-sigma
                for (int rowPortNodeId = 0; rowPortNodeId < portNodeCnt; rowPortNodeId++)
                {
                    int rowCoId = World.PortNode2Coord(UQuantityId, portId, rowPortNodeId);
                    int rowUNodeId = World.Coord2Node(UQuantityId, rowCoId);
                    for (int colPortNodeId = 0; colPortNodeId < portNodeCnt; colPortNodeId++)
                    {
                        int colCoId = World.PortNode2Coord(SigmaQuantityId, portId, colPortNodeId);
                        int colSigmaNodeId = World.Coord2Node(SigmaQuantityId, colCoId);

                        A[rowUNodeId, colSigmaNodeId + sOffset] =
                            -1.0 * BMat[rowPortNodeId, colPortNodeId];
                    }
                }
                // sigma-u
                for (int rowPortNodeId = 0; rowPortNodeId < portNodeCnt; rowPortNodeId++)
                {
                    int rowCoId = World.PortNode2Coord(SigmaQuantityId, portId, rowPortNodeId);
                    int rowSigmaNodeId = World.Coord2Node(SigmaQuantityId, rowCoId);
                    for (int colPortNodeId = 0; colPortNodeId < portNodeCnt; colPortNodeId++)
                    {
                        int colCoId = World.PortNode2Coord(UQuantityId, portId, colPortNodeId);
                        int colUNodeId = World.Coord2Node(UQuantityId, colCoId);

                        A[rowSigmaNodeId + sOffset, colUNodeId] = Fzz[rowPortNodeId, colPortNodeId];
                    }
                }

                // sigma-sigma
                for (int rowPortNodeId = 0; rowPortNodeId < portNodeCnt; rowPortNodeId++)
                {
                    int rowCoId = World.PortNode2Coord(SigmaQuantityId, portId, rowPortNodeId);
                    int rowSigmaNodeId = World.Coord2Node(SigmaQuantityId, rowCoId);
                    for (int colPortNodeId = 0; colPortNodeId < portNodeCnt; colPortNodeId++)
                    {
                        int colCoId = World.PortNode2Coord(SigmaQuantityId, portId, colPortNodeId);
                        int colSigmaNodeId = World.Coord2Node(SigmaQuantityId, colCoId);

                        A[rowSigmaNodeId + sOffset, colSigmaNodeId + sOffset] =
                            rowCoId == colCoId ? 1.0 : 0.0;
                    }
                }

                bool isIncidentPort = (portId == World.GetIncidentPortId(UQuantityId));
                if (isIncidentPort)
                {
                    int incidentModeId = World.GetIncidentModeId(UQuantityId);
                    System.Diagnostics.Debug.Assert(incidentModeId != -1);
                    System.Numerics.Complex beta0 = betas[incidentModeId];
                    System.Numerics.Complex[] uEVec0 = uEVecs[incidentModeId];
                    System.Numerics.Complex[] hSigmaEVec0 = hSigmaEVecs[incidentModeId];
                    System.Numerics.Complex[] Iz;
                    eigenFEM.CalcIncidentVec(omega, beta0, uEVec0, hSigmaEVec0, out Iz);

                    // sigma
                    for (int rowPortNodeId = 0; rowPortNodeId < portNodeCnt; rowPortNodeId++)
                    {
                        int rowCoId = World.PortNode2Coord(SigmaQuantityId, portId, rowPortNodeId);
                        int rowSigmaNodeId = World.Coord2Node(SigmaQuantityId, rowCoId);
                        B[rowSigmaNodeId + sOffset] = Iz[rowPortNodeId];
                    }
                }
            }
        }

        private System.Numerics.Complex[][] CalcS(double omega,
            System.Numerics.Complex[] nodeUonly,
            System.Numerics.Complex[] coordSigmaZX,
            uint portCnt, ElasticSHWaveguide1DEigenFEM[] eigenFEMs)
        {
            var S = new System.Numerics.Complex[portCnt][];

            for (uint portId = 0; portId < portCnt; portId++)
            {
                //===============================================
                // uから求めたσを使用する場合
                System.Numerics.Complex[] portU;
                System.Numerics.Complex[] portHSigmaZX;
                GetPortField(portId, nodeUonly, coordSigmaZX, out portU, out portHSigmaZX);
                //===============================================

                var eigenFEM = eigenFEMs[portId];
                System.Numerics.Complex[] betas = eigenFEM.Betas;
                System.Numerics.Complex[][] uVecs = eigenFEM.UEVecs;
                System.Numerics.Complex[][] hSigmaZXEVecs = eigenFEM.HSigmaZXEVecs;
                int incidentModeId = -1;
                if (World.GetIncidentPortId(UQuantityId) == portId)
                {
                    incidentModeId = (int)World.GetIncidentModeId(UQuantityId);
                }
                System.Numerics.Complex[] S1 = eigenFEM.CalcSMatrix(
                    omega, incidentModeId, betas, uVecs, hSigmaZXEVecs, portU, portHSigmaZX);
                S[portId] = S1;
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
            int uDof = 1;
            int sDof = 1;
            int bcNodeCnt = (int)World.GetPortNodeCount(UQuantityId, portId);
            portU= new System.Numerics.Complex[bcNodeCnt * uDof];
            portHSigmaZX = new System.Numerics.Complex[bcNodeCnt * sDof];
            for (int bcNodeId = 0; bcNodeId < bcNodeCnt; bcNodeId++)
            {
                int coId = World.PortNode2Coord(UQuantityId, portId, bcNodeId);
                int nodeId = World.Coord2Node(UQuantityId, coId);

                portU[bcNodeId * uDof] = nodeUonly[nodeId * uDof];
                portHSigmaZX[bcNodeId * sDof] =
                    System.Numerics.Complex.ImaginaryOne * coordSigmaZX[coId * sDof];
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
                var ma = ma0 as LinearElasticMaterial;
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
