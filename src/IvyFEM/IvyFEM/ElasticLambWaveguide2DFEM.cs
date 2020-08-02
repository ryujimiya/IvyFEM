using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class ElasticLambWaveguide2DFEM : FEM
    {
        public uint UQuantityId { get; private set; } = 0;
        public uint SigmaQuantityId { get; private set; } = 1;

        // Solve
        // input
        public double Frequency { get; set; }
        // output
        public System.Numerics.Complex[] U { get; private set; }
        public System.Numerics.Complex[] CoordU { get; private set; }
        public System.Numerics.Complex[] CoordSigma { get; private set; }
        public System.Numerics.Complex[][] S { get; private set; }
        public ElasticLambWaveguide1DEigenFEM[] EigenFEMs { get; private set; } 

        public ElasticLambWaveguide2DFEM(FEWorld world)
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
            ElasticLambWaveguide1DEigenFEM[] eigenFEMs;
            SetBoundaryCondition(omega, A, B, portCnt, out eigenFEMs);
            System.Diagnostics.Debug.WriteLine("SetBoundaryCondition t = " + (System.Environment.TickCount - t));
            EigenFEMs = eigenFEMs;

            //!!!!!!!!!!!
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
            int uDof = 2;
            int uNodeCnt = (int)World.GetNodeCount(UQuantityId);
            System.Numerics.Complex[] nodeUonly = new System.Numerics.Complex[uNodeCnt * uDof];
            Array.Copy(U, 0, nodeUonly, 0, nodeUonly.Length);

            // uを抽出(座標ベース)
            System.Numerics.Complex[] coordUonly;
            // σを算出(座標ベース)
            System.Numerics.Complex[] coordSigma;
            CalcCoordUCoordSigma(World, Frequency, nodeUonly, out coordUonly, out coordSigma);
            CoordU = coordUonly;
            CoordSigma = coordSigma;

            // Sパラメータを算出
            S = CalcS(omega, nodeUonly, CoordSigma, portCnt, eigenFEMs);
            System.Diagnostics.Debug.WriteLine("CalcS t = " + (System.Environment.TickCount - t));
        }

        private void CalcA(double omega, IvyFEM.Linear.ComplexSparseMatrix A)
        {
            int uDof = 2;

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
                        double kxxVal = (lambda + 2.0 * mu) * sNxNx[row, col] + mu * sNyNy[row, col];
                        double kxyVal = lambda * sNxNy[row, col] + mu * sNyNx[row, col];
                        double kyxVal = lambda * sNyNx[row, col] + mu * sNxNy[row, col];
                        double kyyVal = (lambda + 2.0 * mu) * sNyNy[row, col] + mu * sNxNx[row, col];
                        double mxxVal = rho * sNN[row, col];
                        double myyVal = rho * sNN[row, col];
                        //
                        double axxVal = kxxVal - omega * omega * mxxVal;
                        double axyVal = kxyVal;
                        double ayxVal = kyxVal;
                        double ayyVal = kyyVal - omega * omega * myyVal;

                        A[rowNodeId * uDof, colNodeId * uDof] += axxVal;
                        A[rowNodeId * uDof, colNodeId * uDof + 1] += axyVal;
                        A[rowNodeId * uDof + 1, colNodeId * uDof] += ayxVal;
                        A[rowNodeId * uDof + 1, colNodeId * uDof + 1] += ayyVal;
                    }
                }
            }
        }

        private void SetBoundaryCondition(double omega,
            IvyFEM.Linear.ComplexSparseMatrix A, System.Numerics.Complex[] B,
            uint portCnt, out ElasticLambWaveguide1DEigenFEM[] eigenFEMs)
        {
            int uDof = 2;
            int sDof = 2;

            int uNodeCnt = (int)World.GetNodeCount(UQuantityId);
            int sOffset = uNodeCnt * uDof;

            eigenFEMs = new ElasticLambWaveguide1DEigenFEM[portCnt];
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

                var eigenFEM = new ElasticLambWaveguide1DEigenFEM(World, UQuantityId, portId);
                eigenFEMs[portId] = eigenFEM;

                eigenFEM.Frequency = Frequency;
                eigenFEM.NormalX = normalX;
                eigenFEM.Solve();
                System.Numerics.Complex[] betas = eigenFEM.Betas;
                System.Numerics.Complex[][] hUEVecs = eigenFEM.HUEVecs;
                System.Numerics.Complex[][] hSigmaEVecs = eigenFEM.HSigmaEVecs;

                IvyFEM.Lapack.ComplexMatrix BMat;
                IvyFEM.Lapack.ComplexMatrix Fxx;
                IvyFEM.Lapack.ComplexMatrix Fxy;
                IvyFEM.Lapack.ComplexMatrix Fyx;
                IvyFEM.Lapack.ComplexMatrix Fyy;
                eigenFEM.CalcBoundaryMatrix(
                    omega, betas, hUEVecs, hSigmaEVecs,
                    out BMat, out Fxx, out Fxy, out Fyx, out Fyy);

                // u-sigma
                for (int rowPortNodeId = 0; rowPortNodeId < portNodeCnt; rowPortNodeId++)
                {
                    int rowCoId = World.PortNode2Coord(UQuantityId, portId, rowPortNodeId);
                    int rowUNodeId = World.Coord2Node(UQuantityId, rowCoId);
                    for (int colPortNodeId = 0; colPortNodeId < portNodeCnt; colPortNodeId++)
                    {
                        int colCoId = World.PortNode2Coord(SigmaQuantityId, portId, colPortNodeId);
                        int colSigmaNodeId = World.Coord2Node(SigmaQuantityId, colCoId);

                        A[rowUNodeId * uDof, colSigmaNodeId * sDof + sOffset] =
                            -1.0 * BMat[rowPortNodeId, colPortNodeId];
                        A[rowUNodeId * uDof + 1, colSigmaNodeId * sDof + sOffset + 1] =
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

                        A[rowSigmaNodeId * sDof + sOffset, colUNodeId * uDof] =
                            rowCoId == colCoId ? 1.0 : 0.0;
                        A[rowSigmaNodeId * sDof + sOffset, colUNodeId * uDof + 1] = Fxy[rowPortNodeId, colPortNodeId];
                        A[rowSigmaNodeId * sDof + sOffset + 1, colUNodeId * uDof] = 0.0;
                        A[rowSigmaNodeId * sDof + sOffset + 1, colUNodeId * uDof + 1] = Fyy[rowPortNodeId, colPortNodeId];
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

                        A[rowSigmaNodeId * sDof + sOffset, colSigmaNodeId * sDof + sOffset] =
                            Fxx[rowPortNodeId, colPortNodeId];
                        A[rowSigmaNodeId * sDof + sOffset, colSigmaNodeId * sDof + sOffset + 1] = 0.0;
                        A[rowSigmaNodeId * sDof + sOffset + 1, colSigmaNodeId * sDof + sOffset] =
                            Fyx[rowPortNodeId, colPortNodeId];
                        A[rowSigmaNodeId * sDof + sOffset + 1, colSigmaNodeId * sDof + sOffset + 1] =
                            rowCoId == colCoId ? 1.0 : 0.0;
                    }
                }

                bool isIncidentPort = (portId == World.GetIncidentPortId(UQuantityId));
                if (isIncidentPort)
                {
                    int incidentModeId = World.GetIncidentModeId(UQuantityId);
                    System.Diagnostics.Debug.Assert(incidentModeId != -1);
                    System.Numerics.Complex beta0 = betas[incidentModeId];
                    System.Numerics.Complex[] hUEVec0 = hUEVecs[incidentModeId];
                    System.Numerics.Complex[] hSigmaEVec0 = hSigmaEVecs[incidentModeId];
                    System.Numerics.Complex[] Ix;
                    System.Numerics.Complex[] Iy;
                    eigenFEM.CalcIncidentVec(omega, beta0, hUEVec0, hSigmaEVec0, out Ix, out Iy);

                    // sigma
                    for (int rowPortNodeId = 0; rowPortNodeId < portNodeCnt; rowPortNodeId++)
                    {
                        int rowCoId = World.PortNode2Coord(SigmaQuantityId, portId, rowPortNodeId);
                        int rowSigmaNodeId = World.Coord2Node(SigmaQuantityId, rowCoId);
                        B[rowSigmaNodeId * sDof + sOffset] = Ix[rowPortNodeId];
                        B[rowSigmaNodeId * sDof + sOffset + 1] = Iy[rowPortNodeId];
                    }
                }
            }
        }

        private System.Numerics.Complex[][] CalcS(double omega,
            System.Numerics.Complex[] nodeUonly,
            System.Numerics.Complex[] coordSigma,
            uint portCnt, ElasticLambWaveguide1DEigenFEM[] eigenFEMs)
        {
            var S = new System.Numerics.Complex[portCnt][];

            /*
            int uDof = 2;
            int sDof = 2;
            int sGlobalOffset;
            {
                int uNodeCnt = (int)World.GetNodeCount(UQuantityId);
                sGlobalOffset = uNodeCnt * uDof;
            }
            */
            for (uint portId = 0; portId < portCnt; portId++)
            {
                /*
                //===============================================
                // 方程式を解いて得られたσを利用する場合
                System.Numerics.Complex[] portHU;
                System.Numerics.Complex[] dummyportHSigma;
                GetPortField(portId, nodeUonly, coordSigma, out portHU, out dummyportHSigma);
                //------------------------------------------
                //!!!!!!
                // σ：方程式を解いて得られる値を使用してみる
                int sPortNodeCnt = (int)World.GetPortNodeCount(SigmaQuantityId, portId);
                System.Numerics.Complex[] portHSigma = new System.Numerics.Complex[portHU.Length];
                {
                    for (int sPortNodeId = 0; sPortNodeId < sPortNodeCnt; sPortNodeId++)
                    {
                        int coId = World.PortNode2Coord(SigmaQuantityId, portId, sPortNodeId);
                        int sGlobalNodeId = World.Coord2Node(SigmaQuantityId, coId);
                        // uのポート節点番号に変換
                        int uPortNodeId = World.PortCoord2Node(UQuantityId, portId, coId);
                        if (uPortNodeId == -1)
                        {
                            continue;
                        }

                        // hat σxx = jσxx
                        portHSigma[uPortNodeId * uDof] =
                            System.Numerics.Complex.ImaginaryOne *
                            U[sGlobalNodeId * sDof + sGlobalOffset];
                        portHSigma[uPortNodeId * uDof + 1] =
                            U[sGlobalNodeId * sDof + sGlobalOffset + 1];
                    }
                }
                //------------------------------------------
                //===============================================
                */
                //===============================================
                // uから求めたσを使用する場合
                System.Numerics.Complex[] portHU;
                System.Numerics.Complex[] portHSigma;
                GetPortField(portId, nodeUonly, coordSigma, out portHU, out portHSigma);
                //===============================================

                var eigenFEM = eigenFEMs[portId];
                System.Numerics.Complex[] betas = eigenFEM.Betas;
                System.Numerics.Complex[][] hUEVecs = eigenFEM.HUEVecs;
                System.Numerics.Complex[][] hSigmaEVecs = eigenFEM.HSigmaEVecs;
                int incidentModeId = -1;
                if (World.GetIncidentPortId(UQuantityId) == portId)
                {
                    incidentModeId = (int)World.GetIncidentModeId(UQuantityId);
                }
                System.Numerics.Complex[] S1 = eigenFEM.CalcSMatrix(
                    omega, incidentModeId, betas, hUEVecs, hSigmaEVecs, portHU, portHSigma);
                S[portId] = S1;
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
            int bcNodeCnt = (int)World.GetPortNodeCount(UQuantityId, portId);
            portHU= new System.Numerics.Complex[bcNodeCnt * uDof];
            portHSigma = new System.Numerics.Complex[bcNodeCnt * sDof];
            for (int bcNodeId = 0; bcNodeId < bcNodeCnt; bcNodeId++)
            {
                int coId = World.PortNode2Coord(UQuantityId, portId, bcNodeId);
                int nodeId = World.Coord2Node(UQuantityId, coId);

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
                    System.Numerics.Complex sigmaXXValue = 0.0;
                    System.Numerics.Complex sigmaYXValue = 0.0;
                    for (int kNode = 0; kNode < elemNodeCnt; kNode++)
                    {
                        int kNodeId = nodes[kNode];
                        System.Numerics.Complex uxValue = (kNodeId != -1) ? Uonly[kNodeId * uDof] : 0.0;
                        System.Numerics.Complex uyValue = (kNodeId != -1) ? Uonly[kNodeId * uDof + 1] : 0.0;
                        sigmaXXValue +=
                            (lambda + 2.0 * mu) * Nx[kNode] * uxValue + lambda * Ny[kNode] * uyValue;
                        sigmaYXValue +=
                            mu * Nx[kNode] * uyValue + mu * Ny[kNode] * uxValue;
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
        }
    }
}
