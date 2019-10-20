using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class PCWaveguide2DFEM : FEM
    {
        public uint QuantityId { get; private set; } = 0;

        public IList<PCWaveguidePortInfo> WgPortInfos { get; set; } = null;

        // Solve
        // input
        public double Frequency { get; set; }
        // output
        public System.Numerics.Complex[] Ez { get; private set; }
        public System.Numerics.Complex[][] S { get; private set; }
        public PCWaveguide2DEigenFEM[] EigenFEMs { get; private set; }

        public PCWaveguide2DFEM(FEWorld world)
        {
            World = world;
        }

        public override void Solve()
        {
            int t;
            Ez = null;
            S = null;
            EigenFEMs = null;

            // 角周波数
            double omega = 2.0 * Math.PI * Frequency;
            // 波数
            double k0 = omega / Constants.C0;

            t = System.Environment.TickCount;
            int nodeCnt = (int)World.GetNodeCount(QuantityId);
            var A = new IvyFEM.Linear.ComplexSparseMatrix(nodeCnt, nodeCnt);
            var B = new System.Numerics.Complex[nodeCnt];
            CalcA(k0, A);
            System.Diagnostics.Debug.WriteLine("CalcAB t = " + (System.Environment.TickCount - t));

            t = System.Environment.TickCount;
            uint portCnt = World.GetPortCount(QuantityId);
            PCWaveguide2DEigenFEM[] eigenFEMs;
            SetBoundaryCondition(omega, A, B, portCnt, out eigenFEMs);
            System.Diagnostics.Debug.WriteLine("SetBoundaryCondition t = " + (System.Environment.TickCount - t));
            EigenFEMs = eigenFEMs;

            t = System.Environment.TickCount;
            //----------------------------------
            System.Numerics.Complex[] X;
            Solver.ComplexSolve(out X, A, B);
            Ez = X;
            //----------------------------------
            System.Diagnostics.Debug.WriteLine("Solve t = " + (System.Environment.TickCount - t));

            t = System.Environment.TickCount;
            S = CalcS(omega, Ez, portCnt, eigenFEMs);
            System.Diagnostics.Debug.WriteLine("CalcS t = " + (System.Environment.TickCount - t));
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
                        double a = (1.0 / ma.Muxx) * sNyNy[row, col] + (1.0 / ma.Muyy) * sNxNx[row, col] -
                            (k0 * k0 * ma.Epzz) * sNN[row, col];

                        A[rowNodeId, colNodeId] += (System.Numerics.Complex)a;
                    }
                }
            }
        }

        private void SetBoundaryCondition(double omega,
            IvyFEM.Linear.ComplexSparseMatrix A, System.Numerics.Complex[] B,
            uint portCnt, out PCWaveguide2DEigenFEM[] eigenFEMs)
        {
            eigenFEMs = new PCWaveguide2DEigenFEM[portCnt];
            for (uint portId = 0; portId < portCnt; portId++)
            {
                PCWaveguidePortInfo wgPortInfo = WgPortInfos[(int)portId];
                var eigenFEM = new PCWaveguide2DEigenFEM(World, QuantityId, portId, wgPortInfo);
                eigenFEMs[portId] = eigenFEM;

                eigenFEM.Frequency = Frequency;
                eigenFEM.Solve();
                System.Numerics.Complex[] betas = eigenFEM.Betas;
                System.Numerics.Complex[][] ezEVecs = eigenFEM.BcEVecs;
                System.Numerics.Complex[][] ezXEVecs = eigenFEM.BcFxEVecs;
                int modeCnt = betas.Length;

                IvyFEM.Lapack.ComplexMatrix b = eigenFEM.CalcBoundaryMatrix(omega, betas, ezEVecs, ezXEVecs);
                uint bcIndex = 0;
                IList<int> bcNodes = eigenFEM.WgPortInfo.BcNodess[(int)bcIndex];
                int bcNodeCnt = bcNodes.Count;
                for (int row = 0; row < bcNodeCnt; row++)
                {
                    int rowPortNodeId = bcNodes[row];
                    //int rowCoId = World.PeriodicPortBcNode2Coord(QuantityId, portId, bcIndex, rowPortNodeId);
                    int rowCoId = World.PortNode2Coord(QuantityId, portId, rowPortNodeId);
                    int rowNodeId = World.Coord2Node(QuantityId, rowCoId);
                    for (int col = 0; col < bcNodeCnt; col++)
                    {
                        int colPortNodeId = bcNodes[col];
                        //int colCoId = World.PeriodicPortBcNode2Coord(QuantityId, portId, bcIndex, colPortNodeId);
                        int colCoId = World.PortNode2Coord(QuantityId, portId, colPortNodeId);
                        int colNodeId = World.Coord2Node(QuantityId, colCoId);

                        A[rowNodeId, colNodeId] += b[row, col];
                    }
                }

                bool isIncidentPort = (portId == World.GetIncidentPortId(QuantityId));
                int incidentModeId = World.GetIncidentModeId(QuantityId);
                if (isIncidentPort && incidentModeId < modeCnt)
                {
                    System.Diagnostics.Debug.Assert(incidentModeId != -1);
                    System.Numerics.Complex beta0 = betas[incidentModeId];
                    System.Numerics.Complex[] ezEVec0 = ezEVecs[incidentModeId];
                    System.Numerics.Complex[] ezXEVec0 = ezXEVecs[incidentModeId];
                    System.Numerics.Complex[] I = eigenFEM.CalcIncidentVec(beta0, ezEVec0, ezXEVec0);
                    for (int row = 0; row < bcNodeCnt; row++)
                    {
                        int rowPortNodeId = bcNodes[row];
                        //int rowCoId = World.PeriodicPortBcNode2Coord(QuantityId, portId, bcIndex, rowPortNodeId);
                        int rowCoId = World.PortNode2Coord(QuantityId, portId, rowPortNodeId);
                        int rowNodeId = World.Coord2Node(QuantityId, rowCoId);

                        B[rowNodeId] += I[row];
                    }
                }
            }
        }

        private System.Numerics.Complex[][] CalcS(double omega,
            System.Numerics.Complex[] Ez,
            uint portCnt, PCWaveguide2DEigenFEM[] eigenFEMs)
        {
            var S = new System.Numerics.Complex[portCnt][];
            for (uint portId = 0; portId < portCnt; portId++)
            {
                System.Numerics.Complex[] portEz = GetPortEz(portId, Ez);
                var eigenFEM = eigenFEMs[portId];
                System.Numerics.Complex[] betas = eigenFEM.Betas;
                System.Numerics.Complex[][] ezEVecs = eigenFEM.BcEVecs;
                System.Numerics.Complex[][] ezXEVecs = eigenFEM.BcFxEVecs;
                int incidentModeId = -1;
                if (World.GetIncidentPortId(QuantityId) == portId)
                {
                    incidentModeId = (int)World.GetIncidentModeId(QuantityId);
                }
                System.Numerics.Complex[] S1 = eigenFEM.CalcSMatrix(
                    omega, incidentModeId, betas, ezEVecs, ezXEVecs, portEz);
                S[portId] = S1;
            }
            return S;
        }

        private System.Numerics.Complex[] GetPortEz(uint portId, System.Numerics.Complex[] Ez)
        {
            var eigenFEM = EigenFEMs[portId];

            uint bcIndex = 0;
            IList<int> bcNodes = eigenFEM.WgPortInfo.BcNodess[(int)bcIndex];
            int bcNodeCnt = bcNodes.Count;
            System.Numerics.Complex[] portEz = new System.Numerics.Complex[bcNodeCnt];
            for (int row = 0; row < bcNodeCnt; row++)
            {
                int rowPortNodeId = bcNodes[row];
                //int rowCoId = World.PeriodicPortBcNode2Coord(QuantityId, portId, bcIndex, rowPortNodeId);
                int rowCoId = World.PortNode2Coord(QuantityId, portId, rowPortNodeId);
                int rowNodeId = World.Coord2Node(QuantityId, rowCoId);
                portEz[row] = Ez[rowNodeId];
            }
            return portEz;
        }
    }
}
