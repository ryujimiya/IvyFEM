using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class EMWaveguide2DHPlaneFEM : FEM
    {
        public uint QuantityId { get; private set; } = 0;
        public EMWaveguideType WaveguideType { get; set; } = EMWaveguideType.HPlane2D;
        /// <summary>
        ///  E面導波路の場合の導波管幅a(図面で言うと高さ方向の幅)
        /// </summary>
        public double WaveguideWidthForEPlane { get; set; } = 0;
        /// <summary>
        /// TEモードで実装した式をTMモードに流用するため
        ///   TEモードの場合は μ0
        ///   TMモードの場合は ε0
        /// </summary>
        public double ReplacedMu0 = Constants.Mu0;

        // Solve
        // input
        public double Frequency { get; set; }
        // output
        public System.Numerics.Complex[] Ez { get; private set; }
        public System.Numerics.Complex[][] S { get; private set; }
        public EMWaveguide1DEigenFEM[] EigenFEMs { get; private set; } 

        public EMWaveguide2DHPlaneFEM(FEWorld world)
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
            System.Diagnostics.Debug.WriteLine("CalcA t = " + (System.Environment.TickCount - t));

            t = System.Environment.TickCount;
            uint portCnt = World.GetPortCount(QuantityId);
            EMWaveguide1DEigenFEM[] eigenFEMs;
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

        private void CalcA(double k0, IvyFEM.Linear.ComplexSparseMatrix A)
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
                if (WaveguideType == EMWaveguideType.HPlane2D)
                {
                    maPxx = 1.0 / ma.Muxx;
                    maPyy = 1.0 / ma.Muyy;
                    maQzz = ma.Epzz;
                }
                else if (WaveguideType == EMWaveguideType.EPlane2D)
                {
                    // LSE(TE^z)モード(Ez = 0:紙面に垂直な方向の電界を０)として解析する
                    //   波動方程式の導出でμx = μy  εx = εyを仮定した
                    maPxx = 1.0 / ma.Epxx;
                    maPyy = 1.0 / ma.Epyy;
                    maQzz = ma.Muzz -
                        (Math.PI * Math.PI * ma.Muzz) /
                        (k0 * k0 * ma.Epyy * WaveguideWidthForEPlane * WaveguideWidthForEPlane * ma.Muxx);
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
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
                        double a = maPxx * sNyNy[row, col] + maPyy * sNxNx[row, col] -
                            (k0 * k0 * maQzz) * sNN[row, col];

                        A[rowNodeId, colNodeId] += (System.Numerics.Complex)a;
                    }
                }
            }
        }

        private void SetBoundaryCondition(double omega,
            IvyFEM.Linear.ComplexSparseMatrix A, System.Numerics.Complex[] B,
            uint portCnt, out EMWaveguide1DEigenFEM[] eigenFEMs)
        {
            eigenFEMs = new EMWaveguide1DEigenFEM[portCnt];
            for (uint portId = 0; portId < portCnt; portId++)
            {
                uint portNodeCnt = World.GetPortNodeCount(QuantityId, portId);

                var eigenFEM = new EMWaveguide1DEigenFEM(World, QuantityId, portId);
                eigenFEMs[portId] = eigenFEM;

                eigenFEM.WaveguideType = WaveguideType;
                eigenFEM.WaveguideWidthForEPlane = WaveguideWidthForEPlane;
                eigenFEM.ReplacedMu0 = ReplacedMu0;
                eigenFEM.Frequency = Frequency;
                eigenFEM.Solve();
                System.Numerics.Complex[] betas = eigenFEM.Betas;
                System.Numerics.Complex[][] ezEVecs = eigenFEM.EzEVecs;
                IvyFEM.Lapack.ComplexMatrix b = eigenFEM.CalcBoundaryMatrix(omega, betas, ezEVecs);
                for (int row = 0; row < portNodeCnt; row++)
                {
                    int rowCoId = World.PortNode2Coord(QuantityId, portId, row);
                    int rowNodeId = World.Coord2Node(QuantityId, rowCoId);
                    for (int col = 0; col < portNodeCnt; col++)
                    {
                        int colCoId = World.PortNode2Coord(QuantityId, portId, col);
                        int colNodeId = World.Coord2Node(QuantityId, colCoId);

                        A[rowNodeId, colNodeId] += b[row, col];
                    }
                }

                bool isIncidentPort = (portId == World.GetIncidentPortId(QuantityId));
                if (isIncidentPort)
                {
                    int incidentModeId = World.GetIncidentModeId(QuantityId);
                    System.Diagnostics.Debug.Assert(incidentModeId != -1);
                    System.Numerics.Complex beta0 = betas[incidentModeId];
                    System.Numerics.Complex[] ezEVec0 = ezEVecs[incidentModeId];
                    System.Numerics.Complex[] I = eigenFEM.CalcIncidentVec(beta0, ezEVec0);
                    for (int row = 0; row < portNodeCnt; row++)
                    {
                        int rowCoId = World.PortNode2Coord(QuantityId, portId, row);
                        int rowNodeId = World.Coord2Node(QuantityId, rowCoId);

                        B[rowNodeId] += I[row];
                    }
                }
            }
        }

        private System.Numerics.Complex[][] CalcS(double omega,
            System.Numerics.Complex[] Ez,
            uint portCnt, EMWaveguide1DEigenFEM[] eigenFEMs)
        {
            var S = new System.Numerics.Complex[portCnt][];
            for (uint portId = 0; portId < portCnt; portId++)
            {
                System.Numerics.Complex[] portEz = GetPortEz(portId, Ez);
                var eigenFEM = eigenFEMs[portId];
                System.Numerics.Complex[] betas = eigenFEM.Betas;
                System.Numerics.Complex[][] ezEVecs = eigenFEM.EzEVecs;
                int incidentModeId = -1;
                if (World.GetIncidentPortId(QuantityId) == portId)
                {
                    incidentModeId = (int)World.GetIncidentModeId(QuantityId);
                }
                System.Numerics.Complex[] S1 = eigenFEM.CalcSMatrix(omega, incidentModeId, betas, ezEVecs, portEz);
                S[portId] = S1;
            }
            return S;
        }

        private System.Numerics.Complex[] GetPortEz(uint portId, System.Numerics.Complex[] Ez)
        {
            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, portId);
            System.Numerics.Complex[] portEz= new System.Numerics.Complex[nodeCnt];
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
