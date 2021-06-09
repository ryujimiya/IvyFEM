using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class EMWaveguide3DFEM : FEM
    {
        // 磁界？
        public bool IsMagneticField { get; set; } = false;

        // Solve
        // input
        public double Frequency { get; set; }
        // output
        public System.Numerics.Complex[] E { get; private set; }
        public System.Numerics.Complex[] CoordExyz { get; private set; }
        public System.Numerics.Complex[][] S { get; private set; }
        public EMWaveguide2DEigenFEMForPort[] EigenFEMs { get; private set; } 

        public EMWaveguide3DFEM(FEWorld world)
        {
            World = world;
        }

        public override void Solve()
        {
            int t;
            E = null;
            S = null;
            EigenFEMs = null;

            // 角周波数
            double omega = 2.0 * Math.PI * Frequency;
            // 波数
            double k0 = omega / Constants.C0;

            uint quantityId = 0;

            t = System.Environment.TickCount;
            int edgeNodeCnt = (int)World.GetEdgeNodeCount(quantityId);
            var A = new IvyFEM.Linear.ComplexSparseMatrix(edgeNodeCnt, edgeNodeCnt);
            var B = new System.Numerics.Complex[edgeNodeCnt];
            CalcA(k0, A);
            System.Diagnostics.Debug.WriteLine("CalcA t = " + (System.Environment.TickCount - t));

            t = System.Environment.TickCount;
            uint portCnt = World.GetPortCount(quantityId);
            EMWaveguide2DEigenFEMForPort[] eigenFEMs;
            SetBoundaryCondition(omega, A, B, portCnt, out eigenFEMs);
            System.Diagnostics.Debug.WriteLine("SetBoundaryCondition t = " + (System.Environment.TickCount - t));
            EigenFEMs = eigenFEMs;

            t = System.Environment.TickCount;
            //----------------------------------
            System.Numerics.Complex[] X;
            Solver.ComplexSolve(out X, A, B);
            E = X;
            //----------------------------------
            System.Diagnostics.Debug.WriteLine("Solve t = " + (System.Environment.TickCount - t));

            System.Numerics.Complex[] coordExyz;
            CalcModeCoordExyz(E, out coordExyz);
            CoordExyz = coordExyz;

            t = System.Environment.TickCount;
            S = CalcS(omega, E, portCnt, eigenFEMs);
            System.Diagnostics.Debug.WriteLine("CalcS t = " + (System.Environment.TickCount - t));
        }

        private void CalcA(double k0, IvyFEM.Linear.ComplexSparseMatrix A)
        {
            uint quantityId = 0;

            IList<uint> feIds = World.GetTetrahedronFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = World.GetTetrahedronFE(quantityId, feId);
                uint elemNodeCnt = tetFE.NodeCount;
                uint elemEdgeNodeCnt = tetFE.EdgeCount;
                int[] edgeIds = new int[elemEdgeNodeCnt];
                bool[] isReverses = new bool[elemEdgeNodeCnt];
                int[] edgeNodes = new int[elemEdgeNodeCnt];
                for (int iENode = 0; iENode < elemEdgeNodeCnt; iENode++)
                {
                    int[] coIds = tetFE.EdgeCoordIdss[iENode];
                    bool isReverse;
                    int edgeId = World.GetEdgeIdFromCoords(quantityId, coIds[0], coIds[1], out isReverse);
                    int edgeNodeId = World.Edge2EdgeNode(quantityId, edgeId);

                    edgeIds[iENode] = edgeId;
                    isReverses[iENode] = isReverse;
                    edgeNodes[iENode] = edgeNodeId;
                }

                Material ma0 = World.GetMaterial(tetFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is DielectricMaterial);
                var ma = ma0 as DielectricMaterial;
                double maPxx = 0.0;
                double maPyy = 0.0;
                double maPzz = 0.0;
                double maQxx = 0.0;
                double maQyy = 0.0;
                double maQzz = 0.0;
                if (IsMagneticField)
                {
                    // 磁界
                    maPxx = 1.0 / ma.Epxx;
                    maPyy = 1.0 / ma.Epyy;
                    maPzz = 1.0 / ma.Epzz;
                    maQxx = ma.Muxx;
                    maQyy = ma.Muyy;
                    maQzz = ma.Muzz;
                }
                else
                {
                    // 電界
                    maPxx = 1.0 / ma.Muxx;
                    maPyy = 1.0 / ma.Muyy;
                    maPzz = 1.0 / ma.Muzz;
                    maQxx = ma.Epxx;
                    maQyy = ma.Epyy;
                    maQzz = ma.Epzz;
                }

                IntegrationPoints ip = TetrahedronFE.GetIntegrationPoints(World.TetIntegrationPointCount);//Point5
                for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                {
                    double[] L = ip.Ls[ipPt];
                    double[][] N = tetFE.CalcEdgeN(L);
                    double[][] rotN = tetFE.CalcRotEdgeN(L);

                    double detJ = tetFE.GetDetJacobian(L);
                    double weight = ip.Weights[ipPt];
                    double detJWeight = (1.0 / 6.0) * weight * detJ;

                    for (int row = 0; row < elemEdgeNodeCnt; row++)
                    {
                        int rowEdgeNodeId = edgeNodes[row];
                        if (rowEdgeNodeId == -1)
                        {
                            continue;
                        }
                        double rowEdgeSgn = isReverses[row] ? -1.0 : 1.0;
                        for (int col = 0; col < elemEdgeNodeCnt; col++)
                        {
                            int colEdgeNodeId = edgeNodes[col];
                            if (colEdgeNodeId == -1)
                            {
                                continue;
                            }
                            double colEdgeSgn = isReverses[col] ? -1.0 : 1.0;

                            double kVal = detJWeight * (
                                rotN[row][0] * maPxx * rotN[col][0] +
                                rotN[row][1] * maPyy * rotN[col][1] +
                                rotN[row][2] * maPzz * rotN[col][2]);
                            double mVal = detJWeight * (
                                N[row][0] * maQxx * N[col][0] +
                                N[row][1] * maQyy * N[col][1] +
                                N[row][2] * maQzz * N[col][2]);
                            A[rowEdgeNodeId, colEdgeNodeId] += 
                                rowEdgeSgn * colEdgeSgn * (kVal - k0 * k0 * mVal);
                        }
                    }
                }
            }
        }

        private void SetBoundaryCondition(double omega,
            IvyFEM.Linear.ComplexSparseMatrix A, System.Numerics.Complex[] B,
            uint portCnt, out EMWaveguide2DEigenFEMForPort[] eigenFEMs)
        {
            uint quantityId = 0;

            eigenFEMs = new EMWaveguide2DEigenFEMForPort[portCnt];
            for (uint portId = 0; portId < portCnt; portId++)
            {
                uint portEdgeNodeCnt = World.GetPortEdgeNodeCount(quantityId, portId);

                var eigenFEM = new EMWaveguide2DEigenFEMForPort(World, portId);
                eigenFEMs[portId] = eigenFEM;

                eigenFEM.IsMagneticField = IsMagneticField;
                eigenFEM.Frequency = Frequency;
                eigenFEM.Solve();
                System.Numerics.Complex[] betas = eigenFEM.Betas;
                System.Numerics.Complex[][] etEVecs = eigenFEM.EtEVecs;
                System.Numerics.Complex[][] ezEVecs = eigenFEM.EzEVecs;
                IvyFEM.Lapack.ComplexMatrix b = eigenFEM.CalcBoundaryMatrix(omega, betas, etEVecs, ezEVecs);
                for (int row = 0; row < portEdgeNodeCnt; row++)
                {
                    int rowEdgeId = World.PortEdgeNode2Edge(quantityId, portId, row);
                    int rowNodeId = World.Edge2EdgeNode(quantityId, rowEdgeId);
                    for (int col = 0; col < portEdgeNodeCnt; col++)
                    {
                        int colEdgeId = World.PortEdgeNode2Edge(quantityId, portId, col);
                        int colNodeId = World.Edge2EdgeNode(quantityId, colEdgeId);

                        A[rowNodeId, colNodeId] += -1.0 * b[row, col];
                    }
                }

                bool isIncidentPort = (portId == World.GetIncidentPortId(quantityId));
                if (isIncidentPort)
                {
                    int incidentModeId = World.GetIncidentModeId(quantityId);
                    System.Diagnostics.Debug.Assert(incidentModeId != -1);
                    System.Numerics.Complex beta0 = betas[incidentModeId];
                    System.Numerics.Complex[] etEVec0 = etEVecs[incidentModeId];
                    System.Numerics.Complex[] ezEVec0 = ezEVecs[incidentModeId];
                    System.Numerics.Complex[] I = eigenFEM.CalcIncidentVec(omega, beta0, etEVec0, ezEVec0);
                    for (int row = 0; row < portEdgeNodeCnt; row++)
                    {
                        int rowEdgeId = World.PortEdgeNode2Edge(quantityId, portId, row);
                        int rowNodeId = World.Edge2EdgeNode(quantityId, rowEdgeId);

                        B[rowNodeId] += I[row];
                    }
                }
            }
        }

        private System.Numerics.Complex[][] CalcS(double omega,
            System.Numerics.Complex[] E,
            uint portCnt, EMWaveguide2DEigenFEMForPort[] eigenFEMs)
        {
            uint quantityId = 0;

            var S = new System.Numerics.Complex[portCnt][];
            for (uint portId = 0; portId < portCnt; portId++)
            {
                System.Numerics.Complex[] portEt = GetPortEt(portId, E);
                var eigenFEM = eigenFEMs[portId];
                System.Numerics.Complex[] betas = eigenFEM.Betas;
                System.Numerics.Complex[][] etEVecs = eigenFEM.EtEVecs;
                System.Numerics.Complex[][] ezEVecs = eigenFEM.EzEVecs;
                int incidentModeId = -1;
                if (World.GetIncidentPortId(quantityId) == portId)
                {
                    incidentModeId = World.GetIncidentModeId(quantityId);
                }
                System.Numerics.Complex[] S1 = eigenFEM.CalcSMatrix(
                    omega, incidentModeId, betas, etEVecs, ezEVecs, portEt);
                S[portId] = S1;
            }
            return S;
        }

        private System.Numerics.Complex[] GetPortEt(uint portId, System.Numerics.Complex[] E)
        {
            uint tQuantityId = 0;
            uint zQuantityId = 1;
            int edgeNodeCnt = (int)World.GetPortEdgeNodeCount(tQuantityId, portId);
            int nodeCnt = (int)World.GetPortNodeCount(zQuantityId, portId);
            int sumNodeCnt = edgeNodeCnt + nodeCnt;
            int offset = edgeNodeCnt;

            System.Numerics.Complex[] portEt = new System.Numerics.Complex[edgeNodeCnt];
            for (int row = 0; row < edgeNodeCnt; row++)
            {
                int edgeId = World.PortEdgeNode2Edge(tQuantityId, portId, row);
                int nodeId = World.Edge2EdgeNode(tQuantityId, edgeId);
                portEt[row] = E[nodeId];
            }
            return portEt;
        }

        private void CalcModeCoordExyz(
            System.Numerics.Complex[] E,
            out System.Numerics.Complex[] coordExyz)
        {
            uint quantityId = 0;
            int edgeNodeCnt = (int)World.GetEdgeNodeCount(quantityId);
            int nodeCnt = (int)World.GetNodeCount(quantityId);
            int coCnt = (int)World.GetCoordCount(quantityId);

            int dof = 3; // x, y, z成分
            coordExyz = new System.Numerics.Complex[coCnt * dof];
            int[] coordValueCnt = new int[coCnt];

            IList<uint> feIds = World.GetTetrahedronFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = World.GetTetrahedronFE(quantityId, feId);
                uint elemEdgeNodeCnt = tetFE.EdgeCount;
                uint elemNodeCnt = tetFE.NodeCount;
                int[] edgeIds = new int[elemEdgeNodeCnt];
                bool[] isReverses = new bool[elemEdgeNodeCnt];
                int[] edgeNodes = new int[elemEdgeNodeCnt];
                for (int iENode = 0; iENode < elemEdgeNodeCnt; iENode++)
                {
                    int[] coIds = tetFE.EdgeCoordIdss[iENode];
                    bool isReverse;
                    int edgeId = World.GetEdgeIdFromCoords(quantityId, coIds[0], coIds[1], out isReverse);
                    int edgeNodeId = World.Edge2EdgeNode(quantityId, edgeId);

                    edgeIds[iENode] = edgeId;
                    isReverses[iENode] = isReverse;
                    edgeNodes[iENode] = edgeNodeId;
                }

                // 辺方向成分の値
                System.Numerics.Complex[] et = new System.Numerics.Complex[elemEdgeNodeCnt];
                for (int iENode = 0; iENode < elemEdgeNodeCnt; iENode++)
                {
                    int edgeNodeId = edgeNodes[iENode];
                    if (edgeNodeId == -1)
                    {
                        continue;
                    }
                    double sgn = isReverses[iENode] ? -1.0 : 1.0;
                    et[iENode] = sgn * E[edgeNodeId];
                }

                // 節点のL
                double[][] nodeL = new double[4][]
                {
                    new double[4] { 1.0, 0.0, 0.0, 0.0 },
                    new double[4] { 0.0, 1.0, 0.0, 0.0 },
                    new double[4] { 0.0, 0.0, 1.0, 0.0 },
                    new double[4] { 0.0, 0.0, 0.0, 1.0 }
                };
                System.Diagnostics.Debug.Assert(nodeL.Length == elemNodeCnt);
                // 節点におけるEx,Ey,Ezを求める
                System.Numerics.Complex[][] exyz = new System.Numerics.Complex[elemNodeCnt][];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    double[] L = nodeL[iNode];
                    double[][] N = tetFE.CalcEdgeN(L);

                    exyz[iNode] = new System.Numerics.Complex[dof];
                    for (int kENode = 0; kENode < elemEdgeNodeCnt; kENode++)
                    {
                        for (int idim = 0; idim < dof; idim++)
                        {
                            exyz[iNode][idim] += N[kENode][idim] * et[kENode];
                        }
                    }
                }

                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    //int nodeId = nodes[iNode];
                    int coId = tetFE.NodeCoordIds[iNode];
                    for (int idim = 0; idim < dof; idim++)
                    {
                        coordExyz[coId * dof + idim] += exyz[iNode][idim];
                    }
                    coordValueCnt[coId]++;
                }
            }

            for (int coId = 0; coId < coCnt; coId++)
            {
                for (int idim = 0; idim < dof; idim++)
                {
                    coordExyz[coId * dof + idim] /= coordValueCnt[coId];
                }
            }
        }
    }
}
