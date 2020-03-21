using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class EMWaveguide2DEigenFEM : FEM
    {
        // 磁界？
        public bool IsMagneticField { get; set; } = false;

        public IvyFEM.Lapack.DoubleMatrix Rtt { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Stz { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Szt { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Ttt { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Tzz { get; protected set; } = null;

        // Solve
        // input
        public double Frequency { get; set; }

        // output
        public System.Numerics.Complex[] Betas { get; protected set; }
        public System.Numerics.Complex[][] EtEVecs { get; protected set; }
        public System.Numerics.Complex[][] EzEVecs { get; protected set; }
        public System.Numerics.Complex[][] CoordExEyEVecs { get; protected set; }

        public EMWaveguide2DEigenFEM(FEWorld world)
        {
            World = world;
        }

        private void CalcMatrixs(double k0)
        {
            uint tQuantityId = 0;
            uint zQuantityId = 1;
            int edgeNodeCnt = (int)World.GetEdgeNodeCount(tQuantityId);
            int nodeCnt = (int)World.GetNodeCount(zQuantityId);

            Rtt = new IvyFEM.Lapack.DoubleMatrix(edgeNodeCnt, edgeNodeCnt);
            Stz = new IvyFEM.Lapack.DoubleMatrix(edgeNodeCnt, nodeCnt);
            Szt = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, edgeNodeCnt);
            Ttt = new IvyFEM.Lapack.DoubleMatrix(edgeNodeCnt, edgeNodeCnt);
            Tzz = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);

            IList<uint> tfeIds = World.GetTriangleFEIds(tQuantityId);
            foreach (uint feId in tfeIds)
            {
                // t成分
                TriangleFE tTriFE = World.GetTriangleFE(tQuantityId, feId);
                uint elemEdgeNodeCnt = tTriFE.EdgeCount;
                int[] edgeIds = new int[elemEdgeNodeCnt];
                bool[] isReverses = new bool[elemEdgeNodeCnt];
                int[] edgeNodes = new int[elemEdgeNodeCnt];
                for (int iENode = 0; iENode < elemEdgeNodeCnt; iENode++)
                {
                    int[] coIds = tTriFE.EdgeCoordIdss[iENode];
                    bool isReverse;
                    int edgeId = World.GetEdgeIdFromCoords(tQuantityId, coIds[0], coIds[1], out isReverse);
                    int edgeNodeId = World.Edge2EdgeNode(tQuantityId, edgeId);

                    edgeIds[iENode] = edgeId;
                    isReverses[iENode] = isReverse;
                    edgeNodes[iENode] = edgeNodeId;
                }
                // z成分
                TriangleFE zTriFE = World.GetTriangleFE(zQuantityId, feId);
                uint elemNodeCnt = zTriFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = zTriFE.NodeCoordIds[iNode];
                    int nodeId = World.Coord2Node(zQuantityId, coId);
                    nodes[iNode] = nodeId;
                }

                Material ma0 = World.GetMaterial(tTriFE.MaterialId);
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

                IntegrationPoints ip = TriangleFE.GetIntegrationPoints(World.TriIntegrationPointCount);//Point7
                for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                {
                    double[] L = ip.Ls[ipPt];
                    double[][] tN = tTriFE.CalcEdgeN(L);
                    double[][] tNaz = new double[elemEdgeNodeCnt][];
                    for (int iEdge = 0; iEdge < elemEdgeNodeCnt; iEdge++)
                    {
                        tNaz[iEdge] = new double[] { tN[iEdge][1], -tN[iEdge][0] }; 
                    }
                    double[] rottN = tTriFE.CalcRotEdgeN(L);
                    double[] N = zTriFE.CalcN(L);
                    double[][] Nu = zTriFE.CalcNu(L);
                    double[] Nx = Nu[0];
                    double[] Ny = Nu[1];
                    double[][] gradN = new double[elemNodeCnt][];
                    double[][] gradNaz = new double[elemNodeCnt][];
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        gradN[iNode] = new double[] { Nx[iNode], Ny[iNode] };
                        gradNaz[iNode] = new double[] { gradN[iNode][1], -gradN[iNode][0] };
                    }

                    double detJ = tTriFE.GetDetJacobian(L);
                    double weight = ip.Weights[ipPt];
                    double detJWeight = (1.0 / 2.0) * weight * detJ;

                    // Rtt
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

                            double rttVal = detJWeight * (
                                tNaz[row][0] * maPxx * tNaz[col][0] +
                                tNaz[row][1] * maPyy * tNaz[col][1]
                                );
                            Rtt[rowEdgeNodeId, colEdgeNodeId] += rowEdgeSgn * colEdgeSgn * rttVal;
                        }
                    }
                    // Stz, Szt
                    for (int row = 0; row < elemEdgeNodeCnt; row++)
                    {
                        int rowEdgeNodeId = edgeNodes[row];
                        if (rowEdgeNodeId == -1)
                        {
                            continue;
                        }
                        double rowEdgeSgn = isReverses[row] ? -1.0 : 1.0;
                        for (int col = 0; col < elemNodeCnt; col++)
                        {
                            int colNodeId = nodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double stzVal = detJWeight * (
                                tNaz[row][0] * maPxx * gradNaz[col][0] +
                                tNaz[row][1] * maPyy * gradNaz[col][1]
                                );
                            Stz[rowEdgeNodeId, colNodeId] += rowEdgeSgn * stzVal;
                            Szt[colNodeId, rowEdgeNodeId] += rowEdgeSgn * stzVal;
                        }
                    }
                    // Ttt
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

                            double tttVal = detJWeight * (
                                maPzz * rottN[row] * rottN[col] -
                                k0 * k0 * (
                                tN[row][0] * maQxx * tN[col][0] +
                                tN[row][1] * maQyy * tN[col][1])
                                );
                            Ttt[rowEdgeNodeId, colEdgeNodeId] += rowEdgeSgn * colEdgeSgn * tttVal;
                        }
                    }
                    // Tzz
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

                            double tzzVal = detJWeight * (
                                gradNaz[row][0] * maPxx * gradNaz[col][0] +
                                gradNaz[row][1] * maPyy * gradNaz[col][1] -
                                k0 * k0 * maQzz * N[row] * N[col]
                                );
                            Tzz[rowNodeId, colNodeId] += tzzVal;
                        }
                    }
                }
            }
        }

        public override void Solve()
        {
            Betas = null;
            EtEVecs = null;
            EzEVecs = null;

            uint tQuantityId = 0;
            uint zQuantityId = 1;
            // 角周波数
            double omega = 2.0 * Math.PI * Frequency;
            // 波数
            double k0 = omega / Constants.C0;

            CalcMatrixs(k0);

            int edgeNodeCnt = (int)World.GetEdgeNodeCount(tQuantityId);
            int nodeCnt = (int)World.GetNodeCount(zQuantityId);
            int sumNodeCnt = edgeNodeCnt + nodeCnt;
            int offset = edgeNodeCnt;
            var A = new IvyFEM.Lapack.DoubleMatrix(sumNodeCnt, sumNodeCnt);
            var B = new IvyFEM.Lapack.DoubleMatrix(sumNodeCnt, sumNodeCnt);

            for (int row = 0; row < edgeNodeCnt; row++)
            {
                // tt
                for (int col = 0; col < edgeNodeCnt; col++)
                {
                    A[row, col] = Ttt[row, col];
                    B[row, col] = -1.0 * Rtt[row, col]; 
                }
                // tz
                for (int col = 0; col < nodeCnt; col++)
                {
                    A[row, offset + col] = 0.0;
                    B[row, offset + col] = -1.0 * Stz[row, col];
                }
            }
            for (int row = 0; row < nodeCnt; row++)
            {
                // zt
                for (int col = 0; col < edgeNodeCnt; col++)
                {
                    A[offset + row, col] = 0.0;
                    B[offset + row, col] = -1.0 * Szt[row, col];
                }
                // zz
                for (int col = 0; col < nodeCnt; col++)
                {
                    A[offset + row, offset + col] = 0;
                    B[offset + row, offset + col] = -1.0 * Tzz[row, col];
                }
            }

            System.Numerics.Complex[] eVals = null;
            System.Numerics.Complex[][] eVecs = null;
            int ret = -1;
            try
            {
                ret = IvyFEM.Lapack.Functions.dggev_dirty(A.Buffer, A.RowLength, A.ColumnLength,
                    B.Buffer, B.RowLength, B.ColumnLength,
                    out eVals, out eVecs);
                System.Diagnostics.Debug.Assert(ret == 0);
            }
            catch (InvalidOperationException exception)
            {
                //System.Diagnostics.Debug.Assert(false);
                System.Diagnostics.Debug.WriteLine("!!!!!!!ERROR!!!!!!!!!");
                System.Diagnostics.Debug.WriteLine(exception.Message);
                System.Diagnostics.Debug.WriteLine(exception.StackTrace);
                ret = -1;
            }
            if (ret != 0)
            {
                // fail safe
                int n = A.RowLength;
                eVals = new System.Numerics.Complex[n];
                eVecs = new System.Numerics.Complex[n][];
                for (int iMode = 0; iMode < n; iMode++)
                {
                    eVecs[iMode] = new System.Numerics.Complex[n];
                }
            }

            SortEVals(eVals, eVecs);
            AdjustPhaseAndVecDirEVecs(eVecs);

            int modeCnt = eVals.Length;
            Betas = new System.Numerics.Complex[modeCnt];
            EtEVecs = new System.Numerics.Complex[modeCnt][];
            EzEVecs = new System.Numerics.Complex[modeCnt][];
            CoordExEyEVecs = new System.Numerics.Complex[modeCnt][];
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                System.Numerics.Complex beta = System.Numerics.Complex.Sqrt(eVals[iMode]);
                Betas[iMode] = beta;
                var eVec = eVecs[iMode];
                var et = new System.Numerics.Complex[edgeNodeCnt];
                var ez = new System.Numerics.Complex[nodeCnt];
                for (int edgeNodeId = 0; edgeNodeId < edgeNodeCnt; edgeNodeId++)
                {
                    et[edgeNodeId] = eVec[edgeNodeId];
                }
                for (int nodeId = 0; nodeId < nodeCnt; nodeId++)
                {
                    ez[nodeId] = beta * eVec[offset + nodeId];
                }
                EtEVecs[iMode] = et;
                EzEVecs[iMode] = ez;

                System.Numerics.Complex[] coordExEyEVec;
                CalcModeCoordExy(beta, et, out coordExEyEVec);
                CoordExEyEVecs[iMode] = coordExEyEVec;
            }
        }

        private void SortEVals(System.Numerics.Complex[] eVals, System.Numerics.Complex[][] eVecs)
        {
            int modeCnt = eVals.Length;
            var eValEVecs = new List<KeyValuePair<System.Numerics.Complex, System.Numerics.Complex[]>>();
            for (int i = 0; i < modeCnt; i++)
            {
                eValEVecs.Add(new KeyValuePair<System.Numerics.Complex, System.Numerics.Complex[]>(eVals[i], eVecs[i]));
            }
            eValEVecs.Sort((a, b) =>
            {
                // eVal(β^2) の実部を比較
                double diff = a.Key.Real - b.Key.Real;
                // 降順
                if (diff > 0)
                {
                    return -1;
                }
                else if (diff < 0)
                {
                    return 1;
                }
                return 0;
            });

            for (int i = 0; i < modeCnt; i++)
            {
                eVals[i] = eValEVecs[i].Key;
                eVecs[i] = eValEVecs[i].Value;
            }
        }

        private void AdjustPhaseAndVecDirEVecs(System.Numerics.Complex[][] eVecs)
        {
            uint tQuantityId = 0;
            uint zQuantityId = 1;
            int modeCnt = eVecs.Length;
            int edgeNodeCnt = (int)World.GetEdgeNodeCount(tQuantityId);
            int nodeCnt = (int)World.GetNodeCount(zQuantityId);
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var eVec = eVecs[iMode];
                // t成分の最大値を求める
                System.Numerics.Complex maxValue = new System.Numerics.Complex(0, 0);
                double maxAbs = 0;
                int maxEdgeId = 0;
                for (int iENode = 0; iENode < edgeNodeCnt; iENode++)
                {
                    System.Numerics.Complex value = eVec[iENode];
                    double abs = value.Magnitude;
                    if (abs > maxAbs)
                    {
                        maxAbs = abs;
                        maxValue = value;
                        maxEdgeId = World.EdgeNode2Edge(tQuantityId, iENode);
                    }
                }
                System.Numerics.Complex phase = maxValue / maxAbs;

                // t成分最大値を取る辺の方向
                double vecSgn = 1.0;
                {
                    int[] maxEdgeCoIds = World.GetEdgeCoordIds(tQuantityId, maxEdgeId);
                    int coId1 = maxEdgeCoIds[0];
                    int coId2 = maxEdgeCoIds[1];
                    double[] coord1 = World.GetCoord(tQuantityId, coId1);
                    double[] coord2 = World.GetCoord(tQuantityId, coId2);
                    double[] vec = { coord2[0] - coord1[0], coord2[1] - coord1[1] };
                    double x = vec[0];
                    double y = vec[1];
                    if (Math.Abs(x) >= Math.Abs(y))
                    {
                        vecSgn = x >= 0 ? 1.0 : -1.0; 
                    }
                    else
                    {
                        vecSgn = y >= 0 ? 1.0 : -1.0;
                    }
                }

                for (int i = 0; i < eVec.Length; i++)
                {
                    eVec[i] *= vecSgn / phase;
                }
            }
        }

        private void CalcModeCoordExy(
            System.Numerics.Complex beta, System.Numerics.Complex[] etEVec,
            out System.Numerics.Complex[] coordExEyEVec)
        {
            uint tQuantityId = 0;
            uint zQuantityId = 1;
            int edgeNodeCnt = (int)World.GetEdgeNodeCount(tQuantityId);
            int nodeCnt = (int)World.GetNodeCount(zQuantityId);
            int coCnt = (int)World.GetCoordCount(zQuantityId);

            int dof = 2; // x, y成分
            coordExEyEVec = new System.Numerics.Complex[coCnt * dof];
            int[] coordValueCnt = new int[coCnt];

            IList<uint> tfeIds = World.GetTriangleFEIds(tQuantityId);
            foreach (uint feId in tfeIds)
            {
                // t成分
                TriangleFE tTriFE = World.GetTriangleFE(tQuantityId, feId);
                uint elemEdgeNodeCnt = tTriFE.EdgeCount;
                int[] edgeIds = new int[elemEdgeNodeCnt];
                bool[] isReverses = new bool[elemEdgeNodeCnt];
                int[] edgeNodes = new int[elemEdgeNodeCnt];
                for (int iENode = 0; iENode < elemEdgeNodeCnt; iENode++)
                {
                    int[] coIds = tTriFE.EdgeCoordIdss[iENode];
                    bool isReverse;
                    int edgeId = World.GetEdgeIdFromCoords(tQuantityId, coIds[0], coIds[1], out isReverse);
                    int edgeNodeId = World.Edge2EdgeNode(tQuantityId, edgeId);

                    edgeIds[iENode] = edgeId;
                    isReverses[iENode] = isReverse;
                    edgeNodes[iENode] = edgeNodeId;
                }
                // z成分
                TriangleFE zTriFE = World.GetTriangleFE(zQuantityId, feId);
                uint elemNodeCnt = zTriFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = zTriFE.NodeCoordIds[iNode];
                    int nodeId = World.Coord2Node(zQuantityId, coId);
                    nodes[iNode] = nodeId;
                }

                // t成分の値
                System.Numerics.Complex[] et = new System.Numerics.Complex[elemEdgeNodeCnt];
                for (int iENode = 0; iENode < elemEdgeNodeCnt; iENode++)
                {
                    int edgeNodeId = edgeNodes[iENode];
                    if (edgeNodeId == -1)
                    {
                        continue;
                    }
                    double sgn = isReverses[iENode] ? -1.0 : 1.0;
                    et[iENode] = sgn * etEVec[edgeNodeId];
                }

                // z成分の節点におけるEx,Eyを求める
                System.Numerics.Complex[][] exey = new System.Numerics.Complex[elemNodeCnt][];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    double[] L = zTriFE.GetNodeL(iNode);
                    double[][] tN = tTriFE.CalcEdgeN(L);

                    exey[iNode] = new System.Numerics.Complex[dof];
                    for (int kENode = 0; kENode < elemEdgeNodeCnt; kENode++)
                    {
                        for (int idim = 0; idim < dof; idim++)
                        {
                            exey[iNode][idim] += tN[kENode][idim] * et[kENode];
                        }
                    }
                }

                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    //int nodeId = nodes[iNode];
                    int coId = zTriFE.NodeCoordIds[iNode];
                    for (int idim = 0; idim < dof; idim++)
                    {
                        coordExEyEVec[coId * dof + idim] += exey[iNode][idim];
                    }
                    coordValueCnt[coId]++;
                }
            }

            for (int coId = 0; coId < coCnt; coId++)
            {
                for (int idim = 0; idim < dof; idim++)
                {
                    coordExEyEVec[coId * dof + idim] /= coordValueCnt[coId];
                }
            }
        }
    }
}
