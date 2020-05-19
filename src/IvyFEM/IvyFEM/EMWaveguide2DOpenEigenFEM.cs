using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class EMWaveguide2DOpenEigenFEM : FEM
    {
        // 磁界？
        public bool IsMagneticField { get; set; } = false;

        public IvyFEM.Lapack.DoubleMatrix Rtt { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Stz { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Szt { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Ttt { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Tzz { get; protected set; } = null;
        // for Open
        protected IvyFEM.Lapack.DoubleMatrix Utt = null;
        protected IvyFEM.Lapack.DoubleMatrix Utn = null;
        protected IvyFEM.Lapack.DoubleMatrix Uzn = null;
        protected IvyFEM.Lapack.DoubleMatrix Uzz = null;

        /// <summary>
        /// 最大比誘電率(フィルタリングに用いる)
        /// </summary>
        public double MaxEp { get; set; } = double.MaxValue;
        /// <summary>
        /// クラッドの比誘電率
        /// 外部の比誘電率
        /// </summary>
        public IList<double> PortEp { get; set; } = new List<double>();
        /// <summary>
        /// 減衰定数
        /// </summary>
        public IList<double> DecayParameters { get; set; } = new List<double>();

        private System.Numerics.Complex PrevBeta = 0.0;
        private System.Numerics.Complex[] Prev2Betas = new System.Numerics.Complex[10];
        private System.Numerics.Complex[] PrevEVec = null;

        // Solve
        // input
        public double Frequency { get; set; }

        // output
        public System.Numerics.Complex[] Betas { get; protected set; }
        public System.Numerics.Complex[][] EtEVecs { get; protected set; }
        public System.Numerics.Complex[][] EzEVecs { get; protected set; }
        public System.Numerics.Complex[][] CoordExEyEVecs { get; protected set; }

        public EMWaveguide2DOpenEigenFEM(FEWorld world)
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

        private void CalcOpenBoundary(
            double k0, double[] alphas)
        {
            uint tQuantityId = 0;
            uint zQuantityId = 1;
            int edgeNodeCnt = (int)World.GetEdgeNodeCount(tQuantityId);
            int nodeCnt = (int)World.GetNodeCount(zQuantityId);

            Utt = new IvyFEM.Lapack.DoubleMatrix(edgeNodeCnt, edgeNodeCnt);
            Utn = new IvyFEM.Lapack.DoubleMatrix(edgeNodeCnt, edgeNodeCnt);
            Uzn = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, edgeNodeCnt);
            Uzz = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);

            int portCnt = (int)World.GetPortCount(tQuantityId);
            for (uint portId = 0; portId < portCnt; portId++)
            {
                double alpha = alphas[portId];

                IList<uint> feIds = World.GetPortLineFEIds(tQuantityId, portId);
                foreach (uint feId in feIds)
                {
                    // t成分
                    LineFE tLineFE = World.GetLineFE(tQuantityId, feId);
                    uint elemEdgeNodeCnt = tLineFE.EdgeCount;
                    int[] edgeIds = new int[elemEdgeNodeCnt];
                    bool[] isReverses = new bool[elemEdgeNodeCnt];
                    int[] edgeNodes = new int[elemEdgeNodeCnt];
                    for (int iENode = 0; iENode < elemEdgeNodeCnt; iENode++)
                    {
                        int[] coIds = tLineFE.EdgeCoordIdss[iENode];
                        bool isReverse;
                        int edgeId = World.GetEdgeIdFromCoords(tQuantityId, coIds[0], coIds[1], out isReverse);
                        int edgeNodeId = World.Edge2EdgeNode(tQuantityId, edgeId);

                        edgeIds[iENode] = edgeId;
                        isReverses[iENode] = isReverse;
                        edgeNodes[iENode] = edgeNodeId;
                    }
                    // z成分
                    LineFE zLineFE = World.GetLineFE(zQuantityId, feId);
                    uint elemNodeCnt = zLineFE.NodeCount;
                    int[] nodes = new int[elemNodeCnt];
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int coId = zLineFE.NodeCoordIds[iNode];
                        int nodeId = World.Coord2Node(zQuantityId, coId);
                        nodes[iNode] = nodeId;
                    }
                    // 境界に接する三角形
                    TriangleFE tTriFE;
                    {
                        int coId1 = tLineFE.EdgeCoordIdss[0][0];
                        int coId2 = tLineFE.EdgeCoordIdss[0][1];
                        IList<uint> triFEIds = World.GetTriangleFEIdsFromEdgeCoord(tQuantityId, coId1, coId2);
                        System.Diagnostics.Debug.Assert(triFEIds.Count == 1);
                        uint triFEId = triFEIds[0];
                        tTriFE = World.GetTriangleFE(tQuantityId, triFEId);
                    }
                    uint tTriElemEdgeNodeCnt = tTriFE.EdgeCount;
                    int[] tTriEdgeIds = new int[tTriElemEdgeNodeCnt];
                    bool[] tTriIsReverses = new bool[tTriElemEdgeNodeCnt];
                    int[] tTriEdgeNodes = new int[tTriElemEdgeNodeCnt];
                    for (int iENode = 0; iENode < tTriElemEdgeNodeCnt; iENode++)
                    {
                        int[] coIds = tTriFE.EdgeCoordIdss[iENode];
                        bool isReverse;
                        int edgeId = World.GetEdgeIdFromCoords(tQuantityId, coIds[0], coIds[1], out isReverse);
                        int edgeNodeId = World.Edge2EdgeNode(tQuantityId, edgeId);

                        tTriEdgeIds[iENode] = edgeId;
                        tTriIsReverses[iENode] = isReverse;
                        tTriEdgeNodes[iENode] = edgeNodeId;
                    }

                    Material ma0 = World.GetMaterial(tLineFE.MaterialId);
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

                    double le = zLineFE.GetLineLength();
                    double[] normal = zLineFE.GetNormal();
                    bool isYDirection = true;
                    double normalSgn = 1.0;
                    if (Math.Abs(normal[0]) >= 1.0e-12 && Math.Abs(normal[1]) < 1.0e-12)
                    {
                        // n = ax
                        isYDirection = true;
                        normalSgn = normal[0] >= 0.0 ? 1.0 : -1.0;
                    }
                    else if (Math.Abs(normal[0]) < 1.0e-12 && Math.Abs(normal[1]) >= 1.0e-12)
                    {
                        // n = ay
                        isYDirection = false;
                        normalSgn = normal[1] >= 0.0 ? 1.0 : -1.0;
                    }
                    else
                    {
                        // Not Supported
                        System.Diagnostics.Debug.Assert(false);
                    }
                    IntegrationPoints ip = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);//Point5
                    for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                    {
                        double[] L = ip.Ls[ipPt];
                        double[][] tN = tLineFE.CalcEdgeN(L);
                        double[] N = zLineFE.CalcN(L);
                        double[][] Nu = zLineFE.CalcNu(L);
                        double[] Ns = Nu[0];
                        double[] triL = new double[3] { 0.0, 0.0, 0.0 };
                        {
                            int vertexCo1 = zLineFE.VertexCoordIds[0];
                            int vertexCo2 = zLineFE.VertexCoordIds[1];
                            int iL1 = Array.IndexOf(tTriFE.VertexCoordIds, vertexCo1);
                            int iL2 = Array.IndexOf(tTriFE.VertexCoordIds, vertexCo2);
                            triL[iL1] = L[0];
                            triL[iL2] = L[1];
                        }
                        double[][] tTriN = tTriFE.CalcEdgeN(triL);
                        double[][][] tTriNu = tTriFE.CalcEdgeNu(triL);
                        double[][] tTriNx = tTriNu[0];
                        double[][] tTriNy = tTriNu[1];

                        double weight = ip.Weights[ipPt];
                        double detJWeight = (1.0 / 2.0) * weight * le;

                        // utt
                        for (int row = 0; row < elemEdgeNodeCnt; row++)
                        {
                            int rowEdgeNodeId = edgeNodes[row];
                            bool rowIsReverse = isReverses[row];
                            double rowEdgeSgn = rowIsReverse ? -1.0 : 1.0;
                            if (rowEdgeNodeId == -1)
                            {
                                continue;
                            }
                            for (int col = 0; col < elemEdgeNodeCnt; col++)
                            {
                                int colEdgeNodeId = edgeNodes[col];
                                bool colIsReverse = isReverses[col];
                                double colEdgeSgn = colIsReverse ? -1.0 : 1.0;

                                if (colEdgeNodeId == -1)
                                {
                                    continue;
                                }
                                double uttVal = 0.0;
                                if (isYDirection)
                                {
                                    uttVal = detJWeight * alpha * maPzz * tN[row][1] * tN[col][1];
                                }
                                else
                                {
                                    uttVal = detJWeight * alpha * maPzz * tN[row][0] * tN[col][0];
                                }
                                Utt[rowEdgeNodeId, colEdgeNodeId] += rowEdgeSgn * colEdgeSgn * uttVal;
                            }
                        }

                        // utn
                        for (int row = 0; row < elemEdgeNodeCnt; row++)
                        {
                            int rowEdgeNodeId = edgeNodes[row];
                            bool rowIsReverse = isReverses[row];
                            double rowEdgeSgn = rowIsReverse ? -1.0 : 1.0;
                            if (rowEdgeNodeId == -1)
                            {
                                continue;
                            }

                            /* En = 0とすると良好な結果が得られた
                            /////////////////////////////////////////////////////
                            // n成分
                            for (int kENode = 0; kENode < tTriElemEdgeNodeCnt; kENode++)
                            {
                                int kEdgeNodeId = tTriEdgeNodes[kENode];
                                bool kIsReverse = tTriIsReverses[kENode];
                                double kEdgeSgn = kIsReverse ? -1.0 : 1.0;
                                if (kEdgeNodeId == -1)
                                {
                                    continue;
                                }
                                double utnVal = 0;
                                if (isYDirection)
                                {
                                    utnVal = detJWeight * maPzz * tN[row][1] * normalSgn * tTriNy[kENode][0];
                                }
                                else
                                {
                                    utnVal = detJWeight * maPzz * tN[row][0] * normalSgn * tTriNx[kENode][1];
                                }
                                Utn[rowEdgeNodeId, kEdgeNodeId] += rowEdgeSgn * kEdgeSgn * utnVal;
                            }
                            /////////////////////////////////////////////////////
                            */
                        }

                        // uzn
                        for (int row = 0; row < elemNodeCnt; row++)
                        {
                            int rowNodeId = nodes[row];
                            if (rowNodeId == -1)
                            {
                                continue;
                            }

                            /* En = 0とすると良好な結果が得られた
                            /////////////////////////////////////////////////////
                            // n成分
                            for (int kENode = 0; kENode < tTriElemEdgeNodeCnt; kENode++)
                            {
                                int kEdgeNodeId = tTriEdgeNodes[kENode];
                                bool kIsReverse = tTriIsReverses[kENode];
                                double kEdgeSgn = kIsReverse ? -1.0 : 1.0;
                                if (kEdgeNodeId == -1)
                                {
                                    continue;
                                }
                                double uznVal = 0;
                                if (isYDirection)
                                {
                                    uznVal = -1.0 * detJWeight * maPyy * N[row] * normalSgn * tTriN[kENode][0];
                                }
                                else
                                {
                                    uznVal = -1.0 * detJWeight * maPxx * N[row] * normalSgn * tTriN[kENode][1];
                                }
                                Uzn[rowNodeId, kEdgeNodeId] += kEdgeSgn * uznVal;
                            }
                            /////////////////////////////////////////////////////
                            */
                        }
                        // uzz
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

                                double uzzVal = 0;
                                if (isYDirection)
                                {
                                    uzzVal = detJWeight * alpha * maPyy * N[row] * N[col];
                                }
                                else
                                {
                                    uzzVal = detJWeight * alpha * maPxx * N[row] * N[col];
                                }
                                Uzz[rowNodeId, colNodeId] += uzzVal;
                            }
                        }
                    }
                }
            }
        }

        public override void Solve()
        {
            // 角周波数
            double omega = 2.0 * Math.PI * Frequency;
            // 波数
            double k0 = omega / Constants.C0;

            uint tQuantityId = 0;
            uint zQuantityId = 1;

            //int iterCnt = 1000;
            int iterCnt = 50;
            int iter = 0;
            int portCnt = (int)World.GetPortCount(tQuantityId);
            System.Diagnostics.Debug.Assert(PortEp.Count == portCnt);
            //double claddingEp0 = double.MaxValue;
            //{
            //    for (uint portId = 0; portId < portCnt; portId++)
            //    {
            //        if (PortEp[(int)portId] < claddingEp0)
            //        {
            //            claddingEp0 = PortEp[(int)portId];
            //        }
            //    }
            //}
            //System.Numerics.Complex prevBeta = k0 * Math.Sqrt(claddingEp0);
            //System.Numerics.Complex prevBeta = k0;
            System.Numerics.Complex prevBeta = 0.0;
            double convRatio = 0;
            bool success = false;
            PrevEVec = null;
            for (iter = 0; iter < iterCnt; iter++)
            {
                DecayParameters = new List<double>();
                for (uint portId = 0; portId < portCnt; portId++)
                {
                    double claddingEp = PortEp[(int)portId];
                    System.Numerics.Complex alphaZ2 = prevBeta * prevBeta - k0 * k0 * claddingEp;
                    System.Numerics.Complex alphaZ = System.Numerics.Complex.Sqrt(alphaZ2);
                    double alpha = alphaZ.Real;
                    if (alpha < 0.0)
                    {
                        alpha = 0.0;
                    }
                    DecayParameters.Add(alpha);
                }
                for (int i = 0; i < (Prev2Betas.Length - 1); i++)
                {
                    Prev2Betas[i + 1] = Prev2Betas[i];
                }
                Prev2Betas[0] = PrevBeta;
                PrevBeta = prevBeta;

                success = SolveIter();
                if (!success)
                {
                    // 失敗したらすぐ抜ける
                    break;
                }
                double beta = Betas[0].Real;
                if (Math.Abs(beta) < 1.0e-30)
                {
                    System.Diagnostics.Debug.WriteLine("[ERROR] curIter={0} β={1} exit loop", iter, beta);
                    break;
                }
                convRatio = ((beta - prevBeta) / beta).Magnitude;
                System.Diagnostics.Debug.WriteLine(
                    "curIter= {0} curConvRatio = {1} β/k0 = {2}", iter, convRatio, beta / k0);
                if (convRatio < 1.0e-6)
                {
                    break;
                }
                // 振動するケース:2つ前以前と比較
                {
                    bool isConverge2 = false; 
                    for (int i = 0; i < Prev2Betas.Length; i++)
                    {
                        double convRatio2 = ((beta - Prev2Betas[i]) / beta).Magnitude;
                        if (convRatio2 < 1.0e-6)
                        {
                            System.Diagnostics.Debug.WriteLine(
                                "curIter= {0} curConvRatio2 = {1} β/k0 = {2}", iter, convRatio2, beta / k0);
                            convRatio = convRatio2;
                            isConverge2 = true;
                            break;
                        }
                    }
                    if (isConverge2)
                    {
                        break;
                    }
                }

                prevBeta = beta;
            }
            //System.Diagnostics.Debug.Assert(iter < iterCnt);
            if (iter == iterCnt)
            {
                System.Diagnostics.Debug.WriteLine(
                    "EMWaveguide2DOpenEigenFEM *Not* converge: Frequency = {0} convRatio = {1}", Frequency, convRatio);
            }
            if (!success)
            {
                System.Diagnostics.Debug.WriteLine(
                    "EMWaveguide2DOpenEigenFEM failed: Frequency = {0}", Frequency);
            }
        }

        ////////////////////////////////////////
        private bool SolveIter()
        {
            EtEVecs = null;
            EzEVecs = null;

            bool success = false;

            uint tQuantityId = 0;
            uint zQuantityId = 1;
            // 角周波数
            double omega = 2.0 * Math.PI * Frequency;
            // 波数
            double k0 = omega / Constants.C0;

            CalcMatrixs(k0);

            {
                double[] alphas = DecayParameters.ToArray();
                CalcOpenBoundary(k0, alphas);
            }

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
                    A[row, col] = Ttt[row, col] + Utt[row, col] + Utn[row, col];
                    B[row, col] = -Rtt[row, col];
                }
                // tz
                for (int col = 0; col < nodeCnt; col++)
                {
                    A[row, offset + col] = 0.0;
                    B[row, offset + col] = -Stz[row, col];
                }
            }
            for (int row = 0; row < nodeCnt; row++)
            {
                // zt
                for (int col = 0; col < edgeNodeCnt; col++)
                {
                    A[offset + row, col] = 0.0;
                    B[offset + row, col] = -Szt[row, col] - Uzn[row, col];
                }
                // zz
                for (int col = 0; col < nodeCnt; col++)
                {
                    A[offset + row, offset + col] = 0.0;
                    B[offset + row, offset + col] = -Tzz[row, col] - Uzz[row, col];
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
            success = (ret == 0);

            SortEVals(eVals, eVecs);
            AdjustPhaseAndVecDirEVecs(eVecs);

            // 固有値はβ^2
            double maxEVal = MaxEp < double.MaxValue ? k0 * k0 * MaxEp : double.MaxValue;
            int hitModeIndex = GetSameModeIndex(eVals, eVecs, maxEVal);
            if (hitModeIndex != -1)
            {
                System.Numerics.Complex[] newEVals = new System.Numerics.Complex[1];
                System.Numerics.Complex[][] newEVecs = new System.Numerics.Complex[1][];
                newEVals[0] = eVals[hitModeIndex];
                newEVecs[0] = eVecs[hitModeIndex];
                eVals = newEVals;
                eVecs = newEVecs;
            }
            else
            {
                success = false;
                // fail safe
                int n = A.RowLength;
                int m = 1;
                eVals = new System.Numerics.Complex[m];
                eVecs = new System.Numerics.Complex[m][];
                for (int iMode = 0; iMode < m; iMode++)
                {
                    eVecs[iMode] = new System.Numerics.Complex[n];
                }
            }

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

            return success;
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
            int sumNodeCnt = edgeNodeCnt + nodeCnt;
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
        /////////////////////////////////////////////////////////////////////

        private int GetSameModeIndex(
            System.Numerics.Complex[] eVals, System.Numerics.Complex[][] eVecs, double maxEVals)
        {
            int hitModeIndex = -1;
            if (PrevEVec == null)
            {
                hitModeIndex = -1;
                int modeCnt = eVals.Length;
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    System.Numerics.Complex eVal = eVals[iMode];
                    if (eVal.Real < maxEVals && eVal.Real >= 0.0)
                    {
                        hitModeIndex = iMode;
                        break;
                    }
                }
            }
            else
            {
                int modeCnt = eVals.Length;
                int matLen = eVecs[0].Length;
                hitModeIndex = -1;
                double maxAbsNorm = 0.0;
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    System.Numerics.Complex eVal = eVals[iMode];
                    if (eVal.Real < maxEVals && eVal.Real >= 0.0)
                    {
                        // OK
                    }
                    else
                    {
                        // NG
                        continue;
                    }
                    var eVec = eVecs[iMode];
                    System.Numerics.Complex[] workEVec1 = new System.Numerics.Complex[matLen];
                    PrevEVec.CopyTo(workEVec1, 0);
                    System.Numerics.Complex[] workEVec2 = new System.Numerics.Complex[matLen];
                    eVec.CopyTo(workEVec2, 0);
                    System.Numerics.Complex norm1 =
                        IvyFEM.Lapack.Functions.zdotu(IvyFEM.Lapack.Utils.Conjugate(workEVec1), workEVec1);
                    norm1 = System.Numerics.Complex.Sqrt(norm1);
                    System.Numerics.Complex norm2 =
                        IvyFEM.Lapack.Functions.zdotu(IvyFEM.Lapack.Utils.Conjugate(workEVec2), workEVec2);
                    norm2 = System.Numerics.Complex.Sqrt(norm2);
                    workEVec1 = IvyFEM.Lapack.Functions.zscal(workEVec1, 1.0 / norm1.Magnitude);
                    workEVec2 = IvyFEM.Lapack.Functions.zscal(workEVec2, 1.0 / norm2.Magnitude);
                    System.Numerics.Complex norm12 =
                        IvyFEM.Lapack.Functions.zdotu(IvyFEM.Lapack.Utils.Conjugate(workEVec1), workEVec2);
                    double absNorm12 = norm12.Magnitude;
                    if (absNorm12 > maxAbsNorm)
                    {
                        hitModeIndex = iMode;
                        maxAbsNorm = absNorm12;
                    }
                }
                System.Diagnostics.Debug.WriteLine("hitModeIndex = {0} absNorm = {1}", hitModeIndex, maxAbsNorm);
            }
            if (hitModeIndex != -1)
            {
                PrevEVec = eVecs[hitModeIndex];
            }
            return hitModeIndex;
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
