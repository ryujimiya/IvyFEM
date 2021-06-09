using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class EMWaveguide2DEigenFEMForPort : FEM
    {
        public uint PortId { get; protected set; } = 0;
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

        public EMWaveguide2DEigenFEMForPort(FEWorld world, uint portId)
        {
            World = world;
            PortId = portId;
        }

        private void CalcMatrixs(double k0)
        {
            uint tQuantityId = 0;
            uint zQuantityId = 1;
            int edgeNodeCnt = (int)World.GetPortEdgeNodeCount(tQuantityId, PortId);
            int nodeCnt = (int)World.GetPortNodeCount(zQuantityId, PortId);

            Rtt = new IvyFEM.Lapack.DoubleMatrix(edgeNodeCnt, edgeNodeCnt);
            Stz = new IvyFEM.Lapack.DoubleMatrix(edgeNodeCnt, nodeCnt);
            Szt = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, edgeNodeCnt);
            Ttt = new IvyFEM.Lapack.DoubleMatrix(edgeNodeCnt, edgeNodeCnt);
            Tzz = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);

            IList<uint> tfeIds = World.GetPortTriangleFEIds(tQuantityId, PortId);
            //--------------------
            double[] origin;
            double[] normal;
            double[] xdir;
            {
                var portCondition = World.GetPortConditions(tQuantityId)[(int)PortId];
                IList<int> eIds = portCondition.IntAdditionalParameters;
                System.Diagnostics.Debug.Assert(eIds.Count == 2);
                Mesher3D mesh = World.Mesh as Mesher3D;
                Cad3D cad = mesh.Cad;
                OpenTK.Vector3d originVec = new OpenTK.Vector3d();
                OpenTK.Vector3d[] dirVecs = new OpenTK.Vector3d[2];
                for (int i = 0; i < eIds.Count; i++)
                {
                    uint eId = (uint)eIds[i];
                    uint sVId;
                    uint eVId;
                    cad.GetEdgeVertexId(eId, out sVId, out eVId);
                    OpenTK.Vector3d v1 = cad.GetVertexCoord(sVId);
                    OpenTK.Vector3d v2 = cad.GetVertexCoord(eVId);
                    if (i == 0)
                    {
                        originVec = v1;
                    }
                    var dirVec = v2 - v1;
                    dirVec.Normalize();
                    dirVecs[i] = dirVec;
                
                }
                var normalVec = OpenTK.Vector3d.Cross(dirVecs[0], dirVecs[1]);

                origin = new double[] { originVec.X, originVec.Y, originVec.Z };
                normal = new double[] { normalVec.X, normalVec.Y, normalVec.Z };
                xdir = new double[] { dirVecs[0].X, dirVecs[0].Y, dirVecs[0].Z };
            }
            //--------------------
            foreach (uint feId in tfeIds)
            {
                // t成分
                TriangleFE tTriFE = World.GetTriangleFE(tQuantityId, feId);
                //--------------------
                tTriFE.FixedOrigin3D = origin;
                tTriFE.FixedNormal3D = normal;
                tTriFE.FixedXDir3D = xdir;
                //--------------------
                uint elemEdgeNodeCnt = tTriFE.EdgeCount;
                int[] edgeIds = new int[elemEdgeNodeCnt];
                bool[] isReverses = new bool[elemEdgeNodeCnt];
                int[] edgeNodes = new int[elemEdgeNodeCnt];
                for (int iENode = 0; iENode < elemEdgeNodeCnt; iENode++)
                {
                    int[] coIds = tTriFE.EdgeCoordIdss[iENode];
                    bool isReverse;
                    int edgeId = World.GetEdgeIdFromCoords(tQuantityId, coIds[0], coIds[1], out isReverse);
                    int edgeNodeId = World.PortEdge2EdgeNode(tQuantityId, PortId, edgeId);

                    edgeIds[iENode] = edgeId;
                    isReverses[iENode] = isReverse;
                    edgeNodes[iENode] = edgeNodeId;
                }
                // z成分
                TriangleFE zTriFE = World.GetTriangleFE(zQuantityId, feId);
                //--------------------
                zTriFE.FixedOrigin3D = origin;
                zTriFE.FixedNormal3D = normal;
                zTriFE.FixedXDir3D = xdir;
                //--------------------
                uint elemNodeCnt = zTriFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = zTriFE.NodeCoordIds[iNode];
                    int nodeId = World.PortCoord2Node(zQuantityId, PortId, coId);
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

            int edgeNodeCnt = (int)World.GetPortEdgeNodeCount(tQuantityId, PortId);
            int nodeCnt = (int)World.GetPortNodeCount(zQuantityId, PortId);
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
            System.Numerics.Complex[] betas;
            System.Numerics.Complex[][] etEVecs;
            System.Numerics.Complex[][] ezEVecs;
            GetBetasEtzVecs(omega, eVals, eVecs, out betas, out etEVecs, out ezEVecs);
            Betas = betas;
            EtEVecs = etEVecs;
            EzEVecs = ezEVecs;

            /*
            //---------------------------------------
            //!!!!!!!!!!!!!TEST
            // 伝搬モードのみに絞る
            var propBetas = new List<System.Numerics.Complex>();
            var propEtEVecs = new List<System.Numerics.Complex[]>();
            var propEzEVecs = new List<System.Numerics.Complex[]>();
            for (int iMode = 0; iMode < Betas.Length; iMode++)
            {
                if (Math.Abs(Betas[iMode].Real / k0) < 1.0e-6)
                {
                    continue;
                }
                propBetas.Add(Betas[iMode]);
                propEtEVecs.Add(EtEVecs[iMode]);
                propEzEVecs.Add(EzEVecs[iMode]);
            }
            if (propBetas.Count == 0)
            {
                Betas = new System.Numerics.Complex[] { 0.0 };
                EtEVecs = new System.Numerics.Complex[][] { new System.Numerics.Complex[edgeNodeCnt] };
                EzEVecs = new System.Numerics.Complex[][] { new System.Numerics.Complex[nodeCnt] };
            }
            else
            {
                Betas = propBetas.ToArray();
                EtEVecs = propEtEVecs.ToArray();
                EzEVecs = propEzEVecs.ToArray();
            }
            //---------------------------------------
            */
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
            int edgeNodeCnt = (int)World.GetPortEdgeNodeCount(tQuantityId, PortId);
            int nodeCnt = (int)World.GetPortNodeCount(zQuantityId, PortId);
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
                        maxEdgeId = World.PortEdgeNode2Edge(tQuantityId, PortId, iENode);
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

        private void GetBetasEtzVecs(double omega,
            System.Numerics.Complex[] eVals, System.Numerics.Complex[][] eVecs,
            out System.Numerics.Complex[] Betas,
            out System.Numerics.Complex[][] EtEVecs, out System.Numerics.Complex[][] EzEVecs)
        {
            System.Diagnostics.Debug.Assert(!IsMagneticField);

            var RttZ = new IvyFEM.Lapack.ComplexMatrix(Rtt);
            var SztZ = new IvyFEM.Lapack.ComplexMatrix(Szt);

            uint tQuantityId = 0;
            uint zQuantityId = 1;
            int edgeNodeCnt = (int)World.GetPortEdgeNodeCount(tQuantityId, PortId);
            int nodeCnt = (int)World.GetPortNodeCount(zQuantityId, PortId);
            int sumNodeCnt = edgeNodeCnt + nodeCnt;
            int offset = edgeNodeCnt;

            int modeCnt = eVals.Length;
            Betas = new System.Numerics.Complex[modeCnt];
            EtEVecs = new System.Numerics.Complex[modeCnt][];
            EzEVecs = new System.Numerics.Complex[modeCnt][];
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                System.Numerics.Complex beta = System.Numerics.Complex.Sqrt(eVals[iMode]);
                if (beta.Magnitude < 1.0e-20)
                {
                    beta = -1.0e-20 * System.Numerics.Complex.ImaginaryOne;
                }
                if (beta.Imaginary > 0)
                {
                    beta = System.Numerics.Complex.Conjugate(beta);
                }
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
                // ez = j ezj Nj
                for (int nodeId = 0; nodeId < nodeCnt; nodeId++)
                {
                    ez[nodeId] *= System.Numerics.Complex.ImaginaryOne;
                }

                var work1 = RttZ * et;
                var work2 = IvyFEM.Lapack.Functions.zdotc(et, work1);
                var work3 = SztZ * et;
                var work4 = IvyFEM.Lapack.Functions.zdotc(ez, work3);

                System.Numerics.Complex dominator =
                    System.Numerics.Complex.ImaginaryOne * System.Numerics.Complex.Conjugate(beta) * work2 -
                    work4;
                System.Numerics.Complex d = System.Numerics.Complex.Sqrt(
                    (System.Numerics.Complex.Conjugate(beta) / beta.Magnitude) *
                    (System.Numerics.Complex.ImaginaryOne* omega * Constants.Mu0) /
                    dominator);

                et = IvyFEM.Lapack.Functions.zscal(et, d);
                ez = IvyFEM.Lapack.Functions.zscal(ez, d);

                EtEVecs[iMode] = et;
                EzEVecs[iMode] = ez;
            }
        }

        public IvyFEM.Lapack.ComplexMatrix CalcBoundaryMatrix(
            double omega,
            System.Numerics.Complex[] betas,
            System.Numerics.Complex[][] etEVecs,
            System.Numerics.Complex[][] ezEVecs)
        {
            System.Diagnostics.Debug.Assert(!IsMagneticField);

            var StzZ = new IvyFEM.Lapack.ComplexMatrix(Stz);
            var SztZT = new IvyFEM.Lapack.ComplexMatrix(Szt);
            SztZT.Transpose();
            var RttZ = new IvyFEM.Lapack.ComplexMatrix(Rtt);

            int edgeNodeCnt = etEVecs[0].Length;
            int nodeCnt = ezEVecs[0].Length;
            IvyFEM.Lapack.ComplexMatrix X = new IvyFEM.Lapack.ComplexMatrix(edgeNodeCnt, edgeNodeCnt);
            int modeCnt = betas.Length;
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var beta = betas[iMode];
                var etEVec = etEVecs[iMode];
                var ezEVec = ezEVecs[iMode];

                System.Numerics.Complex[] conjEtEVec = IvyFEM.Lapack.Utils.Conjugate(etEVec);
                System.Numerics.Complex[] conjEzEVec = IvyFEM.Lapack.Utils.Conjugate(ezEVec);

                System.Numerics.Complex[] vec1 = new System.Numerics.Complex[edgeNodeCnt];
                System.Numerics.Complex[] vec2 = new System.Numerics.Complex[edgeNodeCnt];

                ////////////////////////////////////
                // vec2
                var work1 = SztZT * conjEzEVec;
                var work2 = RttZ * conjEtEVec;
                for (int i = 0; i < edgeNodeCnt; i++)
                {
                    vec2[i] = -work1[i] +
                        System.Numerics.Complex.ImaginaryOne * System.Numerics.Complex.Conjugate(beta) * work2[i];
                }

                ///////////////////////////////////
                // vec1
                var work3 = StzZ * ezEVec;
                var work4 = RttZ * etEVec;
                for (int i = 0; i < edgeNodeCnt; i++)
                {
                    vec1[i] =
                        -(1.0 / (System.Numerics.Complex.ImaginaryOne * omega * Constants.Mu0)) *
                        (work3[i] + System.Numerics.Complex.ImaginaryOne * beta * work4[i]);
                }

                // X
                for (int row = 0; row < edgeNodeCnt; row++)
                {
                    for (int col = 0; col < edgeNodeCnt; col++)
                    {
                        System.Numerics.Complex value =
                            (beta.Magnitude / System.Numerics.Complex.Conjugate(beta)) *
                            vec1[row] * vec2[col];
                        X[row, col] += value;
                    }
                }
            }
            return X;
        }

        public System.Numerics.Complex[] CalcIncidentVec(
            double omega,
            System.Numerics.Complex beta0, System.Numerics.Complex[] etEVec0, System.Numerics.Complex[] ezEVec0)
        {
            System.Diagnostics.Debug.Assert(!IsMagneticField);

            var StzZ = new IvyFEM.Lapack.ComplexMatrix(Stz);
            var RttZ = new IvyFEM.Lapack.ComplexMatrix(Rtt);

            int edgeNodeCnt = etEVec0.Length;
            int nodeCnt = ezEVec0.Length;
            System.Numerics.Complex[] I = new System.Numerics.Complex[edgeNodeCnt];
            {
                var beta = beta0;
                var etEVec = etEVec0;
                var ezEVec = ezEVec0;

                System.Numerics.Complex[] vec1 = new System.Numerics.Complex[edgeNodeCnt];

                var work3 = StzZ * ezEVec;
                var work4 = RttZ * etEVec;
                for (int i = 0; i < edgeNodeCnt; i++)
                {
                    vec1[i] = 
                        -(1.0 / (System.Numerics.Complex.ImaginaryOne * omega * Constants.Mu0)) *
                        (work3[i] + System.Numerics.Complex.ImaginaryOne * beta * work4[i]);
                }

                // I
                for (int row = 0; row < edgeNodeCnt; row++)
                {
                    I[row] = -1.0 * System.Numerics.Complex.ImaginaryOne * omega * Constants.Mu0 * 2.0 * vec1[row];
                }
            }
            return I;
        }

        public System.Numerics.Complex[] CalcSMatrix(
            double omega, int incidentModeId,
            System.Numerics.Complex[] betas,
            System.Numerics.Complex[][] etEVecs,
            System.Numerics.Complex[][] ezEVecs,
            System.Numerics.Complex[] Et)
        {
            System.Diagnostics.Debug.Assert(!IsMagneticField);

            var SztZ = new IvyFEM.Lapack.ComplexMatrix(Szt);
            var RttZ = new IvyFEM.Lapack.ComplexMatrix(Rtt);

            uint tQuantityId = 0;
            uint zQuantityId = 1;
            int edgeNodeCnt = (int)World.GetPortEdgeNodeCount(tQuantityId, PortId);
            int nodeCnt = (int)World.GetPortNodeCount(zQuantityId, PortId);
            int sumNodeCnt = edgeNodeCnt + nodeCnt;
            int offset = edgeNodeCnt;

            System.Diagnostics.Debug.Assert(Et.Length == edgeNodeCnt);
            int modeCnt = betas.Length;
            System.Numerics.Complex[] S = new System.Numerics.Complex[modeCnt];
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var beta = betas[iMode];
                var et = etEVecs[iMode];
                var ez = ezEVecs[iMode];

                var work1 = SztZ * Et;
                var work2 = IvyFEM.Lapack.Functions.zdotc(ez, work1);
                var work3 = RttZ * Et;
                var work4 = IvyFEM.Lapack.Functions.zdotc(et, work3);

                System.Numerics.Complex b =
                    (beta.Magnitude / System.Numerics.Complex.Conjugate(beta)) *
                    (1.0 / (System.Numerics.Complex.ImaginaryOne * omega * Constants.Mu0)) *
                    (-work2 +
                    System.Numerics.Complex.ImaginaryOne * System.Numerics.Complex.Conjugate(beta) * work4);
                if (incidentModeId != -1 && incidentModeId == iMode)
                {
                    b = (-1.0) + b;
                }
                S[iMode] = b;
            }

            return S;
        }

        public System.Numerics.Complex[] CalcModeCoordExyzByCoord(
            double[] coord,
            System.Numerics.Complex beta, System.Numerics.Complex[] etEVec, System.Numerics.Complex[] ezEVec)
        {
            uint tQuantityId = 0;
            uint zQuantityId = 1;

            IList<uint> portFEIds = World.GetPortTriangleFEIds(tQuantityId, PortId);
            IList<uint> feIds = World.GetTriangleFEsWithPointInside(tQuantityId, coord);
            int dof = 3;
            System.Numerics.Complex[] exyz = new System.Numerics.Complex[dof];
            int cnt = 0;
            System.Diagnostics.Debug.Assert(feIds.Count != 0);

            foreach (uint feId in feIds)
            {
                if (!portFEIds.Contains(feId))
                {
                    continue;
                }
                // t成分
                TriangleFE tTriFE = World.GetTriangleFE(tQuantityId, feId);
                System.Diagnostics.Debug.Assert(tTriFE.FixedOrigin3D != null);
                uint elemEdgeNodeCnt = tTriFE.EdgeCount;
                int[] edgeIds = new int[elemEdgeNodeCnt];
                bool[] isReverses = new bool[elemEdgeNodeCnt];
                int[] edgeNodes = new int[elemEdgeNodeCnt];
                for (int iENode = 0; iENode < elemEdgeNodeCnt; iENode++)
                {
                    int[] coIds = tTriFE.EdgeCoordIdss[iENode];
                    bool isReverse;
                    int edgeId = World.GetEdgeIdFromCoords(tQuantityId, coIds[0], coIds[1], out isReverse);
                    int edgeNodeId = World.PortEdge2EdgeNode(tQuantityId, PortId, edgeId);

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
                    int nodeId = World.PortCoord2Node(zQuantityId, PortId, coId);
                    nodes[iNode] = nodeId;
                }

                double[] L = zTriFE.Coord2L(coord);
                double[][] tN = tTriFE.CalcEdgeN(L);
                double[] N = zTriFE.CalcN(L);

                // t成分の値
                for (int iENode = 0; iENode < elemEdgeNodeCnt; iENode++)
                {
                    int edgeNodeId = edgeNodes[iENode];
                    if (edgeNodeId == -1)
                    {
                        continue;
                    }
                    double sgn = isReverses[iENode] ? -1.0 : 1.0;
                    for (int iDof = 0; iDof < 2; iDof++)
                    {
                        exyz[iDof] += sgn * etEVec[edgeNodeId] * tN[iENode][iDof];
                    }
                }
                // z成分の値
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int nodeId = nodes[iNode];
                    if (nodeId == -1)
                    {
                        continue;
                    }
                    int iDof = 2;
                    exyz[iDof] += ezEVec[nodeId] * N[iNode];
                }
                cnt++;
            }
            for (int iDof = 0; iDof < dof; iDof++)
            {
                exyz[iDof] /= cnt;
            }
            return exyz;
        }
    }
}
