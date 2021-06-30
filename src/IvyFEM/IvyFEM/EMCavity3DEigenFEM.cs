using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class EMCavity3DEigenFEM : FEM
    {
        // 磁界？
        public bool IsMagneticField { get; set; } = false;
        // 固有値ソルバー
        public bool IsSymmetricBandSolver { get; set; } = true;

        // Solve
        // output
        public System.Numerics.Complex[] Frequencys { get; protected set; }
        public System.Numerics.Complex[][] EVecs { get; protected set; }
        public System.Numerics.Complex[][] CoordExyzEVecs { get; protected set; }

        public EMCavity3DEigenFEM(FEWorld world)
        {
            World = world;
        }

        private void CalcMatrixs(IvyFEM.Lapack.DoubleMatrix K, IvyFEM.Lapack.DoubleMatrix M)
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

                    // K, M
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
                            K[rowEdgeNodeId, colEdgeNodeId] += rowEdgeSgn * colEdgeSgn * kVal;

                            double mVal = detJWeight * (
                                N[row][0] * maQxx * N[col][0] +
                                N[row][1] * maQyy * N[col][1] +
                                N[row][2] * maQzz * N[col][2]);
                            M[rowEdgeNodeId, colEdgeNodeId] += rowEdgeSgn * colEdgeSgn * mVal;
                        }
                    }
                }
            }
        }

        public override void Solve()
        {
            Frequencys = null;
            EVecs = null;

            uint quantityId = 0;

            int edgeNodeCnt = (int)World.GetEdgeNodeCount(quantityId);
            var K = new IvyFEM.Lapack.DoubleMatrix(edgeNodeCnt, edgeNodeCnt);
            var M = new IvyFEM.Lapack.DoubleMatrix(edgeNodeCnt, edgeNodeCnt);

            CalcMatrixs(K, M);

            System.Numerics.Complex[] eVals;
            System.Numerics.Complex[][] eVecs;
            if (IsSymmetricBandSolver)
            {
                //--------------------------------------------------
                var A = K;
                var B = M;
                double[] doubleEVals;
                double[][] doubleEVecs;
                SolveDoubleSymmetricDefiniteBandGeneralizedEigen(A, B, out doubleEVals, out doubleEVecs);

                // 互換性のため複素数にする
                eVals = new System.Numerics.Complex[doubleEVals.Length];
                eVecs = new System.Numerics.Complex[doubleEVals.Length][];
                for (int iMode = 0; iMode < doubleEVals.Length; iMode++)
                {
                    eVals[iMode] = new System.Numerics.Complex(doubleEVals[iMode], 0.0);
                    eVecs[iMode] = new System.Numerics.Complex[doubleEVecs[iMode].Length];
                    for (int nodeIdB = 0; nodeIdB < doubleEVecs[iMode].Length; nodeIdB++)
                    {
                        eVecs[iMode][nodeIdB] = new System.Numerics.Complex(doubleEVecs[iMode][nodeIdB], 0.0);
                    }
                }
                //--------------------------------------------------
            }
            else
            {
                //--------------------------------------------------
                var A = K;
                var B = M;
                SolveDoubleGeneralizedEigen(A, B, out eVals, out eVecs);
                //--------------------------------------------------
            }

            SortEVals(eVals, eVecs);
            AdjustPhaseAndVecDirEVecs(eVecs);

            //--------------------------------------------
            int modeCnt = eVals.Length;
            Frequencys = new System.Numerics.Complex[modeCnt];
            EVecs = new System.Numerics.Complex[modeCnt][];
            CoordExyzEVecs = new System.Numerics.Complex[modeCnt][];
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                System.Numerics.Complex k0 = System.Numerics.Complex.Sqrt(eVals[iMode]);
                System.Numerics.Complex omega = k0 * Constants.C0;
                System.Numerics.Complex freq = omega / (2.0 * Math.PI);
                Frequencys[iMode] = freq;
                EVecs[iMode] = eVecs[iMode];

                System.Numerics.Complex[] coordExyzEVec;
                CalcModeCoordExyz(eVecs[iMode], out coordExyzEVec);
                CoordExyzEVecs[iMode] = coordExyzEVec;
            }
            //--------------------------------------------
        }

        private void SolveDoubleGeneralizedEigen(
            IvyFEM.Lapack.DoubleMatrix A, IvyFEM.Lapack.DoubleMatrix B,
            out System.Numerics.Complex[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            eVals = null;
            eVecs = null;
            
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
        }

        private void SolveDoubleSymmetricDefiniteBandGeneralizedEigen(
            IvyFEM.Lapack.DoubleMatrix A, IvyFEM.Lapack.DoubleMatrix B,
            out double[] eVals, out double[][] eVecs)
        {
            eVals = null;
            eVecs = null;

            int matLen = A.RowLength;
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            System.Diagnostics.Debug.Assert(B.RowLength == B.ColumnLength);
            System.Diagnostics.Debug.Assert(A.RowLength == B.RowLength);

            // Bのバンド幅を縮小
            double[] dummyVec1 = new double[matLen];
            IvyFEM.Linear.DoubleSparseMatrix sparseB1 = new IvyFEM.Linear.DoubleSparseMatrix(B);
            IvyFEM.Linear.DoubleSparseMatrix sparseB2;
            double[] dummyVec2;
            int[] indexs;
            {
                bool successOrder = IvyFEM.Linear.Utils.OrderToDoubleBandMatrix(
                    out sparseB2, out dummyVec2, out indexs, sparseB1, dummyVec1);
                System.Diagnostics.Debug.Assert(successOrder);
                if (!successOrder)
                {
                    System.Diagnostics.Debug.Assert(false);
                    return;
                }
            }
            // AはBの並び替えにあわせる
            IvyFEM.Linear.DoubleSparseMatrix sparseA2 =
                new IvyFEM.Linear.DoubleSparseMatrix(matLen, matLen);
            for (int row = 0; row < matLen; row++)
            {
                for (int col = 0; col < matLen; col++)
                {
                    sparseA2[row, col] = A[indexs[row], indexs[col]];
                }
            }

            // 行列要素check
            {
                double th = 1.0e-12;
                for (int i = 0; i < matLen; i++)
                {
                    // [B]の正定値行列チェック
                    System.Diagnostics.Debug.Assert(sparseB2[i, i] > 0);
                    for (int j = i; j < matLen; j++)
                    {
                        // [A]は対称行列
                        System.Diagnostics.Debug.Assert(
                            Math.Abs(sparseA2[i, j] - sparseA2[j, i]) < th);
                        // [B]は対称行列
                        System.Diagnostics.Debug.Assert(
                            Math.Abs(sparseB2[i, j] - sparseB2[j, i]) < th);
                    }
                }
            }
            var symmetricBandA = new IvyFEM.Lapack.DoubleSymmetricBandMatrix(sparseA2);
            var symmetricBandB = new IvyFEM.Lapack.DoubleSymmetricBandMatrix(sparseB2);

            double[][] eVecs2;
            int ret = IvyFEM.Lapack.Functions.dsbgv(
                symmetricBandA.Buffer,
                symmetricBandA.RowLength, symmetricBandA.ColumnLength, symmetricBandA.SuperdiaLength,
                symmetricBandB.Buffer,
                symmetricBandB.RowLength, symmetricBandB.ColumnLength, symmetricBandB.SuperdiaLength,
                out eVals, out eVecs2);
            if (ret != 0)
            {
                // fail safe
                int n = A.RowLength;
                eVals = new double[n];
                eVecs = new double[n][];
                for (int iMode = 0; iMode < n; iMode++)
                {
                    eVecs[iMode] = new double[n];
                }
            }

            // eVecs2の各ベクトルを元の番号順に戻す
            int modeCnt = eVals.Length;
            eVecs = new double[modeCnt][];
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                double[] eVec2 = eVecs2[iMode];
                double[] eVec = new double[matLen];
                eVecs[iMode] = eVec;
                for (int i = 0; i < matLen; i++)
                {
                    eVec[indexs[i]] = eVec2[i];
                }
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
                // eVal(k0^2) の実部を比較
                double diff = a.Key.Real - b.Key.Real;
                // 昇順
                if (diff > 0)
                {
                    return 1;
                }
                else if (diff < 0)
                {
                    return -1;
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
            uint quantityId = 0;
            int modeCnt = eVecs.Length;
            int edgeNodeCnt = (int)World.GetEdgeNodeCount(quantityId);
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var eVec = eVecs[iMode];
                // 最大値を求める
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
                        maxEdgeId = World.EdgeNode2Edge(quantityId, iENode);
                    }
                }
                System.Numerics.Complex phase = maxValue / maxAbs;

                // 最大値を取る辺の方向
                double vecSgn = 1.0;
                {
                    int[] maxEdgeCoIds = World.GetEdgeCoordIds(quantityId, maxEdgeId);
                    int coId1 = maxEdgeCoIds[0];
                    int coId2 = maxEdgeCoIds[1];
                    double[] coord1 = World.GetCoord(quantityId, coId1);
                    double[] coord2 = World.GetCoord(quantityId, coId2);
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

        private void CalcModeCoordExyz(
            System.Numerics.Complex[] eEVec,
            out System.Numerics.Complex[] coordExyzEVec)
        {
            uint quantityId = 0;
            int edgeNodeCnt = (int)World.GetEdgeNodeCount(quantityId);
            int nodeCnt = (int)World.GetNodeCount(quantityId);
            int coCnt = (int)World.GetCoordCount(quantityId);

            int dof = 3; // x, y, z成分
            coordExyzEVec = new System.Numerics.Complex[coCnt * dof];
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
                    et[iENode] = sgn * eEVec[edgeNodeId];
                }

                // 節点のL
                double[][] nodeL = null;
                if (elemNodeCnt == 4)
                {
                    nodeL = new double[4][]
                    {
                        new double[4] { 1.0, 0.0, 0.0, 0.0 },
                        new double[4] { 0.0, 1.0, 0.0, 0.0 },
                        new double[4] { 0.0, 0.0, 1.0, 0.0 },
                        new double[4] { 0.0, 0.0, 0.0, 1.0 }
                    };
                }
                else if (elemNodeCnt == 10)
                {
                    nodeL = new double[10][]
                    {
                        new double[4] { 1.0, 0.0, 0.0, 0.0 },
                        new double[4] { 0.0, 1.0, 0.0, 0.0 },
                        new double[4] { 0.0, 0.0, 1.0, 0.0 },
                        new double[4] { 0.0, 0.0, 0.0, 1.0 },
                        new double[4] { 0.5, 0.5, 0.0, 0.0 },
                        new double[4] { 0.0, 0.5, 0.5, 0.0 },
                        new double[4] { 0.5, 0.0, 0.5, 0.0 },
                        new double[4] { 0.5, 0.0, 0.0, 0.5 },
                        new double[4] { 0.0, 0.5, 0.0, 0.5 },
                        new double[4] { 0.0, 0.0, 0.5, 0.5 }
                    };
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
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
                        coordExyzEVec[coId * dof + idim] += exyz[iNode][idim];
                    }
                    coordValueCnt[coId]++;
                }
            }

            for (int coId = 0; coId < coCnt; coId++)
            {
                for (int idim = 0; idim < dof; idim++)
                {
                    coordExyzEVec[coId * dof + idim] /= coordValueCnt[coId];
                }
            }
        }
    }
}
