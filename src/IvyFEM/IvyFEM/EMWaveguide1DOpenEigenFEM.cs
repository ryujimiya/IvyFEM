﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class EMWaveguide1DOpenEigenFEM : FEM
    {
        public uint QuantityId { get; private set; } = 0;
        public uint PortId { get; private set; } = 0;
        /// <summary>
        /// TEモードで実装した式をTMモードに流用するため
        ///   TEモードの場合は μ0
        ///   TMモードの場合は ε0
        /// </summary>
        public double ReplacedMu0 { get; set; } = Constants.Mu0;
        /// <summary>
        /// クラッドの比誘電率
        /// </summary>
        public double CladdingEp { get; set; } = 0.0;
        /// <summary>
        /// 減衰定数
        /// </summary>
        public double DecayParameter { get; set; } = 0.0;

        public IvyFEM.Lapack.DoubleMatrix Txx { get; private set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Ryy { get; private set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Uzz { get; private set; } = null;

        // Solve
        // Input
        public double Frequency { get; set; }
        // Output
        public System.Numerics.Complex[] Betas { get; private set; }
        public System.Numerics.Complex[][] EzEVecs { get; private set; }

        public EMWaveguide1DOpenEigenFEM(FEWorld world, uint quantityId, uint portId)
        {
            World = world;
            QuantityId = quantityId;
            PortId = portId;
        }

        private void CalcMatrixs()
        {
            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            IList<uint> feIds = World.GetPortLineFEIds(QuantityId, PortId);

            Txx = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            Ryy = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            Uzz = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);

            foreach (uint feId in feIds)
            {
                LineFE lineFE = World.GetLineFE(QuantityId, feId);
                uint elemNodeCnt = lineFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = lineFE.NodeCoordIds[iNode];
                    int nodeId = World.PortCoord2Node(QuantityId, PortId, coId);
                    nodes[iNode] = nodeId;
                }

                Material ma0 = World.GetMaterial(lineFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is DielectricMaterial);
                var ma = ma0 as DielectricMaterial;

                double maPxx = 1.0 / ma.Muxx;
                double maPyy = 1.0 / ma.Muyy;
                double maQzz = ma.Epzz;

                double[,] sNN = lineFE.CalcSNN();
                double[,] sNyNy = lineFE.CalcSNxNx();
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
                        double txxVal = maPxx * sNyNy[row, col];
                        double ryyVal = maPyy * sNN[row, col];
                        double uzzVal = maQzz * sNN[row, col];

                        Txx[rowNodeId, colNodeId] += txxVal;
                        Ryy[rowNodeId, colNodeId] += ryyVal;
                        Uzz[rowNodeId, colNodeId] += uzzVal;
                    }
                }
            }
        }

        private void CalcOpenBoundary(double k0, IvyFEM.Lapack.DoubleMatrix A, IvyFEM.Lapack.DoubleMatrix B)
        {
            double alpha = DecayParameter;
            if (Math.Abs(alpha) < 1.0e-60)
            {
                return;
            }
            IList<int> cornerCoIds = new List<int>(){ -1, -1 };
            {
                IList<PortCondition> portcontitions = World.GetPortConditions(QuantityId);
                PortCondition portcondition = portcontitions[(int)PortId];
                IList<uint> eIds = portcondition.EIds;

                uint eId1 = eIds[0];
                Edge2D e1 = World.Mesh.Cad2D.GetEdge(eId1);
                uint vId1 = e1.GetVertexId(true); // 始点
                IList<int> coIds1 = World.GetCoordIdsFromCadId(QuantityId, vId1, CadElementType.Vertex);
                System.Diagnostics.Debug.Assert(coIds1.Count == 1);
                int coId1 = coIds1[0];
                cornerCoIds[0] = coId1;

                uint eId2 = eIds[eIds.Count - 1];
                Edge2D e2 = World.Mesh.Cad2D.GetEdge(eId2);
                uint vId2 = e2.GetVertexId(false); // 終点
                IList<int> coIds2 = World.GetCoordIdsFromCadId(QuantityId, vId2, CadElementType.Vertex);
                System.Diagnostics.Debug.Assert(coIds2.Count == 1);
                int coId2 = coIds2[0];
                cornerCoIds[1] = coId2;
            }
            IList<uint> feIds = World.GetPortLineFEIds(QuantityId, PortId);
            foreach (uint feId in feIds)
            {
                LineFE lineFE = World.GetLineFE(QuantityId, feId);
                uint elemNodeCnt = lineFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = lineFE.NodeCoordIds[iNode];
                    int nodeId = World.PortCoord2Node(QuantityId, PortId, coId);
                    nodes[iNode] = nodeId;
                }

                Material ma0 = World.GetMaterial(lineFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is DielectricMaterial);
                var ma = ma0 as DielectricMaterial;

                double maPxx = 1.0 / ma.Muxx;
                double maPyy = 1.0 / ma.Muyy;
                double maQzz = ma.Epzz;

                for (int row = 0; row < elemNodeCnt; row++)
                {
                    int rowCoId = lineFE.NodeCoordIds[row];
                    int rowNodeId = nodes[row];
                    if (rowNodeId == -1)
                    {
                        continue;
                    }
                    for (int col = 0; col < elemNodeCnt; col++)
                    {
                        int colCoId = lineFE.NodeCoordIds[col];
                        int colNodeId = nodes[col];
                        if (colNodeId == -1)
                        {
                            continue;
                        }

                        if ((rowCoId == colCoId) &&
                            cornerCoIds.Contains(rowCoId))
                        {
                            // 無限要素
                            double txxVal = maPxx * (alpha / 2.0);
                            double ryyVal = maPyy * (1.0 / (2.0 * alpha));
                            double uzzVal = maQzz * (1.0 / (2.0 * alpha));
                            A[rowNodeId, colNodeId] += (k0 * k0) * uzzVal - txxVal;
                            B[rowNodeId, colNodeId] += ryyVal;
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

            int iterCnt = 1000;
            int iter = 0;
            double prevBeta = k0 * Math.Sqrt(CladdingEp);
            double convRatio = 0;
            for (iter = 0; iter < iterCnt; iter++)
            {
                double alpha2 = prevBeta * prevBeta - k0 * k0 * CladdingEp;
                if (alpha2 < 0.0)
                {
                    alpha2 = 0.0;
                }
                DecayParameter = Math.Sqrt(alpha2);

                SolveIter();

                double beta = Betas[0].Real;
                convRatio = Math.Abs((beta - prevBeta) / beta);
                if (convRatio < 1.0e-6)
                {
                    break;
                }
                prevBeta = beta;
            }
            System.Diagnostics.Debug.Assert(iter < iterCnt);
        }

        private void SolveIter()
        {
            Betas = null;
            EzEVecs = null;

            // 角周波数
            double omega = 2.0 * Math.PI * Frequency;
            // 波数
            double k0 = omega / Constants.C0;

            CalcMatrixs();

            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            var A = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            var B = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            for (int col = 0; col < nodeCnt; col++)
            {
                for (int row = 0; row < nodeCnt; row++)
                {
                    A[row, col] = (k0 * k0) * Uzz[row, col] - Txx[row, col];
                    B[row, col] = Ryy[row, col];
                }
            }

            CalcOpenBoundary(k0, A, B);

            System.Numerics.Complex[] eVals;
            System.Numerics.Complex[][] eVecs;
            int ret = IvyFEM.Lapack.Functions.dggev(A.Buffer, A.RowLength, A.ColumnLength,
                B.Buffer, B.RowLength, B.ColumnLength,
                out eVals, out eVecs);
            System.Diagnostics.Debug.Assert(ret == 0);

            SortEVals(eVals, eVecs);
            AdjustPhaseEVecs(eVecs);
            GetBetasEzVecs(omega, eVals, eVecs);
            Betas = eVals;
            EzEVecs = eVecs;
        }

        private void SortEVals(System.Numerics.Complex[] eVals, System.Numerics.Complex[][] eVecs)
        {
            int modeCnt = eVals.Length;
            var eValEVecs = new List<KeyValuePair<System.Numerics.Complex, System.Numerics.Complex[]>>();
            for (int i = 0; i < modeCnt;  i++)
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

        private void AdjustPhaseEVecs(System.Numerics.Complex[][] eVecs)
        {
            int modeCnt = eVecs.Length;
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var eVec = eVecs[iMode];
                int nodeCnt = eVec.Length;
                System.Numerics.Complex maxValue = new System.Numerics.Complex(0, 0);
                double maxAbs = 0;
                for (int iNode = 0; iNode < nodeCnt; iNode++)
                {
                    System.Numerics.Complex value = eVec[iNode];
                    double abs = value.Magnitude;
                    if (abs > maxAbs)
                    {
                        maxAbs = abs;
                        maxValue = value;
                    }
                }
                System.Numerics.Complex phase = maxValue / maxAbs;

                for (int iNode = 0; iNode < nodeCnt; iNode++)
                {
                    eVec[iNode] /= phase;
                }
            }
        }

        private void GetBetasEzVecs(
            double omega, System.Numerics.Complex[] eVals, System.Numerics.Complex[][] eVecs)
        {
            int modeCnt = eVals.Length;
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var eVal = eVals[iMode];
                var beta = System.Numerics.Complex.Sqrt(eVal);
                if (beta.Imaginary > 0)
                {
                    beta = System.Numerics.Complex.Conjugate(beta);
                }
                eVals[iMode] = beta;

                var eVec = eVecs[iMode];
                var RyyZ = (IvyFEM.Lapack.ComplexMatrix)Ryy;
                var work = RyyZ * eVec;
                var work2 = IvyFEM.Lapack.Functions.zdotc(eVec, work);
                System.Numerics.Complex d = 
                    System.Numerics.Complex.Sqrt(
                    (System.Numerics.Complex)(omega * ReplacedMu0) /
                    (((System.Numerics.Complex)beta.Magnitude) * work2));

                eVec = IvyFEM.Lapack.Functions.zscal(eVec, d);
                eVecs[iMode] = eVec;
            }
        }

        public IvyFEM.Lapack.ComplexMatrix CalcBoundaryMatrix(
            double omega, System.Numerics.Complex[] betas, System.Numerics.Complex[][] ezEVecs)
        {
            int nodeCnt = ezEVecs[0].Length;
            IvyFEM.Lapack.ComplexMatrix X = new Lapack.ComplexMatrix(nodeCnt, nodeCnt);

            int modeCnt = betas.Length;
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var beta = betas[iMode];
                var ezEVec = ezEVecs[iMode];
                var RyyZ = (IvyFEM.Lapack.ComplexMatrix)Ryy;
                var vec1 = RyyZ * ezEVec;
                var vec2 = RyyZ * IvyFEM.Lapack.Utils.Conjugate(ezEVec);

                for (int row = 0; row < nodeCnt; row++)
                {
                    for (int col = 0; col < nodeCnt; col++)
                    {
                        System.Numerics.Complex value =
                            (System.Numerics.Complex.ImaginaryOne / (omega * ReplacedMu0)) *
                            beta * beta.Magnitude *
                            vec1[row] * vec2[col];
                        X[row, col] += value;
                    }
                }
            }
            return X;
        }

        public System.Numerics.Complex[] CalcIncidentVec(
            System.Numerics.Complex beta0, System.Numerics.Complex[] ezEVec0)
        {
            System.Numerics.Complex[] I = null;

            var RyyZ = (IvyFEM.Lapack.ComplexMatrix)Ryy;
            var vec1 = RyyZ * ezEVec0;
            var a1 = System.Numerics.Complex.ImaginaryOne * 2.0 * beta0;
            vec1 = IvyFEM.Lapack.Functions.zscal(vec1, a1);
            I = vec1;
            return I;
        }

        public System.Numerics.Complex[] CalcSMatrix(double omega, int incidentModeId,
            System.Numerics.Complex[] betas, System.Numerics.Complex[][] ezEVecs,
            System.Numerics.Complex[] Ez)
        {
            int modeCnt = betas.Length;
            System.Numerics.Complex[] S = new System.Numerics.Complex[modeCnt];

            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var beta = betas[iMode];
                var ezEVec = ezEVecs[iMode];
                var RyyZ = (IvyFEM.Lapack.ComplexMatrix)Ryy;
                var vec1 = RyyZ * IvyFEM.Lapack.Utils.Conjugate(ezEVec);
                System.Numerics.Complex work1 = IvyFEM.Lapack.Functions.zdotu(vec1, Ez);
                System.Numerics.Complex b = (beta.Magnitude / (omega * ReplacedMu0)) * work1;
                if (incidentModeId != -1 && incidentModeId == iMode)
                {
                    b = (-1.0) + b;
                }
                S[iMode] = b;
            }

            return S;
        }
    }
}