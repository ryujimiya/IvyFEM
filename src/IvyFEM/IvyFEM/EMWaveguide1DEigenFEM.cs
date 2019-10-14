using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class EMWaveguide1DEigenFEM : FEM
    {
        public uint QuantityId { get; private set; } = 0;
        public uint PortId { get; private set; } = 0;
        public EMWaveguideType WaveguideType { get; set; } = EMWaveguideType.HPlane2D;
        public double WaveguideWidthForEPlane { get; set; } = 0;
        /// <summary>
        /// TEモードで実装した式をTMモードに流用するため
        ///   TEモードの場合は μ0
        ///   TMモードの場合は ε0
        /// </summary>
        public double ReplacedMu0 { get; set; } = Constants.Mu0;

        public IvyFEM.Lapack.DoubleMatrix Txx { get; private set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Ryy { get; private set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Uzz { get; private set; } = null;

        // Solve
        // Input
        public double Frequency { get; set; }
        // Output
        public System.Numerics.Complex[] Betas { get; private set; }
        public System.Numerics.Complex[][] EzEVecs { get; private set; }

        public EMWaveguide1DEigenFEM(FEWorld world, uint quantityId, uint portId)
        {
            World = world;
            QuantityId = quantityId;
            PortId = portId;
        }

        private uint GetMaId0()
        {
            IList<uint> feIds = World.GetPortLineFEIds(QuantityId, PortId);
            uint feId0 = feIds[0];
            LineFE lineFE = World.GetLineFE(QuantityId, feId0);
            uint maId0 = lineFE.MaterialId;
            return maId0;
        }

        private double GetEpRatioForEPlane(double k0)
        {
            uint maId0 = GetMaId0();
            DielectricMaterial ma0 = World.GetMaterial(maId0) as DielectricMaterial;

            // εyの式
            double barEp = (
                ma0.Epyy -
                Math.PI * Math.PI / (k0 * k0 * WaveguideWidthForEPlane * WaveguideWidthForEPlane * ma0.Muxx)
                ) /
                ma0.Epyy;
            if (Math.Abs(barEp) < Constants.PrecisionLowerLimit)
            {
                barEp = barEp >= 0 ? Constants.PrecisionLowerLimit : -Constants.PrecisionLowerLimit;
            }
            return barEp;
        }

        private void CalcMatrixs(double k0)
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

        public override void Solve()
        {
            Betas = null;
            EzEVecs = null;

            // 角周波数
            double omega = 2.0 * Math.PI * Frequency;
            // 波数
            double k0 = omega / Constants.C0;

            CalcMatrixs(k0);

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

            // 質量行列の正定値マトリクスチェック
            if (WaveguideType == EMWaveguideType.EPlane2D)
            {
                // 比誘電率を置き換えているため、正定値行列にならないことがある（おそらくTE10モードが減衰モードのとき）
                bool isPositiveDefinite = true;
                for (int i = 0; i < B.RowLength; i++)
                {
                    // 対角成分の符号チェック
                    if (B[i, i] < 0.0)
                    {
                        isPositiveDefinite = false;
                        break;
                    }
                }
                if (!isPositiveDefinite)
                {
                    for (int i = 0; i < A.RowLength; i++)
                    {
                        for (int j =0; j < A.ColumnLength; j++)
                        {
                            A[i, j] = -1.0 * A[i, j];
                        }
                    }
                    for (int i = 0; i < B.RowLength; i++)
                    {
                        for (int j = 0; j < B.ColumnLength; j++)
                        {
                            B[i, j] = -1.0 * B[i, j];
                        }
                    }
                }
            }

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
            // 波数
            double k0 = omega / Constants.C0;

            uint maId0 = GetMaId0();
            DielectricMaterial ma0 = World.GetMaterial(maId0) as DielectricMaterial;
            double epRatio = 0; // for EPlane
            if (WaveguideType == EMWaveguideType.HPlane2D)
            {

            }
            else if (WaveguideType == EMWaveguideType.EPlane2D)
            {
                epRatio = GetEpRatioForEPlane(k0);
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }

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
                System.Numerics.Complex d = 0;
                if (WaveguideType == EMWaveguideType.HPlane2D)
                {
                    d = System.Numerics.Complex.Sqrt(
                        (System.Numerics.Complex)(omega * ReplacedMu0) /
                        (((System.Numerics.Complex)beta.Magnitude) * work2));
                }
                else if (WaveguideType == EMWaveguideType.EPlane2D)
                {
                    d = System.Numerics.Complex.Sqrt(
                        omega * Constants.Ep0 * 
                        System.Numerics.Complex.Conjugate(epRatio * ma0.Epyy * (1.0 / ma0.Muyy)) /
                        (((System.Numerics.Complex)beta.Magnitude) * work2));
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }

                eVec = IvyFEM.Lapack.Functions.zscal(eVec, d);
                eVecs[iMode] = eVec;
            }
        }

        public IvyFEM.Lapack.ComplexMatrix CalcBoundaryMatrix(
            double omega, System.Numerics.Complex[] betas, System.Numerics.Complex[][] ezEVecs)
        {
            // 波数
            double k0 = omega / Constants.C0;

            uint maId0 = GetMaId0();
            DielectricMaterial ma0 = World.GetMaterial(maId0) as DielectricMaterial;
            double barEp = 0; // for EPlane
            if (WaveguideType == EMWaveguideType.HPlane2D)
            {

            }
            else if (WaveguideType == EMWaveguideType.EPlane2D)
            {
                barEp = GetEpRatioForEPlane(k0);
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }

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
                        System.Numerics.Complex value = 0;
                        if (WaveguideType == EMWaveguideType.HPlane2D)
                        {
                            value = (System.Numerics.Complex.ImaginaryOne / (omega * ReplacedMu0)) *
                                beta * beta.Magnitude *
                                vec1[row] * vec2[col];

                        }
                        else if (WaveguideType == EMWaveguideType.EPlane2D)
                        {
                            value = (System.Numerics.Complex.ImaginaryOne / 
                                (omega * Constants.Ep0 * 
                                System.Numerics.Complex.Conjugate(barEp * ma0.Epyy * (1.0 / ma0.Muyy))))
                                * beta * beta.Magnitude *
                                vec1[row] * vec2[col];
                        }
                        else
                        {
                            System.Diagnostics.Debug.Assert(false);
                        }
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
            // 波数
            double k0 = omega / Constants.C0;

            uint maId0 = GetMaId0();
            DielectricMaterial ma0 = World.GetMaterial(maId0) as DielectricMaterial;
            double epRatio = 0; // for EPlane
            if (WaveguideType == EMWaveguideType.HPlane2D)
            {

            }
            else if (WaveguideType == EMWaveguideType.EPlane2D)
            {
                epRatio = GetEpRatioForEPlane(k0);
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }

            int modeCnt = betas.Length;
            System.Numerics.Complex[] S = new System.Numerics.Complex[modeCnt];

            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var beta = betas[iMode];
                var ezEVec = ezEVecs[iMode];
                var RyyZ = (IvyFEM.Lapack.ComplexMatrix)Ryy;
                var vec1 = RyyZ * IvyFEM.Lapack.Utils.Conjugate(ezEVec);
                System.Numerics.Complex work1 = IvyFEM.Lapack.Functions.zdotu(vec1, Ez);
                System.Numerics.Complex b = 0;
                if (WaveguideType == EMWaveguideType.HPlane2D)
                {

                    b = (beta.Magnitude / (omega * ReplacedMu0)) * work1;
                }
                else if (WaveguideType == EMWaveguideType.EPlane2D)
                {
                    b = (beta.Magnitude / (omega * Constants.Ep0 *
                        System.Numerics.Complex.Conjugate(ma0.Epyy * epRatio * (1.0 / ma0.Muyy))
                        )) * work1;

                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
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
