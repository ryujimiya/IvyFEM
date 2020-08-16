using OpenTK.Graphics.ES11;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class ElasticSHWaveguide1DEigenFEM : FEM
    {
        public uint QuantityId { get; protected set; } = 0;
        public uint PortId { get; protected set; } = 0;

        public IvyFEM.Lapack.DoubleMatrix Rzz { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Tzz { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Mzz { get; protected set; } = null;
        //
        public IvyFEM.Lapack.DoubleMatrix SNN { get; protected set; } = null;

        // Solve
        // Input
        public double NormalX { get; set; } = 1.0;
        public double Frequency { get; set; }
        // Output
        public System.Numerics.Complex[] Betas { get; protected set; }
        public System.Numerics.Complex[][] UEVecs { get; protected set; } // uz
        public System.Numerics.Complex[][] HSigmaZXEVecs { get; protected set; } // (hat)sigmazx

        public ElasticSHWaveguide1DEigenFEM(FEWorld world, uint quantityId, uint portId)
        {
            World = world;
            QuantityId = quantityId;
            PortId = portId;
        }

        private void PostSolve()
        {
            int modeCnt = Betas.Length;
            // 応力σzxを求める
            HSigmaZXEVecs = new System.Numerics.Complex[modeCnt][];
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var beta = Betas[iMode];
                var uVec = UEVecs[iMode];
                HSigmaZXEVecs[iMode] = CalcStressValues(beta, uVec);
            }
        }

        // 応力σzxを求める
        private System.Numerics.Complex[] CalcStressValues(System.Numerics.Complex beta, System.Numerics.Complex[] uVec)
        {
            // 角周波数
            double omega = 2.0 * Math.PI * Frequency;

            uint uDof = 1;
            uint sDof = 1;
            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            System.Diagnostics.Debug.Assert(uVec.Length == nodeCnt * uDof);

            var sigmaZXVec = new System.Numerics.Complex[nodeCnt * sDof];

            IList<uint> feIds = World.GetPortLineFEIds(QuantityId, PortId);

            Dictionary<int, double> nodeLineLenSum = new Dictionary<int, double>();

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
                System.Diagnostics.Debug.Assert(ma0 is LinearElasticMaterial);
                var ma = ma0 as LinearElasticMaterial;
                double rho = ma.MassDensity;
                double lambda = ma.LameLambda;
                double mu = ma.LameMu;

                double lineLen = lineFE.GetLineLength();
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int iNodeId = nodes[iNode];
                    if (iNodeId == -1)
                    {
                        continue;
                    }
                    double[] L = lineFE.GetNodeL(iNode);
                    double[] N = lineFE.CalcN(L);
                    double[][] Nu = lineFE.CalcNu(L);
                    double[] Ny = Nu[0];

                    System.Numerics.Complex sigmaZX = 0;
                    for (int kNode = 0; kNode < elemNodeCnt; kNode++)
                    {
                        int kNodeId = nodes[kNode];
                        if (kNodeId == -1)
                        {
                            continue;
                        }
                        System.Numerics.Complex uz = uVec[kNodeId];
                        sigmaZX +=
                            -System.Numerics.Complex.ImaginaryOne * beta * mu * N[kNode] * uz;
                    }
                    sigmaZXVec[iNodeId * sDof] += sigmaZX * lineLen;

                    if (nodeLineLenSum.ContainsKey(iNodeId))
                    {
                        nodeLineLenSum[iNodeId] += lineLen;
                    }
                    else
                    {
                        nodeLineLenSum[iNodeId] = lineLen;
                    }
                }
            }

            System.Diagnostics.Debug.Assert(nodeLineLenSum.Count == nodeCnt);
            for (int i = 0; i < nodeCnt; i++)
            {
                double lineLenSum = nodeLineLenSum[i];
                sigmaZXVec[i] /= lineLenSum;
            }

            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            // 規格化
            for (int i = 0; i < nodeCnt; i++)
            {
                System.Numerics.Complex sigmaZX = sigmaZXVec[i];
                //!!!!!!!!!!!!!!!!!!
                //------------------
                sigmaZX *= omega;
                //------------------

                sigmaZXVec[i] = sigmaZX;
            }
            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            // hatσzx = jσzxに変換
            for (int i = 0; i < nodeCnt; i++)
            {
                // σzx
                sigmaZXVec[i] *= System.Numerics.Complex.ImaginaryOne;
            }

            return sigmaZXVec;
        }

        private void CalcMatrixs()
        {
            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            IList<uint> feIds = World.GetPortLineFEIds(QuantityId, PortId);

            Rzz = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            Tzz = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            Mzz = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            //
            SNN = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);

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
                System.Diagnostics.Debug.Assert(ma0 is LinearElasticMaterial);
                var ma = ma0 as LinearElasticMaterial;
                double lambda = ma.LameLambda;
                double mu = ma.LameMu;
                double rho = ma.MassDensity;

                double[,] sNN = lineFE.CalcSNN();
                double[,] sNyNy = lineFE.CalcSNxNx();
                //double[,] sNyN = lineFE.CalcSNxN();
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
                        double rzzVal = mu * sNN[row, col];
                        double tzzVal = mu * sNyNy[row, col];
                        double mzzVal = rho * sNN[row, col];
                        double sNNVal = sNN[row, col];

                        Rzz[rowNodeId, colNodeId] += rzzVal;
                        Tzz[rowNodeId, colNodeId] += tzzVal;
                        Mzz[rowNodeId, colNodeId] += mzzVal;
                        //
                        SNN[rowNodeId, colNodeId] += sNNVal;
                    }
                }
            }
        }

        public override void Solve()
        {
            _Solve();
            PostSolve();
        }

        private void _Solve()
        {
            Betas = null;
            UEVecs = null;

            // 角周波数
            double omega = 2.0 * Math.PI * Frequency;

            CalcMatrixs();

            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            int uDof = 1;
            var AMat = new IvyFEM.Lapack.DoubleMatrix(nodeCnt * uDof, nodeCnt * uDof);
            var BMat = new IvyFEM.Lapack.DoubleMatrix(nodeCnt * uDof, nodeCnt * uDof);
            for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
            {
                for (int jNodeId = 0; jNodeId < nodeCnt; jNodeId++)
                {
                    //
                    AMat[iNodeId, jNodeId] = Tzz[iNodeId, jNodeId] - omega * omega * Mzz[iNodeId, jNodeId];

                    //
                    BMat[iNodeId, jNodeId] = -1.0 * Rzz[iNodeId, jNodeId];
                }
            }

            uint maxQuantityId = QuantityId;
            bool isDoubleSize = false;
            int portId = (int)PortId; // ポート番号指定
            DoubleSetFixedCadsCondtionForEigen(AMat, BMat, maxQuantityId, isDoubleSize, portId);

            System.Diagnostics.Debug.WriteLine("solve eigen");
            System.Numerics.Complex[] eVals = null;
            System.Numerics.Complex[][] eVecs = null;
            int ret = -1;
            try
            {
                ret = IvyFEM.Lapack.Functions.dggev_dirty(AMat.Buffer, AMat.RowLength, AMat.ColumnLength,
                    BMat.Buffer, BMat.RowLength, BMat.ColumnLength,
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
                int n = AMat.RowLength;
                eVals = new System.Numerics.Complex[n];
                eVecs = new System.Numerics.Complex[n][];
                for (int iMode = 0; iMode < n; iMode++)
                {
                    eVecs[iMode] = new System.Numerics.Complex[n];
                }
            }

            // βに変換
            {
                int modeCnt = eVals.Length;
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    System.Numerics.Complex eVal = eVals[iMode];
                    // λ = β^2と置いた場合(β =√λ)
                    System.Numerics.Complex beta = System.Numerics.Complex.Sqrt(eVal);
                    if (beta.Imaginary > 0)
                    {
                        beta = System.Numerics.Complex.Conjugate(beta);
                    }
                    eVals[iMode] = beta;
                }
            }

            SortEVals(eVals, eVecs);
            System.Numerics.Complex[] propEVals;
            System.Numerics.Complex[][] propEVecs;
            GetPropagationModes(eVals, eVecs, out propEVals, out propEVecs);
            System.Diagnostics.Debug.Assert(propEVals.Length > 0);//DEBUG
            GetBetasUVecs(omega, propEVals, propEVecs);
            AdjustPhaseEVecs(propEVecs);

            Betas = propEVals;
            UEVecs = propEVecs;
        }

        private void SortEVals(System.Numerics.Complex[] betas, System.Numerics.Complex[][] eVecs)
        {
            int modeCnt = betas.Length;
            var eValEVecs = new List<KeyValuePair<System.Numerics.Complex, System.Numerics.Complex[]>>();
            for (int i = 0; i < modeCnt; i++)
            {
                eValEVecs.Add(
                    new KeyValuePair<System.Numerics.Complex, System.Numerics.Complex[]>(betas[i], eVecs[i]));
            }
            eValEVecs.Sort((a, b) =>
            {
                System.Numerics.Complex betaA = a.Key;
                System.Numerics.Complex betaB = b.Key;
                int cmp = 0;
                //double th = 1.0e-12;// これだと小さすぎる場合がある
                double th = 1.0e-10;
                if (Math.Abs(betaA.Real) < th && Math.Abs(betaB.Real) < th)
                {
                    double cmpfi = Math.Abs(betaA.Imaginary) - Math.Abs(betaB.Imaginary);
                    // 昇順
                    if (Math.Abs(cmpfi) < 1.0e-15)
                    {
                        cmp = 0;
                    }
                    else if (cmpfi > 0)
                    {
                        cmp = 1;
                    }
                    else
                    {
                        cmp = -1;
                    }
                }
                else if (Math.Abs(betaA.Real) < th && Math.Abs(betaB.Real) >= th)
                {
                    if (betaB.Real >= 0)
                    {
                        cmp = 1;
                    }
                    else
                    {
                        cmp = -1;
                    }
                }
                else if (Math.Abs(betaA.Real) >= th && Math.Abs(betaB.Real) < th)
                {
                    if (betaA.Real >= 0)
                    {
                        cmp = -1;
                    }
                    else
                    {
                        cmp = 1;
                    }
                }
                else
                {
                    double cmpf = betaA.Real - betaB.Real;
                    if (Math.Abs(cmpf) < 1.0e-15)
                    {
                        cmp = 0;
                    }
                    else
                    {
                        if (cmpf > 0)
                        {
                            cmp = -1;
                        }
                        else
                        {
                            cmp = 1;
                        }
                    }
                }
                return cmp;
            });

            for (int i = 0; i < modeCnt; i++)
            {
                betas[i] = eValEVecs[i].Key;
                eVecs[i] = eValEVecs[i].Value;
            }
        }

        // 伝搬モード一覧の取得
        private void GetPropagationModes(
            System.Numerics.Complex[] betas, System.Numerics.Complex[][] eVecs,
            out System.Numerics.Complex[] propBetas, out System.Numerics.Complex[][] propEVecs)
        {
            propBetas = null;
            propEVecs = null;

            IList<System.Numerics.Complex> retBetas = new List<System.Numerics.Complex>();
            IList<System.Numerics.Complex[]> retEVecs = new List<System.Numerics.Complex[]>();
            int modeCnt = eVecs.Length;
            //double th = 1.0e-12; // これだと小さすぎるときがある
            double th = 1.0e-10;
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                System.Numerics.Complex beta = betas[iMode];
                System.Numerics.Complex[] eVec = eVecs[iMode];
                /*
                if (IsPlaneWave(eVec))
                {
                    // 平面波は除外する
                    System.Diagnostics.Debug.WriteLine("!!!!! skip plane wave. iMode=" + iMode);
                    continue;
                }
                */
                //if (Math.Abs(beta.Imaginary) < th && beta.Real > th) // 正の方向の伝搬モード
                //通常は正負ペアだが、負だけのときがある?
                if (Math.Abs(beta.Imaginary) < th)
                {
                    retBetas.Add(beta);
                    System.Numerics.Complex[] copyEVec = new System.Numerics.Complex[eVec.Length];
                    eVec.CopyTo(copyEVec, 0);
                    retEVecs.Add(copyEVec);
                }
            }
            {
                IList<System.Numerics.Complex> retBetas2 = new List<System.Numerics.Complex>();
                IList<System.Numerics.Complex[]> retEVecs2 = new List<System.Numerics.Complex[]>();

                // 伝搬モードのペアを除外する
                System.Numerics.Complex prevBeta = double.MaxValue; // ありえない値をセット
                for (int iMode = 0; iMode < retBetas.Count; iMode++)
                {
                    System.Numerics.Complex beta = retBetas[iMode];
                    System.Numerics.Complex[] eVec = retEVecs[iMode];
                    if (Math.Abs(Math.Abs(beta.Real) - Math.Abs(prevBeta.Real)) >= 1.0e-12)
                    {
                        // 異なる固有値
                        // OK
                        if (beta.Real < 0)
                        {
                            // 固有値を反転させる
                            beta *= -1.0;
                        }
                        retBetas2.Add(beta);
                        retEVecs2.Add(eVec);
                    }
                    else
                    {
                        // 正負の異なる同じ固有値
                        // skip
                    }
                    prevBeta = beta;
                }
                // 置き換える
                retBetas = retBetas2;
                retEVecs = retEVecs2;
            }
            propBetas = retBetas.ToArray();
            propEVecs = retEVecs.ToArray();
        }

        private bool IsPlaneWave(System.Numerics.Complex[] eVec)
        {
            bool isPlaneWave = true;
            int nodeCnt = eVec.Length;
            System.Numerics.Complex uz0 = eVec[0];
            for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
            {
                System.Numerics.Complex uz = eVec[iNodeId];
                if ((uz - uz0).Magnitude >= 1.0e-12)
                {
                    isPlaneWave = false;
                    break;
                }
            }
            return isPlaneWave;
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
                for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
                {
                    System.Numerics.Complex value = eVec[iNodeId];
                    double abs = value.Magnitude;
                    if (abs > maxAbs)
                    {
                        maxAbs = abs;
                        maxValue = value;
                    }
                }
                System.Numerics.Complex phase =
                    Math.Abs(maxAbs) >= 1.0e-12 ? maxValue / maxAbs : 1.0;

                for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
                {
                    eVec[iNodeId] /= phase;

                    //1.0e-12が小さすぎる場合がある
                    //System.Diagnostics.Debug.Assert(Math.Abs(eVec[iNodeId].Imaginary) < 1.0e-12); // 実数
                    System.Diagnostics.Debug.Assert(Math.Abs(eVec[iNodeId].Imaginary) < 1.0e-10); // 実数
                }
            }
        }

        private void GetBetasUVecs(
            double omega, System.Numerics.Complex[] eVals, System.Numerics.Complex[][] eVecs)
        {
            // 一様媒質とする
            uint maId0;
            {
                IList<uint> feIds = World.GetPortLineFEIds(QuantityId, PortId);
                uint feId0 = feIds[0];
                LineFE lineFE = World.GetLineFE(QuantityId, feId0);
                maId0 = lineFE.MaterialId;
            }
            var ma0 = World.GetMaterial(maId0);
            System.Diagnostics.Debug.Assert(ma0 is LinearElasticMaterial);
            var ma = ma0 as LinearElasticMaterial;
            double rho = ma.MassDensity;
            double lambda = ma.LameLambda;
            double mu = ma.LameMu;

            int modeCnt = eVals.Length;
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var beta = eVals[iMode];

                var sNN = new IvyFEM.Lapack.ComplexMatrix(SNN);

                var eVec = eVecs[iMode];
                int nodeCnt = eVec.Length;
                var uzVec = eVec;
                System.Numerics.Complex power = 0;
                //
                var work1 = sNN * uzVec;
                var work2 = IvyFEM.Lapack.Functions.zdotc(uzVec, work1);
                power += -System.Numerics.Complex.ImaginaryOne * beta * mu * work2;
                //
                power *= System.Numerics.Complex.ImaginaryOne * omega;

                //------------------------------------------
                // 速度ベースで規格化する場合
                //------------------------------------------
                System.Numerics.Complex d =
                    System.Numerics.Complex.Sqrt(
                        omega * omega * System.Numerics.Complex.Conjugate(beta) /
                        (beta.Magnitude * power));
                eVec = IvyFEM.Lapack.Functions.zscal(eVec, (1.0 / omega) * d);
                //------------------------------------------

                eVecs[iMode] = eVec;
            }
        }

        public void CalcBoundaryMatrix(
            double omega, 
            System.Numerics.Complex[] betas,
            System.Numerics.Complex[][] uEVecs,
            System.Numerics.Complex[][] hSigmaEVecs,
            out IvyFEM.Lapack.ComplexMatrix B,
            out IvyFEM.Lapack.ComplexMatrix Fzz)
        {
            // 一様媒質とする
            uint maId0;
            {
                IList<uint> feIds = World.GetPortLineFEIds(QuantityId, PortId);
                uint feId0 = feIds[0];
                LineFE lineFE = World.GetLineFE(QuantityId, feId0);
                maId0 = lineFE.MaterialId;
            }
            var ma0 = World.GetMaterial(maId0);
            System.Diagnostics.Debug.Assert(ma0 is LinearElasticMaterial);
            var ma = ma0 as LinearElasticMaterial;
            double rho = ma.MassDensity;
            double lambda = ma.LameLambda;
            double mu = ma.LameMu;

            int nodeCnt = uEVecs[0].Length;

            int modeCnt = betas.Length;

            B = new Lapack.ComplexMatrix(nodeCnt, nodeCnt);
            Fzz = new Lapack.ComplexMatrix(nodeCnt, nodeCnt);

            var sNN = new IvyFEM.Lapack.ComplexMatrix(SNN);

            //
            B = IvyFEM.Lapack.ComplexMatrix.Scal(sNN, 1.0);

            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var beta = betas[iMode];
                var fVec = uEVecs[iMode];
                var gVec = hSigmaEVecs[iMode];

                for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
                {
                    for (int jNodeId = 0; jNodeId < nodeCnt; jNodeId++)
                    {
                        for (int kNodeId = 0; kNodeId < nodeCnt; kNodeId++)
                        {
                            System.Numerics.Complex sNkNj = sNN[kNodeId, jNodeId];
                            //
                            System.Numerics.Complex fzi = fVec[iNodeId];
                            System.Numerics.Complex fzk = fVec[kNodeId];
                            System.Numerics.Complex gzxi = gVec[iNodeId];
                            System.Numerics.Complex gzxk = gVec[kNodeId];

                            // 速度ベースで規格化する場合
                            //------------------------------------------
                            // hat F
                            Fzz[iNodeId, jNodeId] +=
                                (1.0 / System.Numerics.Complex.ImaginaryOne) *
                                (1.0 / omega) *
                                gzxi * gzxk * sNkNj;
                            //------------------------------------------
                        }
                    }
                }
            }
        }

        public void CalcIncidentVec(
            double omega,
            System.Numerics.Complex beta0,
            System.Numerics.Complex[] uEVec0,
            System.Numerics.Complex[] hSigmaEVec0,
            out System.Numerics.Complex[] Iz)
        {
            int nodeCnt = uEVec0.Length;

            Iz = new System.Numerics.Complex[nodeCnt];

            for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
            {
                System.Numerics.Complex fz = uEVec0[iNodeId];
                System.Numerics.Complex gzx = hSigmaEVec0[iNodeId];

                // 速度ベースで規格化する場合
                //------------------------------------------
                Iz[iNodeId] = (1.0 / System.Numerics.Complex.ImaginaryOne) * 2.0 * gzx;
                //------------------------------------------
            }
        }

        public System.Numerics.Complex[] CalcSMatrix(
            double omega,
            int incidentModeId,
            System.Numerics.Complex[] betas,
            System.Numerics.Complex[][] uEVecs,
            System.Numerics.Complex[][] hSigmaZXEVecs,
            System.Numerics.Complex[] u,
            System.Numerics.Complex[] hSigmaZX)
        {

            // 一様媒質とする
            uint maId0;
            {
                IList<uint> feIds = World.GetPortLineFEIds(QuantityId, PortId);
                uint feId0 = feIds[0];
                LineFE lineFE = World.GetLineFE(QuantityId, feId0);
                maId0 = lineFE.MaterialId;
            }
            var ma0 = World.GetMaterial(maId0);
            System.Diagnostics.Debug.Assert(ma0 is LinearElasticMaterial);
            var ma = ma0 as LinearElasticMaterial;
            double rho = ma.MassDensity;
            double lambda = ma.LameLambda;
            double mu = ma.LameMu;

            int nodeCnt = uEVecs[0].Length;
            System.Diagnostics.Debug.Assert(u.Length == uEVecs[0].Length);
            System.Diagnostics.Debug.Assert(hSigmaZX.Length == hSigmaZXEVecs[0].Length);

            int modeCnt = betas.Length;

            System.Numerics.Complex[] S = new System.Numerics.Complex[modeCnt];

            var sNN = new IvyFEM.Lapack.ComplexMatrix(SNN);

            var _uz = new System.Numerics.Complex[nodeCnt];
            var _conjugateUz = new System.Numerics.Complex[nodeCnt];
            var _hSigmaZX = new System.Numerics.Complex[nodeCnt];
            var _conjugatehSigmaZX = new System.Numerics.Complex[nodeCnt];
            for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
            {
                _uz[iNodeId] = u[iNodeId];
                _conjugateUz[iNodeId] = System.Numerics.Complex.Conjugate(u[iNodeId]);
                _hSigmaZX[iNodeId] = hSigmaZX[iNodeId];
                _conjugatehSigmaZX[iNodeId] = System.Numerics.Complex.Conjugate(hSigmaZX[iNodeId]);
            }

            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var beta = betas[iMode];
                var fVec = uEVecs[iMode];
                var gVec = hSigmaZXEVecs[iMode];

                var fzVec = new System.Numerics.Complex[nodeCnt];
                var gzxVec = new System.Numerics.Complex[nodeCnt];

                for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
                {
                    fzVec[iNodeId] = fVec[iNodeId];
                    gzxVec[iNodeId] = gVec[iNodeId];
                }

                System.Numerics.Complex b = 0;
                // 速度ベースで規格化した場合
                //-----------------------------------
                var work1 = sNN * _uz;
                System.Numerics.Complex work2 =(1.0 / omega) * IvyFEM.Lapack.Functions.zdotu(gzxVec, work1);
                b += work2;
                //-----------------------------------

                //!!!!!!!!!!!!!!
                b = (-1.0 * NormalX) * b;
                //!!!!!!!!!!!!!!
                if (incidentModeId != -1 && incidentModeId == iMode)
                {
                    //---------------------
                    b = (-1.0) + b;
                    //---------------------
                }
                S[iMode] = b;
            }

            return S;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////
        // for FirstOrderABC
        public void CalcBoundaryMatrixForFirstOrderABC(
            double omega,
            System.Numerics.Complex beta0,
            out IvyFEM.Lapack.ComplexMatrix Bzz,
            out IvyFEM.Lapack.ComplexMatrix sNN)
        {
            System.Diagnostics.Debug.Assert(Math.Abs(beta0.Imaginary) < 1.0e-10); // 伝搬モード

            // 一様媒質とする
            uint maId0;
            {
                IList<uint> feIds = World.GetPortLineFEIds(QuantityId, PortId);
                uint feId0 = feIds[0];
                LineFE lineFE = World.GetLineFE(QuantityId, feId0);
                maId0 = lineFE.MaterialId;
            }
            var ma0 = World.GetMaterial(maId0);
            System.Diagnostics.Debug.Assert(ma0 is LinearElasticMaterial);
            var ma = ma0 as LinearElasticMaterial;
            double rho = ma.MassDensity;
            double lambda = ma.LameLambda;
            double mu = ma.LameMu;

            sNN = new IvyFEM.Lapack.ComplexMatrix(SNN);


            // 波の速度
            //-------------------------------------
            //double vp = Math.Sqrt((lambda + 2.0 * mu) / rho);
            double vs = Math.Sqrt(mu / rho);
            //-------------------------------------

            // S waveのABC
            //------------------------------------------------------------------
            //
            System.Numerics.Complex czz =
                -1.0 * System.Numerics.Complex.ImaginaryOne * omega * rho * vs;
            Bzz = IvyFEM.Lapack.ComplexMatrix.Scal(sNN, czz);
            //------------------------------------------------------------------
        }

        /////////////////////////////////////////////////////////////////////////////////
        // ABC / PML 共用
        public System.Numerics.Complex CalcModeAmp(
            double omega,
            System.Numerics.Complex beta,
            System.Numerics.Complex[] uEVec,
            System.Numerics.Complex[] hSigmaZXEVec,
            System.Numerics.Complex[] u,
            System.Numerics.Complex[] hSigmaZX)
        {

            // 一様媒質とする
            uint maId0;
            {
                IList<uint> feIds = World.GetPortLineFEIds(QuantityId, PortId);
                uint feId0 = feIds[0];
                LineFE lineFE = World.GetLineFE(QuantityId, feId0);
                maId0 = lineFE.MaterialId;
            }
            var ma0 = World.GetMaterial(maId0);
            System.Diagnostics.Debug.Assert(ma0 is LinearElasticMaterial);
            var ma = ma0 as LinearElasticMaterial;
            double rho = ma.MassDensity;
            double lambda = ma.LameLambda;
            double mu = ma.LameMu;

            int nodeCnt = uEVec.Length;
            System.Diagnostics.Debug.Assert(u.Length == uEVec.Length);
            System.Diagnostics.Debug.Assert(hSigmaZX.Length == hSigmaZXEVec.Length);

            System.Numerics.Complex amp = 0.0;

            var sNN = new IvyFEM.Lapack.ComplexMatrix(SNN);

            var _uz = new System.Numerics.Complex[nodeCnt];
            var _conjugateUz = new System.Numerics.Complex[nodeCnt];
            var _hSigmaZX = new System.Numerics.Complex[nodeCnt];
            var _conjugatehSigmaZX = new System.Numerics.Complex[nodeCnt];
            for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
            {
                _uz[iNodeId] = u[iNodeId];
                _conjugateUz[iNodeId] = System.Numerics.Complex.Conjugate(u[iNodeId]);
                _hSigmaZX[iNodeId] = hSigmaZX[iNodeId];
                _conjugatehSigmaZX[iNodeId] = System.Numerics.Complex.Conjugate(hSigmaZX[iNodeId]);
            }

            {
                var fVec = uEVec;
                var gVec = hSigmaZXEVec;

                var fzVec = new System.Numerics.Complex[nodeCnt];
                var gzxVec = new System.Numerics.Complex[nodeCnt];

                for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
                {
                    fzVec[iNodeId] = fVec[iNodeId];
                    gzxVec[iNodeId] = gVec[iNodeId];
                }

                System.Numerics.Complex b = 0;
                // 速度ベースで規格化した場合
                //-----------------------------------
                var work1 = sNN * _uz;
                System.Numerics.Complex work2 = (1.0 / omega) * IvyFEM.Lapack.Functions.zdotu(gzxVec, work1);
                b += work2;
                //-----------------------------------

                /*問題による?ので呼び出し側でやる
                //------------------------
                amp = NormalX * b;
                //------------------------
                */
                //OK
                //------------------------
                amp = b;
                //------------------------
            }

            return amp;
        }

        /////////////////////////////////////////////////////////////////////////////////
        // PML
        public void CalcBoundaryMatrixForPML(
            out IvyFEM.Lapack.ComplexMatrix sNN)
        {
            // 一様媒質とする
            uint maId0;
            {
                IList<uint> feIds = World.GetPortLineFEIds(QuantityId, PortId);
                uint feId0 = feIds[0];
                LineFE lineFE = World.GetLineFE(QuantityId, feId0);
                maId0 = lineFE.MaterialId;
            }
            var ma0 = World.GetMaterial(maId0);
            System.Diagnostics.Debug.Assert(ma0 is LinearElasticMaterial);
            var ma = ma0 as LinearElasticMaterial;
            double rho = ma.MassDensity;
            double lambda = ma.LameLambda;
            double mu = ma.LameMu;

            sNN = new IvyFEM.Lapack.ComplexMatrix(SNN);
        }
    }
}
