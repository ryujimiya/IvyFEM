using OpenTK.Graphics.ES11;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class ElasticLambWaveguide1DEigenFEM : FEM
    {
        public uint QuantityId { get; protected set; } = 0;
        public uint PortId { get; protected set; } = 0;

        public IvyFEM.Lapack.DoubleMatrix Rxx { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Ryy { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Sxy { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Syx { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Txx { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Tyy { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Mxx { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Myy { get; protected set; } = null;
        //
        public IvyFEM.Lapack.DoubleMatrix SNN { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix SNNy { get; protected set; } = null;

        // Solve
        // Input
        public double NormalX { get; set; } = 1.0;
        public double Frequency { get; set; }
        // Output
        public System.Numerics.Complex[] Betas { get; protected set; }
        public System.Numerics.Complex[][] HUEVecs { get; protected set; } // ux, uy
        public System.Numerics.Complex[][] HSigmaEVecs { get; protected set; } // (hat)sigmax, sigmay

        public ElasticLambWaveguide1DEigenFEM(FEWorld world, uint quantityId, uint portId)
        {
            World = world;
            QuantityId = quantityId;
            PortId = portId;
        }

        private void PostSolve()
        {
            int modeCnt = Betas.Length;
            // 応力σxx,σyxを求める
            HSigmaEVecs = new System.Numerics.Complex[modeCnt][];
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var beta = Betas[iMode];
                var hUVec = HUEVecs[iMode];
                HSigmaEVecs[iMode] = CalcStressValues(beta, hUVec);
            }
        }

        // 応力σxx,σyxを求める
        private System.Numerics.Complex[] CalcStressValues(System.Numerics.Complex beta, System.Numerics.Complex[] hUVec)
        {
            // 角周波数
            double omega = 2.0 * Math.PI * Frequency;

            uint uDof = 2;
            uint uyOffset = 1;
            uint sDof = 2;
            uint syxOffset = 1;
            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            System.Diagnostics.Debug.Assert(hUVec.Length == nodeCnt * uDof);

            var sigmaVec = new System.Numerics.Complex[nodeCnt * sDof];

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

                    System.Numerics.Complex sigmaXX = 0;
                    System.Numerics.Complex sigmaYX = 0;
                    for (int kNode = 0; kNode < elemNodeCnt; kNode++)
                    {
                        int kNodeId = nodes[kNode];
                        if (kNodeId == -1)
                        {
                            continue;
                        }
                        // hatから元の量に変換
                        // uy = j hat uy
                        System.Numerics.Complex ux = hUVec[kNodeId * uDof];
                        System.Numerics.Complex uy =
                            System.Numerics.Complex.ImaginaryOne * hUVec[kNodeId * uDof + uyOffset];
                        sigmaXX +=
                            -System.Numerics.Complex.ImaginaryOne * beta * (lambda + 2.0 * mu) * N[kNode] * ux +
                            lambda * Ny[kNode] * uy;
                        sigmaYX +=
                            -System.Numerics.Complex.ImaginaryOne * beta * mu * N[kNode] * uy +
                            mu * Ny[kNode] * ux;
                    }
                    sigmaVec[iNodeId * sDof] += sigmaXX * lineLen;
                    sigmaVec[iNodeId * sDof + syxOffset] += sigmaYX * lineLen;

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
                sigmaVec[i * sDof] /= lineLenSum;
                sigmaVec[i * sDof + syxOffset] /= lineLenSum;
            }

            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            // 規格化
            for (int i = 0; i < nodeCnt; i++)
            {
                System.Numerics.Complex sigmaXX = sigmaVec[i * sDof];
                System.Numerics.Complex sigmaYX = sigmaVec[i * sDof + syxOffset];
                //!!!!!!!!!!!!!!!!!!
                /*
                //------------------
                // OKの式
                sigmaXX *= 1.0;
                sigmaYX *= 1.0;
                //------------------
                */
                //------------------
                sigmaXX *= beta.Magnitude;
                sigmaYX *= beta.Magnitude;
                //------------------

                sigmaVec[i * sDof] = sigmaXX;
                sigmaVec[i * sDof + syxOffset] = sigmaYX;
            }
            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            // hatσxx = jσxxに変換
            for (int i = 0; i < nodeCnt; i++)
            {
                // σxx
                sigmaVec[i * uDof] *= System.Numerics.Complex.ImaginaryOne;
            }

            return sigmaVec;
        }

        private void CalcMatrixs()
        {
            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            IList<uint> feIds = World.GetPortLineFEIds(QuantityId, PortId);

            Rxx = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            Ryy = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            Sxy = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            Syx = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            Txx = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            Tyy = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            Mxx = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            Myy = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            //
            SNN = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);
            SNNy = new IvyFEM.Lapack.DoubleMatrix(nodeCnt, nodeCnt);

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
                double[,] sNyN = lineFE.CalcSNxN();
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
                        double rxxVal = (lambda + 2.0 * mu) * sNN[row, col];
                        double txxVal = mu * sNyNy[row, col];
                        double sxyVal = lambda * sNyN[col, row] - mu * sNyN[row, col]; // sNNy sNyN
                        double mxxVal = rho * sNN[row, col];
                        double syxVal = -lambda * sNyN[row, col] + mu * sNyN[col, row]; // sNyN sNNy
                        double tyyVal = (lambda + 2.0 * mu) * sNyNy[row, col];
                        double ryyVal = mu * sNN[row, col];
                        double myyVal = rho * sNN[row, col];
                        //
                        double sNNVal = sNN[row, col];
                        double sNNyVal = sNyN[col, row]; // sNNy

                        Rxx[rowNodeId, colNodeId] += rxxVal;
                        Ryy[rowNodeId, colNodeId] += ryyVal;
                        Sxy[rowNodeId, colNodeId] += sxyVal;
                        Syx[rowNodeId, colNodeId] += syxVal;
                        Txx[rowNodeId, colNodeId] += txxVal;
                        Tyy[rowNodeId, colNodeId] += tyyVal;
                        Mxx[rowNodeId, colNodeId] += mxxVal;
                        Myy[rowNodeId, colNodeId] += myyVal;
                        //
                        SNN[rowNodeId, colNodeId] += sNNVal;
                        SNNy[rowNodeId, colNodeId] += sNNyVal;
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
            HUEVecs = null;

            // 角周波数
            double omega = 2.0 * Math.PI * Frequency;

            CalcMatrixs();

            int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);
            int uDof = 2;
            int uyOffset = 1;
            var AMat = new IvyFEM.Lapack.DoubleMatrix(nodeCnt * uDof, nodeCnt * uDof);
            var BMat = new IvyFEM.Lapack.DoubleMatrix(nodeCnt * uDof, nodeCnt * uDof);
            var CMat = new IvyFEM.Lapack.DoubleMatrix(nodeCnt * uDof, nodeCnt * uDof);
            for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
            {
                for (int jNodeId = 0; jNodeId < nodeCnt; jNodeId++)
                {
                    //
                    AMat[iNodeId * uDof, jNodeId * uDof] = -Rxx[iNodeId, jNodeId];
                    AMat[iNodeId * uDof + uyOffset, jNodeId * uDof + uyOffset] = -Ryy[iNodeId, jNodeId];

                    //
                    BMat[iNodeId * uDof, jNodeId * uDof + uyOffset] = Sxy[iNodeId, jNodeId];
                    BMat[iNodeId * uDof + uyOffset, jNodeId * uDof] = Syx[iNodeId, jNodeId];

                    //
                    CMat[iNodeId * uDof, jNodeId * uDof] = Txx[iNodeId, jNodeId] - omega * omega * Mxx[iNodeId, jNodeId];
                    CMat[iNodeId * uDof + uyOffset, jNodeId * uDof + uyOffset] =
                        Tyy[iNodeId, jNodeId] - omega * omega * Myy[iNodeId, jNodeId];
                }
            }

            // 非線形固有値問題
            //  [C] + λ[B] + λ^2[A]{Φ}= {0}

            // [A]の逆行列が存在する場合のみ有効な方法
            System.Diagnostics.Debug.WriteLine("calc [A]-1");
            // [A]の逆行列を求める
            var invAMat = IvyFEM.Lapack.DoubleMatrix.Inverse(AMat);
            System.Diagnostics.Debug.WriteLine("calc [A]-1[C]");
            // [A]-1[C]
            var invACMat = invAMat * CMat;
            System.Diagnostics.Debug.WriteLine("calc [A]-1[B]");
            // [A]-1[B]
            var invABMat = invAMat * BMat;

            int nodeCnt1 = AMat.RowLength;
            var AToSolve = new IvyFEM.Lapack.DoubleMatrix(nodeCnt1 * 2, nodeCnt1 * 2);
            System.Diagnostics.Debug.WriteLine("set [A]");
            for (int i = 0; i < nodeCnt1; i++)
            {
                for (int j = 0; j < nodeCnt1; j++)
                {
                    AToSolve[i, j] = 0.0;
                    AToSolve[i, j + nodeCnt1] = i == j ? 1.0 : 0.0;
                    //  { [K] - jβ[C] - β^2[M] }{Φ}= {0}
                    // λ = -jβと置いた場合
                    AToSolve[i + nodeCnt1, j] = -1.0 * invACMat[i, j];
                    AToSolve[i + nodeCnt1, j + nodeCnt1] = -1.0 * invABMat[i, j];
                }
            }

            uint maxQuantityId = QuantityId;
            bool isDoubleSize = true; // 2倍のサイズ
            int portId = (int)PortId; // ポート番号指定
            DoubleSetFixedCadsCondtionForEigen(AToSolve, null, maxQuantityId, isDoubleSize, portId);

            System.Diagnostics.Debug.WriteLine("solve eigen");
            System.Numerics.Complex[] eVals = null;
            System.Numerics.Complex[][] eVecs = null;
            {
                // 固有値問題を解く
                System.Numerics.Complex[][] eVecs0 = null;
                int ret = -1;
                try
                {
                    /*
                    ret = IvyFEM.Lapack.Functions.dgeev(
                        AToSolve.Buffer, AToSolve.RowLength, AToSolve.ColumnLength,
                        out eVals, out eVecs0);
                    */
                    ret = IvyFEM.Lapack.Functions.dgeev_dirty(
                        AToSolve.Buffer, AToSolve.RowLength, AToSolve.ColumnLength,
                        out eVals, out eVecs0);
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
                    int n = AToSolve.RowLength;
                    eVals = new System.Numerics.Complex[n];
                    eVecs0 = new System.Numerics.Complex[n][];
                    for (int iMode = 0; iMode < n; iMode++)
                    {
                        eVecs0[iMode] = new System.Numerics.Complex[n];
                    }
                }

                // 縮小前の固有ベクトルを求める
                int modeCnt = eVecs0.Length;
                eVecs = new System.Numerics.Complex[modeCnt][];
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    System.Numerics.Complex[] eVec0 = eVecs0[iMode];
                    System.Numerics.Complex[] eVec = new System.Numerics.Complex[nodeCnt * uDof];
                    eVecs[iMode] = eVec;
                    for (int i = 0; i < nodeCnt1; i++)
                    {
                        eVec[i] = eVec0[i];
                    }
                }
            }

            // βに変換
            {
                int modeCnt = eVals.Length;
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    System.Numerics.Complex eVal = eVals[iMode];
                    // λ = αと置いた場合(β = -jα= -jλ)
                    System.Numerics.Complex beta = -1.0 * System.Numerics.Complex.ImaginaryOne * eVal;
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
            // hat uy = -juy に変換
            {
                int modeCnt = propEVals.Length;
                for (int iMode = 0; iMode < modeCnt; iMode++)
                {
                    System.Numerics.Complex[] eVec = propEVecs[iMode];
                    int n = eVec.Length / uDof;
                    for (int iNodeId = 0; iNodeId < n; iNodeId++)
                    {
                        // uy
                        eVec[iNodeId * uDof + 1] *= -1.0 * System.Numerics.Complex.ImaginaryOne;
                    }
                }
            }

            Betas = propEVals;
            HUEVecs = propEVecs;
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
                if (IsPlaneWave(eVec))
                {
                    // 平面波は除外する
                    System.Diagnostics.Debug.WriteLine("!!!!! skip plane wave. iMode=" + iMode);
                    continue;
                }
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
            int uDof = 2;
            int nodeCnt = eVec.Length / uDof;
            System.Numerics.Complex ux0 = eVec[0];
            System.Numerics.Complex uy0 = eVec[1];
            for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
            {
                System.Numerics.Complex ux = eVec[iNodeId * uDof];
                System.Numerics.Complex uy = eVec[iNodeId * uDof + 1];
                if ((ux - ux0).Magnitude >= 1.0e-12)
                {
                    isPlaneWave = false;
                    break;
                }
                if ((uy - uy0).Magnitude >= 1.0e-12)
                {
                    isPlaneWave = false;
                    break;
                }
            }
            return isPlaneWave;
        }

        private void AdjustPhaseEVecs(System.Numerics.Complex[][] eVecs)
        {
            int uDof = 2;
            int uyOffset = 1;
            int modeCnt = eVecs.Length;
            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var eVec = eVecs[iMode];
                int nodeCnt = eVec.Length / uDof;
                System.Numerics.Complex maxValue = new System.Numerics.Complex(0, 0);
                double maxAbs = 0;
                for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
                {
                    System.Numerics.Complex valueX = eVec[iNodeId * uDof];
                    System.Numerics.Complex valueY = eVec[iNodeId * uDof + uyOffset];
                    // uxの位相を基準にする
                    double abs = valueX.Magnitude;
                    if (abs > maxAbs)
                    {
                        maxAbs = abs;
                        maxValue = valueX;
                    }
                }
                System.Numerics.Complex phase =
                    Math.Abs(maxAbs) >= 1.0e-12 ? maxValue / maxAbs : 1.0;

                for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
                {
                    eVec[iNodeId * uDof] /= phase;
                    eVec[iNodeId * uDof + uyOffset] /= phase;

                    //1.0e-12が小さすぎる場合がある
                    //System.Diagnostics.Debug.Assert(Math.Abs(eVec[iNodeId * uDof].Imaginary) < 1.0e-12); // 実数
                    //System.Diagnostics.Debug.Assert(Math.Abs(eVec[iNodeId * uDof + uyOffset].Real) < 1.0e-12); // 純虚数
                    System.Diagnostics.Debug.Assert(Math.Abs(eVec[iNodeId * uDof].Imaginary) < 1.0e-10); // 実数
                    System.Diagnostics.Debug.Assert(Math.Abs(eVec[iNodeId * uDof + uyOffset].Real) < 1.0e-10); // 純虚数
                }
            }
        }

        private void GetBetasUVecs(
            double omega, System.Numerics.Complex[] eVals, System.Numerics.Complex[][] eVecs)
        {
            int uDof = 2;
            int uyOffset = 1;

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
                var sNNy = new IvyFEM.Lapack.ComplexMatrix(SNNy);

                var eVec = eVecs[iMode];
                int nodeCnt = eVec.Length / uDof;
                var uxVec = new System.Numerics.Complex[nodeCnt];
                var uyVec = new System.Numerics.Complex[nodeCnt];
                for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
                {
                    uxVec[iNodeId] = eVec[iNodeId * uDof];
                    uyVec[iNodeId] = eVec[iNodeId * uDof + uyOffset];
                }
                System.Numerics.Complex power = 0;
                //
                var work1 = sNN * uxVec;
                var work2 = IvyFEM.Lapack.Functions.zdotc(uxVec, work1);
                power += -System.Numerics.Complex.ImaginaryOne * beta * (lambda + 2.0 * mu) * work2;
                //
                var work3 = sNNy * uyVec;
                var work4 = IvyFEM.Lapack.Functions.zdotc(uxVec, work3);
                power += lambda * work4;
                //
                var work5 = sNN * uyVec;
                var work6 = IvyFEM.Lapack.Functions.zdotc(uyVec, work5);
                power += -System.Numerics.Complex.ImaginaryOne * beta * mu * work6;
                //
                var work7 = sNNy * uxVec;
                var work8 = IvyFEM.Lapack.Functions.zdotc(uyVec, work7);
                power += mu * work8;
                //
                power *= System.Numerics.Complex.ImaginaryOne * omega;

                //------------------------------------------
                // 速度ベースで規格化する場合
                /*
                //------------------------------------------
                // OKの式
                System.Numerics.Complex d =
                    System.Numerics.Complex.Sqrt(
                        Math.Sqrt(beta.Magnitude) * System.Numerics.Complex.Conjugate(beta) /
                        (omega * beta.Magnitude * power));
                eVec = IvyFEM.Lapack.Functions.zscal(eVec, 1.0 * d);
                //------------------------------------------
                */
                //------------------------------------------
                System.Numerics.Complex d =
                    System.Numerics.Complex.Sqrt(
                        omega * omega * omega * System.Numerics.Complex.Conjugate(beta) /
                        (beta.Magnitude * beta.Magnitude * power));
                eVec = IvyFEM.Lapack.Functions.zscal(eVec, (1.0 / omega) * d);
                //------------------------------------------

                //------------------------------------------
                /*
                // uとσで同じ規格化をする場合
                System.Numerics.Complex d =
                    System.Numerics.Complex.Sqrt(
                        System.Numerics.Complex.ImaginaryOne * omega * System.Numerics.Complex.Conjugate(beta) /
                        (beta.Magnitude * power));

                eVec = IvyFEM.Lapack.Functions.zscal(eVec, d);
                */
                //------------------------------------------

                eVecs[iMode] = eVec;
            }
        }

        public void CalcBoundaryMatrix(
            double omega, 
            System.Numerics.Complex[] betas,
            System.Numerics.Complex[][] hUEVecs,
            System.Numerics.Complex[][] hSigmaEVecs,
            out IvyFEM.Lapack.ComplexMatrix B,
            out IvyFEM.Lapack.ComplexMatrix Fxx,
            out IvyFEM.Lapack.ComplexMatrix Fxy,
            out IvyFEM.Lapack.ComplexMatrix Fyx,
            out IvyFEM.Lapack.ComplexMatrix Fyy)
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

            int uDof = 2;
            int uyOffset = 1;
            int sDof = 2;
            int syxOffset = 1;
            int nodeCnt = hUEVecs[0].Length / uDof;

            int modeCnt = betas.Length;

            //B = new Lapack.ComplexMatrix(nodeCnt, nodeCnt);
            Fxx = new Lapack.ComplexMatrix(nodeCnt, nodeCnt);
            Fxy = new Lapack.ComplexMatrix(nodeCnt, nodeCnt);
            Fyx = new Lapack.ComplexMatrix(nodeCnt, nodeCnt);
            Fyy = new Lapack.ComplexMatrix(nodeCnt, nodeCnt);

            var sNN = new IvyFEM.Lapack.ComplexMatrix(SNN);

            //
            B = IvyFEM.Lapack.ComplexMatrix.Scal(sNN, 1.0);

            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var beta = betas[iMode];
                var fVec = hUEVecs[iMode];
                var gVec = hSigmaEVecs[iMode];

                for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
                {
                    for (int jNodeId = 0; jNodeId < nodeCnt; jNodeId++)
                    {
                        for (int kNodeId = 0; kNodeId < nodeCnt; kNodeId++)
                        {
                            System.Numerics.Complex sNkNj = sNN[kNodeId, jNodeId];
                            //
                            System.Numerics.Complex fxi = fVec[iNodeId * uDof];
                            System.Numerics.Complex fxk = fVec[kNodeId * uDof];
                            System.Numerics.Complex gyxi = gVec[iNodeId * sDof + syxOffset];
                            System.Numerics.Complex gyxk = gVec[kNodeId * sDof + syxOffset];

                            // 速度ベースで規格化する場合
                            /* OKの式
                            //--------------------------------------------
                            // hat F
                            Fxx[iNodeId, jNodeId] +=
                                System.Numerics.Complex.ImaginaryOne *
                                1.0 *
                                fxi * fxk * sNkNj;
                            Fxy[iNodeId, jNodeId] +=
                                -1.0 * System.Numerics.Complex.ImaginaryOne *
                                1.0 *
                                fxi * gyxk * sNkNj;
                            Fyx[iNodeId, jNodeId] +=
                                System.Numerics.Complex.ImaginaryOne *
                                1.0 *
                                gyxi * fxk * sNkNj;
                            Fyy[iNodeId, jNodeId] +=
                                -1.0 * System.Numerics.Complex.ImaginaryOne *
                                1.0 *
                                gyxi * gyxk * sNkNj;
                            //--------------------------------------------
                            */
                            //--------------------------------------------
                            // hat F
                            Fxx[iNodeId, jNodeId] +=
                                System.Numerics.Complex.ImaginaryOne *
                                1.0 *
                                fxi * fxk * sNkNj;
                            Fxy[iNodeId, jNodeId] +=
                                -1.0 * System.Numerics.Complex.ImaginaryOne *
                                (1.0 / omega) *
                                fxi * gyxk * sNkNj;
                            Fyx[iNodeId, jNodeId] +=
                                System.Numerics.Complex.ImaginaryOne *
                                1.0 *
                                gyxi * fxk * sNkNj;
                            Fyy[iNodeId, jNodeId] +=
                                -1.0 * System.Numerics.Complex.ImaginaryOne *
                                (1.0 / omega) *
                                gyxi * gyxk * sNkNj;
                            //--------------------------------------------
                        }
                    }
                }
            }
        }

        public void CalcIncidentVec(
            double omega,
            System.Numerics.Complex beta0,
            System.Numerics.Complex[] hUEVec0,
            System.Numerics.Complex[] hSigmaEVec0,
            out System.Numerics.Complex[] Ix,
            out System.Numerics.Complex[] Iy)
        {
            int uDof = 2;
            int uyOffset = 1;
            int sDof = 2;
            int syxOffset = 1;
            int nodeCnt = hUEVec0.Length / uDof;

            Ix = new System.Numerics.Complex[nodeCnt];
            Iy = new System.Numerics.Complex[nodeCnt];

            for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
            {
                System.Numerics.Complex fx = hUEVec0[iNodeId * uDof];
                System.Numerics.Complex gyx = hSigmaEVec0[iNodeId * sDof + syxOffset];

                // 速度ベースで規格化する場合
                /*
                // OKの式
                //------------------------------
                Ix[iNodeId] = 2.0 * fx;
                Iy[iNodeId] = 2.0 * gyx;
                //------------------------------
                */
                //------------------------------
                Ix[iNodeId] = 2.0 * fx;
                Iy[iNodeId] = 2.0 * gyx;
                //------------------------------
            }
        }

        public System.Numerics.Complex[] CalcSMatrix(
            double omega,
            int incidentModeId,
            System.Numerics.Complex[] betas,
            System.Numerics.Complex[][] hUEVecs,
            System.Numerics.Complex[][] hSigmaEVecs,
            System.Numerics.Complex[] hU,
            System.Numerics.Complex[] hSigma)
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

            int uDof = 2;
            int uyOffset = 1;
            int sDof = 2;
            int syxOffset = 1;
            int nodeCnt = hUEVecs[0].Length / uDof;
            System.Diagnostics.Debug.Assert(hU.Length == hUEVecs[0].Length);
            System.Diagnostics.Debug.Assert(hSigma.Length == hSigmaEVecs[0].Length);

            int modeCnt = betas.Length;

            System.Numerics.Complex[] S = new System.Numerics.Complex[modeCnt];

            var sNN = new IvyFEM.Lapack.ComplexMatrix(SNN);

            var hUx = new System.Numerics.Complex[nodeCnt];
            var hUy = new System.Numerics.Complex[nodeCnt];
            var conjugatehUx = new System.Numerics.Complex[nodeCnt];
            var conjugatehUy = new System.Numerics.Complex[nodeCnt];
            var hSigmaXX = new System.Numerics.Complex[nodeCnt];
            var hSigmaYX = new System.Numerics.Complex[nodeCnt];
            var conjugatehSigmaXX = new System.Numerics.Complex[nodeCnt];
            var conjugatehSigmaYX = new System.Numerics.Complex[nodeCnt];
            for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
            {
                hUx[iNodeId] = hU[iNodeId * uDof];
                hUy[iNodeId] = hU[iNodeId * uDof + uyOffset];
                conjugatehUx[iNodeId] = System.Numerics.Complex.Conjugate(hU[iNodeId * uDof]);
                conjugatehUy[iNodeId] = System.Numerics.Complex.Conjugate(hU[iNodeId * uDof + uyOffset]);
                hSigmaXX[iNodeId] = hSigma[iNodeId * sDof];
                hSigmaYX[iNodeId] = hSigma[iNodeId * sDof + syxOffset];
                conjugatehSigmaXX[iNodeId] = System.Numerics.Complex.Conjugate(hSigma[iNodeId * sDof]);
                conjugatehSigmaYX[iNodeId] = System.Numerics.Complex.Conjugate(hSigma[iNodeId * sDof + syxOffset]);
            }

            for (int iMode = 0; iMode < modeCnt; iMode++)
            {
                var beta = betas[iMode];
                var fVec = hUEVecs[iMode];
                var gVec = hSigmaEVecs[iMode];

                var fxVec = new System.Numerics.Complex[nodeCnt];
                var fyVec = new System.Numerics.Complex[nodeCnt];
                var gxxVec = new System.Numerics.Complex[nodeCnt];
                var gyxVec = new System.Numerics.Complex[nodeCnt];

                for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
                {
                    fxVec[iNodeId] = fVec[iNodeId * uDof];
                    fyVec[iNodeId] = fVec[iNodeId * uDof + uyOffset];
                    gxxVec[iNodeId] = gVec[iNodeId * uDof];
                    gyxVec[iNodeId] = gVec[iNodeId * uDof + syxOffset];
                }

                System.Numerics.Complex b = 0;
                // 速度ベースで規格化した場合
                /*OKの式
                /////////////////////
                //-----------------------------------
                var work1 = sNN * hSigmaXX;
                System.Numerics.Complex work2 = IvyFEM.Lapack.Functions.zdotu(fxVec, work1);
                b += work2;
                //-----------------------------------

                //-----------------------------------
                var work3 = sNN * hUy;
                System.Numerics.Complex work4 =IvyFEM.Lapack.Functions.zdotu(gyxVec, work3);
                b += work4;
                //-----------------------------------
                /////////////////////
                */

                /////////////////////
                //-----------------------------------
                var work1 = sNN * hSigmaXX;
                System.Numerics.Complex work2 = IvyFEM.Lapack.Functions.zdotu(fxVec, work1);
                b += work2;
                //-----------------------------------

                //-----------------------------------
                var work3 = sNN * hUy;
                System.Numerics.Complex work4 = (1.0 / omega) * IvyFEM.Lapack.Functions.zdotu(gyxVec, work3);
                b += work4;
                //-----------------------------------
                /////////////////////

                //!!!!!!!!!!!!!!
                b = NormalX * b;
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
            out IvyFEM.Lapack.ComplexMatrix Bxx,
            out IvyFEM.Lapack.ComplexMatrix Bxy,
            out IvyFEM.Lapack.ComplexMatrix Byx,
            out IvyFEM.Lapack.ComplexMatrix Byy,
            out IvyFEM.Lapack.ComplexMatrix sNN,
            out IvyFEM.Lapack.ComplexMatrix sNNy)
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

            //int uDof = 2;
            //int uyOffset = 1;
            //int sDof = 2;
            //int syxOffset = 1;
            //int nodeCnt = (int)World.GetPortNodeCount(QuantityId, PortId);

            //Bxx = new Lapack.ComplexMatrix(nodeCnt, nodeCnt);
            //Byy = new Lapack.ComplexMatrix(nodeCnt, nodeCnt);

            sNN = new IvyFEM.Lapack.ComplexMatrix(SNN);

            sNNy = new IvyFEM.Lapack.ComplexMatrix(SNNy);

            // 波の速度
            //-------------------------------------
            double vp = Math.Sqrt((lambda + 2.0 * mu) / rho);
            double vs = Math.Sqrt(mu / rho);
            //-------------------------------------
            /*
            //-------------------------------------
            double vp = (beta0.Real / omega) * ((lambda + 2.0 * mu) / rho);
            double vs = (beta0.Real / omega) * (mu / rho);
            //-------------------------------------
            */

            // P-S waveのABC
            //------------------------------------------------------------------
            //
            System.Numerics.Complex cxx =
                -1.0 * System.Numerics.Complex.ImaginaryOne * omega * rho * vp;
            Bxx = IvyFEM.Lapack.ComplexMatrix.Scal(sNN, cxx);
            //
            Bxy = new IvyFEM.Lapack.ComplexMatrix(sNNy.RowLength, sNNy.ColumnLength);
            //
            Byx = new IvyFEM.Lapack.ComplexMatrix(sNNy.RowLength, sNNy.ColumnLength);
            //
            System.Numerics.Complex cyy =
                -1.0 * System.Numerics.Complex.ImaginaryOne * omega * rho * vs;
            Byy = IvyFEM.Lapack.ComplexMatrix.Scal(sNN, cyy);
            //------------------------------------------------------------------
        }

        public System.Numerics.Complex CalcModeAmp(
            double omega,
            System.Numerics.Complex beta,
            System.Numerics.Complex[] hUEVec,
            System.Numerics.Complex[] hSigmaEVec,
            System.Numerics.Complex[] hU,
            System.Numerics.Complex[] hSigma)
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

            int uDof = 2;
            int uyOffset = 1;
            int sDof = 2;
            int syxOffset = 1;
            int nodeCnt = hUEVec.Length / uDof;
            System.Diagnostics.Debug.Assert(hU.Length == hUEVec.Length);
            System.Diagnostics.Debug.Assert(hSigma.Length == hSigmaEVec.Length);

            System.Numerics.Complex amp = 0.0;

            var sNN = new IvyFEM.Lapack.ComplexMatrix(SNN);

            var hUx = new System.Numerics.Complex[nodeCnt];
            var hUy = new System.Numerics.Complex[nodeCnt];
            var conjugatehUx = new System.Numerics.Complex[nodeCnt];
            var conjugatehUy = new System.Numerics.Complex[nodeCnt];
            var hSigmaXX = new System.Numerics.Complex[nodeCnt];
            var hSigmaYX = new System.Numerics.Complex[nodeCnt];
            var conjugatehSigmaXX = new System.Numerics.Complex[nodeCnt];
            var conjugatehSigmaYX = new System.Numerics.Complex[nodeCnt];
            for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
            {
                hUx[iNodeId] = hU[iNodeId * uDof];
                hUy[iNodeId] = hU[iNodeId * uDof + uyOffset];
                conjugatehUx[iNodeId] = System.Numerics.Complex.Conjugate(hU[iNodeId * uDof]);
                conjugatehUy[iNodeId] = System.Numerics.Complex.Conjugate(hU[iNodeId * uDof + uyOffset]);
                hSigmaXX[iNodeId] = hSigma[iNodeId * sDof];
                hSigmaYX[iNodeId] = hSigma[iNodeId * sDof + syxOffset];
                conjugatehSigmaXX[iNodeId] = System.Numerics.Complex.Conjugate(hSigma[iNodeId * sDof]);
                conjugatehSigmaYX[iNodeId] = System.Numerics.Complex.Conjugate(hSigma[iNodeId * sDof + syxOffset]);
            }

            {
                var fVec = hUEVec;
                var gVec = hSigmaEVec;

                var fxVec = new System.Numerics.Complex[nodeCnt];
                var fyVec = new System.Numerics.Complex[nodeCnt];
                var gxxVec = new System.Numerics.Complex[nodeCnt];
                var gyxVec = new System.Numerics.Complex[nodeCnt];

                for (int iNodeId = 0; iNodeId < nodeCnt; iNodeId++)
                {
                    fxVec[iNodeId] = fVec[iNodeId * uDof];
                    fyVec[iNodeId] = fVec[iNodeId * uDof + uyOffset];
                    gxxVec[iNodeId] = gVec[iNodeId * uDof];
                    gyxVec[iNodeId] = gVec[iNodeId * uDof + syxOffset];
                }

                System.Numerics.Complex b = 0;
                // 速度ベースで規格化した場合
                /*
                // 元の式
                /////////////////////
                //-----------------------------------
                var work1 = sNN * hSigmaXX;
                System.Numerics.Complex work2 = IvyFEM.Lapack.Functions.zdotu(fxVec, work1);
                b += work2;
                //-----------------------------------

                //-----------------------------------
                var work3 = sNN * hUy;
                System.Numerics.Complex work4 = (1.0 / omega) * IvyFEM.Lapack.Functions.zdotu(gyxVec, work3);
                b += work4;
                //-----------------------------------
                /////////////////////
                */
                /*
                //hUyにconjugateを用いる
                /////////////////////
                //-----------------------------------
                var work1 = sNN * hSigmaXX;
                System.Numerics.Complex work2 = IvyFEM.Lapack.Functions.zdotu(fxVec, work1);
                b += work2;
                //-----------------------------------

                //-----------------------------------
                var work3 = sNN * conjugatehUy;
                System.Numerics.Complex work4 = (1.0 / omega) * IvyFEM.Lapack.Functions.zdotu(gyxVec, work3);
                b += work4;
                //-----------------------------------
                /////////////////////
                */
                /////////////////////
                //-----------------------------------
                var work1 = sNN * hSigmaXX;
                System.Numerics.Complex work2 = IvyFEM.Lapack.Functions.zdotu(fxVec, work1);
                b += work2;
                //-----------------------------------

                //-----------------------------------
                var work3 = sNN * hUy;
                System.Numerics.Complex work4 = (1.0 / omega) * IvyFEM.Lapack.Functions.zdotu(gyxVec, work3);
                b += work4;
                //-----------------------------------
                /////////////////////

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
            out IvyFEM.Lapack.ComplexMatrix sNN,
            out IvyFEM.Lapack.ComplexMatrix sNNy)
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

            sNNy = new IvyFEM.Lapack.ComplexMatrix(SNNy);
        }
    }
}
