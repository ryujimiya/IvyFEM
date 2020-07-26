using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class EddyCurrent2DTDFEM : FEM
    {
        public double TimeStep { get; private set; } = 0;
        public double NewmarkBeta { get; private set; } = 1.0 / 4.0;
        public double NewmarkGamma { get; private set; } = 1.0 / 2.0;
        public uint ValueId { get; private set; } = 0;
        public uint PrevValueId { get; private set; } = 0;

        // output
        public double[] U { get; private set; } = null;
        public double[] CoordB { get; protected set; } = null; // 磁束密度

        public EddyCurrent2DTDFEM(FEWorld world, double timeStep,
            double newmarkBeta, double newmarkGamma,
            uint valueId, uint prevValueId)
        {
            World = world;
            TimeStep = timeStep;
            NewmarkBeta = newmarkBeta;
            NewmarkGamma = newmarkGamma;
            ValueId = valueId;
            PrevValueId = prevValueId;
        }

        public void UpdateFieldValuesTimeDomain()
        {
            UpdateFieldValuesNewmarkBetaTimeDomain(
                U, ValueId, PrevValueId,
                TimeStep,
                NewmarkBeta, NewmarkGamma);
        }

        private void CalcAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint quantityId = 0;
            double dt = TimeStep;
            double beta = NewmarkBeta;
            double gamma = NewmarkGamma;
            var FV = World.GetFieldValue(ValueId);
            IList<uint> feIds = World.GetTriangleFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = World.GetTriangleFE(quantityId, feId);
                Material ma0 = World.GetMaterial(triFE.MaterialId);

                int[] coIds = triFE.NodeCoordIds;
                uint elemNodeCnt = triFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = coIds[iNode];
                    int nodeId = World.Coord2Node(quantityId, coId);
                    nodes[iNode] = nodeId;
                }

                System.Diagnostics.Debug.Assert(ma0 is EddyCurrentMaterial);
                var ma = ma0 as EddyCurrentMaterial;
                double mur = ma.Mu;
                double sigma = ma.Sigma;
                double gradPhi = ma.GradPhi;
                double mu = mur * Constants.Mu0;

                double[] sN = triFE.CalcSN();
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

                        int colCoId = coIds[col];
                        double[] u = FV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                        double[] vel = FV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                        double[] acc = FV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);

                        double k = (1.0 / mu) * (sNxNx[row, col] + sNyNy[row, col]);
                        double m = sigma * sNN[row, col];

                        A[rowNodeId, colNodeId] +=
                            k + (gamma / (beta * dt)) * m;
                        B[rowNodeId] +=
                            m * (
                            (gamma / (beta * dt)) * u[0] -
                            (1.0 - gamma / beta) * vel[0] -
                            dt * (1.0 - gamma / (2.0 * beta)) * acc[0]
                            );
                    }
                }

                for (int row = 0; row < elemNodeCnt; row++)
                {
                    int rowNodeId = nodes[row];
                    if (rowNodeId == -1)
                    {
                        continue;
                    }
                    B[rowNodeId] += -sigma *gradPhi * sN[row];
                }
            }
        }

        private void PostSolve()
        {
            // 磁束密度を求める
            CoordB = CalcMagneticDensityValues(U);
        }

        // 磁束密度を計算する
        private double[] CalcMagneticDensityValues(double[] AVec)
        {
            uint quantityId = 0;
            int nodeCnt = (int)World.GetNodeCount(quantityId);
            System.Diagnostics.Debug.Assert(AVec.Length == nodeCnt);
            int coCnt = (int)World.GetCoordCount(quantityId);

            int bDof = 2;
            double[] BVec = new double[coCnt * bDof];

            IList<uint> feIds = World.GetTriangleFEIds(quantityId);

            Dictionary<int, double> coordAreaSum = new Dictionary<int, double>();

            foreach (uint feId in feIds)
            {
                TriangleFE triFE = World.GetTriangleFE(quantityId, feId);
                uint elemNodeCnt = triFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = triFE.NodeCoordIds[iNode];
                    int nodeId = World.Coord2Node(quantityId, coId);
                    nodes[iNode] = nodeId;
                }

                Material ma0 = World.GetMaterial(triFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is EddyCurrentMaterial);
                var ma = ma0 as EddyCurrentMaterial;
                double mur = ma.Mu;
                double sigma = ma.Sigma;
                double gradPhi = ma.GradPhi;
                double mu = mur * Constants.Mu0;

                double A = triFE.GetArea();
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int iNodeId = nodes[iNode];
                    int iCoId = triFE.NodeCoordIds[iNode];
                    double[] L = triFE.GetNodeL(iNode);
                    double[][] Nu = triFE.CalcNu(L);
                    double[] Nx = Nu[0];
                    double[] Ny = Nu[1];

                    double fx = 0.0;
                    double fy = 0.0;
                    for (int kNode = 0; kNode < elemNodeCnt; kNode++)
                    {
                        int kNodeId = nodes[kNode];
                        if (kNodeId == -1)
                        {
                            continue;
                        }
                        double f = AVec[kNodeId];
                        fx += f * Nx[kNode];
                        fy += f * Ny[kNode];
                    }
                    BVec[iCoId * bDof] += fy * A;
                    BVec[iCoId * bDof + 1] += -fx * A;

                    if (coordAreaSum.ContainsKey(iNodeId))
                    {
                        coordAreaSum[iNodeId] += A;
                    }
                    else
                    {
                        coordAreaSum[iNodeId] = A;
                    }
                }
            }

            System.Diagnostics.Debug.Assert(coordAreaSum.Count == coCnt);
            for (int coId = 0; coId < coCnt; coId++)
            {
                double areaSum = coordAreaSum[coId];
                BVec[coId * bDof] /= areaSum;
                BVec[coId * bDof + 1] /= areaSum;
            }
            return BVec;
        }

        public override void Solve()
        {
            uint quantityId = 0;
            int nodeCnt = (int)World.GetNodeCount(quantityId);
            U = new double[nodeCnt];

            var A = new IvyFEM.Linear.DoubleSparseMatrix(nodeCnt, nodeCnt);
            var B = new double[nodeCnt];

            CalcAB(A, B);

            DoubleSetFixedCadsCondtion(A, B);

            double[] X;
            Solver.DoubleSolve(out X, A, B);
            U = X;

            PostSolve();
        }
    }
}
