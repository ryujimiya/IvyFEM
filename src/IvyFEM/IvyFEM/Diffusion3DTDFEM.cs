using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Diffusion3DTDFEM : FEM
    {
        public double TimeStep { get; private set; } = 0;
        public double NewmarkBeta { get; private set; } = 1.0 / 4.0;
        public double NewmarkGamma { get; private set; } = 1.0 / 2.0;
        public uint ValueId { get; private set; } = 0;
        public uint PrevValueId { get; private set; } = 0;

        // output
        public double[] U { get; private set; } = null;

        public Diffusion3DTDFEM(FEWorld world, double timeStep,
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
            IList<uint> feIds = World.GetTetrahedronFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = World.GetTetrahedronFE(quantityId, feId);
                Material ma0 = World.GetMaterial(tetFE.MaterialId);

                int[] coIds = tetFE.NodeCoordIds;
                uint elemNodeCnt = tetFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = coIds[iNode];
                    int nodeId = World.Coord2Node(quantityId, coId);
                    nodes[iNode] = nodeId;
                }

                System.Diagnostics.Debug.Assert(ma0 is DiffusionMaterial);
                var ma = ma0 as DiffusionMaterial;
                double rho = ma.MassDensity;
                double cap = ma.Capacity;
                double lambda = ma.DiffusionCoef;
                double f = ma.F;

                double[] sN = tetFE.CalcSN();
                double[,] sNN = tetFE.CalcSNN();
                double[,][,] sNuNv = tetFE.CalcSNuNv();
                double[,] sNxNx = sNuNv[0, 0];
                double[,] sNyNy = sNuNv[1, 1];
                double[,] sNzNz = sNuNv[2, 2];

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

                        double k = lambda * (sNxNx[row, col] + sNyNy[row, col] + sNzNz[row, col]);
                        double m = rho * cap * sNN[row, col];

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
                    B[rowNodeId] += f * sN[row];
                }
            }
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
        }
    }
}
