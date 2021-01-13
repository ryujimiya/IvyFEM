using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Poisson3DFEM : FEM
    {
        // output
        public double[] U { get; private set; } = null;

        public Poisson3DFEM(FEWorld world)
        {
            World = world;
        }

        private void CalcAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint quantityId = 0;
            IList<uint> feIds = World.GetTetrahedronFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = World.GetTetrahedronFE(quantityId, feId);
                uint elemNodeCnt = tetFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = tetFE.NodeCoordIds[iNode];
                    int nodeId = World.Coord2Node(quantityId, coId);
                    nodes[iNode] = nodeId;
                }

                Material ma0 = World.GetMaterial(tetFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is PoissonMaterial);
                var ma = ma0 as PoissonMaterial;
                double k = ma.Alpha;
                double f = ma.F;

                double[,] sNN = tetFE.CalcSNN();
                double[,][,] sNuNv = tetFE.CalcSNuNv();
                double[,] sNxNx = sNuNv[0, 0];
                double[,] sNyNy = sNuNv[1, 1];
                double[,] sNzNz = sNuNv[2, 2];
                double[] sN = tetFE.CalcSN();
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
                        double a = k * (sNxNx[row, col] + sNyNy[row, col] + sNzNz[row, col]);

                        A[rowNodeId, colNodeId] += a;
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
