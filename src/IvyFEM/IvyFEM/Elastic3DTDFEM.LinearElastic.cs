using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic3DTDFEM
    {
        protected void CalcLinearElasticElementAB(
            uint feId, IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            {
                uint quantityId0 = 0;
                TetrahedronFE workTetFE = World.GetTetrahedronFE(quantityId0, feId);
                Material workMa0 = World.GetMaterial(workTetFE.MaterialId);
                if (!(workMa0 is LinearElasticMaterial))
                {
                    return;
                }
            }

            uint quantityId = 0;
            int nodeCnt = (int)World.GetNodeCount(quantityId);
            System.Diagnostics.Debug.Assert(World.GetDof(quantityId) == 3);
            int dof = 3;

            double dt = TimeStep;
            double beta = NewmarkBeta;
            double gamma = NewmarkGamma;
            var FV = World.GetFieldValue(UValueId);

            TetrahedronFE tetFE = World.GetTetrahedronFE(quantityId, feId);
            Material ma0 = World.GetMaterial(tetFE.MaterialId);
            System.Diagnostics.Debug.Assert(ma0 is LinearElasticMaterial);

            int[] coIds = tetFE.NodeCoordIds;
            uint elemNodeCnt = tetFE.NodeCount;
            int[] nodes = new int[elemNodeCnt];
            for (int iNode = 0; iNode < elemNodeCnt; iNode++)
            {
                int coId = coIds[iNode];
                int nodeId = World.Coord2Node(quantityId, coId);
                nodes[iNode] = nodeId;
            }

            var ma = ma0 as LinearElasticMaterial;
            double lambda = ma.LameLambda;
            double mu = ma.LameMu;
            double rho = ma.MassDensity;
            double[] g = { ma.GravityX, ma.GravityY, ma.GravityZ };

            double[] sN = tetFE.CalcSN();
            double[,] sNN = tetFE.CalcSNN();
            double[,][,] sNuNv = tetFE.CalcSNuNv();
            double[,] sNxNx = sNuNv[0, 0];
            double[,] sNxNy = sNuNv[0, 1];
            double[,] sNxNz = sNuNv[0, 2];
            double[,] sNyNx = sNuNv[1, 0];
            double[,] sNyNy = sNuNv[1, 1];
            double[,] sNyNz = sNuNv[1, 2];
            double[,] sNzNx = sNuNv[2, 0];
            double[,] sNzNy = sNuNv[2, 1];
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

                    double[,] k = new double[dof, dof];
                    double[,] m = new double[dof, dof];
                    k[0, 0] = (lambda + mu) * sNxNx[row, col] +
                        mu * (sNxNx[row, col] + sNyNy[row, col] + sNzNz[row, col]);
                    k[0, 1] = lambda * sNxNy[row, col] + mu * sNyNx[row, col];
                    k[0, 2] = lambda * sNxNz[row, col] + mu * sNzNx[row, col];
                    k[1, 0] = lambda * sNyNx[row, col] + mu * sNxNy[row, col];
                    k[1, 1] = (lambda + mu) * sNyNy[row, col] +
                        mu * (sNxNx[row, col] + sNyNy[row, col] + sNzNz[row, col]);
                    k[1, 2] = lambda * sNyNz[row, col] + mu * sNzNy[row, col];
                    k[2, 0] = lambda * sNzNx[row, col] + mu * sNxNz[row, col];
                    k[2, 1] = lambda * sNzNy[row, col] + mu * sNyNz[row, col];
                    k[2, 2] = (lambda + mu) * sNzNz[row, col] +
                        mu * (sNxNx[row, col] + sNyNy[row, col] + sNzNz[row, col]);

                    m[0, 0] = rho * sNN[row, col];
                    m[0, 1] = 0.0;
                    m[0, 2] = 0.0;
                    m[1, 0] = 0.0;
                    m[1, 1] = rho * sNN[row, col];
                    m[1, 2] = 0.0;
                    m[2, 0] = 0.0;
                    m[2, 1] = 0.0;
                    m[2, 2] = rho * sNN[row, col];

                    for (int rowDof = 0; rowDof < dof; rowDof++)
                    {
                        for (int colDof = 0; colDof < dof; colDof++)
                        {
                            A[rowNodeId * dof + rowDof, colNodeId * dof + colDof] +=
                                (1.0 / (beta * dt * dt)) * m[rowDof, colDof] +
                                k[rowDof, colDof];

                            B[rowNodeId * dof + rowDof] +=
                                m[rowDof, colDof] * (
                                (1.0 / (beta * dt * dt)) * u[colDof] +
                                (1.0 / (beta * dt)) * vel[colDof] +
                                (1.0 / (2.0 * beta) - 1.0) * acc[colDof]);
                        }
                    }
                }
            }

            for (int row = 0; row < elemNodeCnt; row++)
            {
                int rowNodeId = nodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                for (int rowDof = 0; rowDof < dof; rowDof++)
                {
                    B[rowNodeId * dof + rowDof] += rho * g[rowDof] * sN[row];
                }
            }
        }
    }
}
