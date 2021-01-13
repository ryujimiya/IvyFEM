﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class AdvectionDiffusion3DTDFEM : FEM
    {
        public double TimeStep { get; private set; } = 0;
        public double NewmarkBeta { get; private set; } = 1.0 / 4.0;
        public double NewmarkGamma { get; private set; } = 1.0 / 2.0;
        uint ValueId { get; set; } = 0;
        uint PrevValueId { get; set; } = 0;

        public uint VeloValueId { get; set; } = 0;

        // output
        public double[] U { get; private set; } = null;

        public AdvectionDiffusion3DTDFEM(FEWorld world, uint veloValueId,
            double timeStep,
            double newmarkBeta, double newmarkGamma,
            uint valueId, uint prevValueId)
        {
            World = world;
            VeloValueId = veloValueId;
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
            var veloFV = World.GetFieldValue(VeloValueId);
            IList<uint> feIds = World.GetTetrahedronFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = World.GetTetrahedronFE(quantityId, feId);
                uint elemNodeCnt = tetFE.NodeCount;
                int[] coIds = tetFE.NodeCoordIds;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = coIds[iNode];
                    int nodeId = World.Coord2Node(quantityId, coId);
                    nodes[iNode] = nodeId;
                }

                Material ma0 = World.GetMaterial(tetFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is DiffusionMaterial);
                var ma = ma0 as DiffusionMaterial;
                double nu = ma.DiffusionCoef;
                double f = ma.F;
                double[][] velos = new double[elemNodeCnt][];
                for (int iNode=0; iNode < elemNodeCnt; iNode++)
                {
                    int nodeId = nodes[iNode];
                    if (nodeId == -1)
                    {
                        continue;
                    }
                    int coId = World.Node2Coord(quantityId, nodeId);
                    double[] velo = veloFV.GetDoubleValue(coId, FieldDerivativeType.Value);
                    velos[iNode] = velo;
                }

                double tau;
                {
                    double[] aveVelo = {
                        (velos[0][0] + velos[1][0] + velos[2][0] + velos[3][0]) / 4.0,
                        (velos[0][1] + velos[1][1] + velos[2][1] + velos[3][1]) / 4.0,
                        (velos[0][2] + velos[1][2] + velos[2][2] + velos[3][2]) / 4.0
                    };
                    double veloNorm =
                        Math.Sqrt(aveVelo[0] * aveVelo[0] + aveVelo[1] * aveVelo[1] + aveVelo[2] * aveVelo[2]);
                    double[] veloDir = { aveVelo[0] / veloNorm, aveVelo[1] / veloNorm, aveVelo[2] / veloNorm };
                    double h;
                    {
                        double[] L = new double[] { 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0 };
                        double[][] Nu = tetFE.CalcNu(L);
                        double[] Nx = Nu[0];
                        double[] Ny = Nu[1];
                        double[] Nz = Nu[2];
                        double tmp = 0;
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            tmp += Math.Abs(veloDir[0] * Nx[iNode] + veloDir[1] * Ny[iNode] + veloDir[2] * Nz[iNode]);
                        }
                        h = 2.0 / tmp;
                    }
                    double lambda = (1.0 / 2.0) * veloNorm * h / nu;

                    //tau = (1.0 / 2.0) * (h / veloNorm) * (1.0 / Math.Tanh(lambda) - 1.0 / lambda);
                    if (nu > 1.0e-20)
                    {
                        if (lambda < 3.0)
                        {
                            tau = (1.0 / 2.0) * h  / veloNorm * lambda / 3.0;
                        }
                        else
                        {
                            tau = (1.0 / 2.0) * h / veloNorm;
                        }
                    }
                    else
                    {
                        tau = (1.0 / 2.0) * h / veloNorm;
                    }
                }

                double[] sN = tetFE.CalcSN();
                IntegrationPoints ip = TetrahedronFE.GetIntegrationPoints(World.TetIntegrationPointCount);//Point5
                for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                {
                    double[] L = ip.Ls[ipPt];
                    double[] N = tetFE.CalcN(L);
                    double[][] Nu = tetFE.CalcNu(L);
                    double[] Nx = Nu[0];
                    double[] Ny = Nu[1];
                    double[] Nz = Nu[2];
                    double detJ = tetFE.GetDetJacobian(L);
                    double weight = ip.Weights[ipPt];
                    double detJWeight = (1.0 / 6.0) * weight * detJ;

                    double[] velo = new double[3];
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        for (int iDof = 0; iDof < 3; iDof++)
                        {
                            velo[iDof] += N[iNode] * velos[iNode][iDof];
                        }
                    }

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

                            double k = detJWeight * nu * (Nx[row] * Nx[col] + Ny[row] * Ny[col] + Nz[row] * Nz[col]);
                            double c = detJWeight * 
                                (N[row] + tau * (velo[0] * Nx[row] + velo[1] * Ny[row] + velo[2] * Nz[row])) *
                                (velo[0] * Nx[col] + velo[1] * Ny[col] + velo[2] * Nz[col]);
                            double m = detJWeight * N[row] * N[col];

                            A[rowNodeId, colNodeId] += k + c + (gamma / (beta * dt)) * m;
                            B[rowNodeId] +=
                                m * (
                                (gamma / (beta * dt)) * u[0] -
                                (1.0 - gamma / beta) * vel[0] -
                                dt * (1.0 - gamma / (2.0 * beta)) * acc[0]
                                );
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