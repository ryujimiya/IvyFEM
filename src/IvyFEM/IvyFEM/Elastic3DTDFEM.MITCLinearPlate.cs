﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Animation;

namespace IvyFEM
{
    partial class Elastic3DTDFEM
    {
        protected IvyFEM.Lapack.DoubleMatrix CalcMITCLinearPlateKe(
            TriangleFE d1TriFE, double h, OpenTK.Vector3d[] xPts,
            IvyFEM.Lapack.DoubleMatrix Cb, IvyFEM.Lapack.DoubleMatrix Cs, double kappa)
        {
            return Elastic3DFEMUtils.CalcMITCLinearPlateKe(d1TriFE, h, xPts, Cb, Cs, kappa);
        }

        protected IvyFEM.Lapack.DoubleMatrix CalcMITCLinearPlateMe(
            TriangleFE d1TriFE, double h, OpenTK.Vector3d[] xPts, double rho)
        {
            return Elastic3DFEMUtils.CalcMITCLinearPlateMe(d1TriFE, h, xPts, rho);
        }

        protected void CalcMITCLinearPlateElementAB(
            uint feId, IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            {
                uint quantityId0 = 0;
                TriangleFE workTriFE = World.GetTriangleFE(quantityId0, feId);
                uint workMaId = workTriFE.MaterialId;
                if (!World.IsMaterialId(workMaId))
                {
                    return;
                }
                Material workMa0 = World.GetMaterial(workMaId);
                if (!(workMa0 is MITCLinearPlateMaterial))
                {
                    return;
                }
            }

            uint dQuantityId = 0; // displacement
            uint rQuantityId = 1; // rotation
            System.Diagnostics.Debug.Assert(World.GetDof(dQuantityId) == 3);
            System.Diagnostics.Debug.Assert(World.GetDof(rQuantityId) == 3);
            int dDof = 3; // u v w
            int rDof = 3; //θx θy θz
            int dNodeCnt = (int)World.GetNodeCount(dQuantityId);
            int rNodeCnt = (int)World.GetNodeCount(rQuantityId);
            int rOffset = dNodeCnt * dDof;

            double dt = TimeStep;
            double beta = NewmarkBeta;
            double gamma = NewmarkGamma;
            System.Diagnostics.Debug.Assert(UValueIds.Count == 2);
            var dFV = World.GetFieldValue(UValueIds[0]);
            var rFV = World.GetFieldValue(UValueIds[1]);

            TriangleFE dTriFE = World.GetTriangleFE(dQuantityId, feId);
            TriangleFE rTriFE = World.GetTriangleFE(rQuantityId, feId);
            uint maId = dTriFE.MaterialId;
            if (!World.IsMaterialId(maId))
            {
                return;
            }
            Material ma0 = World.GetMaterial(maId);
            System.Diagnostics.Debug.Assert(ma0 is MITCLinearPlateMaterial);

            int[] dCoIds = dTriFE.NodeCoordIds;
            uint dElemNodeCnt = dTriFE.NodeCount;
            // 暫定で要素マトリクスをこのメソッド内に直接記述
            System.Diagnostics.Debug.Assert(dElemNodeCnt == 3); // 1次要素
            int[] dNodes = new int[dElemNodeCnt];
            for (int iNode = 0; iNode < dElemNodeCnt; iNode++)
            {
                int coId = dCoIds[iNode];
                int nodeId = World.Coord2Node(dQuantityId, coId);
                dNodes[iNode] = nodeId;
            }
            int[] rCoIds = rTriFE.NodeCoordIds;
            uint rElemNodeCnt = rTriFE.NodeCount;
            // 暫定で要素マトリクスをこのメソッド内に直接記述
            System.Diagnostics.Debug.Assert(rElemNodeCnt == 3); // 1次要素
            int[] rNodes = new int[rElemNodeCnt];
            for (int iNode = 0; iNode < rElemNodeCnt; iNode++)
            {
                int coId = rCoIds[iNode];
                int nodeId = World.Coord2Node(rQuantityId, coId);
                rNodes[iNode] = nodeId;
            }
            OpenTK.Vector3d[] dPts = new OpenTK.Vector3d[dElemNodeCnt];
            for (int iNode = 0; iNode < dElemNodeCnt; iNode++)
            {
                int coId = dCoIds[iNode];
                double[] coord = World.GetCoord(dQuantityId, coId);
                dPts[iNode] = new OpenTK.Vector3d(coord[0], coord[1], coord[2]);
            }

            var ma = ma0 as MITCLinearPlateMaterial;
            double h = ma.Thickness;
            double rho = ma.MassDensity;
            double E = ma.Young;
            double nu = ma.Poisson;
            double kappa = ma.ShearCorrectionFactor;
            IvyFEM.Lapack.DoubleMatrix Cb = new IvyFEM.Lapack.DoubleMatrix(3, 3);
            {
                double k1 = E / (1.0 - nu * nu);
                Cb[0, 0] = k1 * 1.0;
                Cb[0, 1] = k1 * nu;
                Cb[0, 2] = k1 * 0.0;
                Cb[1, 1] = k1 * 1.0;
                Cb[1, 2] = k1 * 0.0;
                Cb[2, 2] = k1 * (1.0 / 2.0) * (1.0 - nu);
                for (int i = 0; i < 3; i++)
                {
                    for (int j = i + 1; j < 3; j++)
                    {
                        Cb[j, i] = Cb[i, j];
                    }
                }
            }
            IvyFEM.Lapack.DoubleMatrix Cs = new IvyFEM.Lapack.DoubleMatrix(2, 2);
            {
                double k1 = E / (2.0 * (1.0 + nu));
                Cs[0, 0] = k1 * 1.0;
                Cs[0, 1] = k1 * 0.0;
                Cs[1, 0] = k1 * 0.0;
                Cs[1, 1] = k1 * 1.0;
            }

            //------------------------------------------
            // local
            var Ke = CalcMITCLinearPlateKe(dTriFE, h, dPts, Cb, Cs, kappa);
            var Me = CalcMITCLinearPlateMe(dTriFE, h, dPts, rho);
            // Note: グローバル座標へ変換済み

            // local dof
            int localDof = 6;
            System.Diagnostics.Debug.Assert(localDof == (dDof + rDof));
            int localROffset = dDof;
            // displacement
            for (int row = 0; row < dElemNodeCnt; row++)
            {
                int rowNodeId = dNodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                // displacement
                for (int col = 0; col < dElemNodeCnt; col++)
                {
                    int colNodeId = dNodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    int colCoId = dCoIds[col];
                    double[] u = dFV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                    double[] vel = dFV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                    double[] acc = dFV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);

                    for (int rowDof = 0; rowDof < dDof; rowDof++)
                    {
                        for (int colDof = 0; colDof < dDof; colDof++)
                        {
                            double kValue = Ke[row * localDof + rowDof, col * localDof + colDof];
                            double mValue = Me[row * localDof + rowDof, col * localDof + colDof];
                            A[rowNodeId * dDof + rowDof, colNodeId * dDof + colDof] +=
                                (1.0 / (beta * dt * dt)) * mValue +
                                kValue;

                            B[rowNodeId * dDof + rowDof] +=
                                mValue * (
                                (1.0 / (beta * dt * dt)) * u[colDof] +
                                (1.0 / (beta * dt)) * vel[colDof] +
                                (1.0 / (2.0 * beta) - 1.0) * acc[colDof]);
                        }
                    }
                }
                // rotation
                for (int col = 0; col < rElemNodeCnt; col++)
                {
                    int colNodeId = rNodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    int colCoId = rCoIds[col];
                    double[] u = rFV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                    double[] vel = rFV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                    double[] acc = rFV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);

                    for (int rowDof = 0; rowDof < dDof; rowDof++)
                    {
                        for (int colDof = 0; colDof < rDof; colDof++)
                        {
                            double kValue = Ke[row * localDof + rowDof, col * localDof + colDof + localROffset];
                            double mValue = Me[row * localDof + rowDof, col * localDof + colDof + localROffset];
                            A[rowNodeId * dDof + rowDof, colNodeId * rDof + colDof + rOffset] +=
                                (1.0 / (beta * dt * dt)) * mValue +
                                kValue;

                            B[rowNodeId * dDof + rowDof] +=
                                mValue * (
                                (1.0 / (beta * dt * dt)) * u[colDof] +
                                (1.0 / (beta * dt)) * vel[colDof] +
                                (1.0 / (2.0 * beta) - 1.0) * acc[colDof]);
                        }
                    }
                }
            }
            // rotation
            for (int row = 0; row < rElemNodeCnt; row++)
            {
                int rowNodeId = rNodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                // displacement
                for (int col = 0; col < dElemNodeCnt; col++)
                {
                    int colNodeId = dNodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    int colCoId = dCoIds[col];
                    double[] u = dFV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                    double[] vel = dFV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                    double[] acc = dFV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);

                    for (int rowDof = 0; rowDof < rDof; rowDof++)
                    {
                        for (int colDof = 0; colDof < dDof; colDof++)
                        {
                            double kValue = Ke[row * localDof + rowDof + localROffset, col * localDof + colDof];
                            double mValue = Me[row * localDof + rowDof + localROffset, col * localDof + colDof];
                            A[rowNodeId * rDof + rowDof + rOffset, colNodeId * dDof + colDof] +=
                                (1.0 / (beta * dt * dt)) * mValue +
                                kValue;

                            B[rowNodeId * rDof + rowDof + rOffset] +=
                                mValue * (
                                (1.0 / (beta * dt * dt)) * u[colDof] +
                                (1.0 / (beta * dt)) * vel[colDof] +
                                (1.0 / (2.0 * beta) - 1.0) * acc[colDof]);
                        }
                    }
                }
                // rotation
                for (int col = 0; col < rElemNodeCnt; col++)
                {
                    int colNodeId = rNodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    int colCoId = rCoIds[col];
                    double[] u = rFV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                    double[] vel = rFV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                    double[] acc = rFV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);

                    for (int rowDof = 0; rowDof < rDof; rowDof++)
                    {
                        for (int colDof = 0; colDof < rDof; colDof++)
                        {
                            double kValue =
                                Ke[row * localDof + rowDof + localROffset, col * localDof + colDof + localROffset];
                            double mValue =
                                Me[row * localDof + rowDof + localROffset, col * localDof + colDof + localROffset];
                            A[rowNodeId * rDof + rowDof + rOffset, colNodeId * rDof + colDof + rOffset] +=
                                (1.0 / (beta * dt * dt)) * mValue +
                                kValue;

                            B[rowNodeId * rDof + rowDof + rOffset] +=
                                mValue * (
                                (1.0 / (beta * dt * dt)) * u[colDof] +
                                (1.0 / (beta * dt)) * vel[colDof] +
                                (1.0 / (2.0 * beta) - 1.0) * acc[colDof]);
                        }
                    }
                }
            }
        }
    }
}
