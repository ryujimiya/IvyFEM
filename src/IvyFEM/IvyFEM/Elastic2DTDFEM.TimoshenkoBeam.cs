using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DTDFEM
    {
        protected IvyFEM.Lapack.DoubleMatrix CalcTimoshenkoBeamKl(
            uint dElemNodeCnt, uint rElemNodeCnt, int dDof, int rDof,
            LineFE dLineFE, LineFE rLineFE,
            double E, double G, double kappa, double Ae, double Iz)
        {
            return Elastic2DFEMUtils.CalcTimoshenkoBeamKl(
                dElemNodeCnt, rElemNodeCnt, dDof, rDof,
                dLineFE, rLineFE,
                E, G, kappa, Ae, Iz);
        }

        protected IvyFEM.Lapack.DoubleMatrix CalcTimoshenkoBeamMl(
            uint dElemNodeCnt, uint rElemNodeCnt, int dDof, int rDof,
            LineFE dLineFE, LineFE rLineFE,
            double rho, double Ae, double Iz)
        {
            return Elastic2DFEMUtils.CalcTimoshenkoBeamMl(
                dElemNodeCnt, rElemNodeCnt, dDof, rDof,
                dLineFE, rLineFE,
                rho, Ae, Iz);
        }

        protected void CalcTimoshenkoBeamElementABForLine(
            uint feId, IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            {
                // Note: feIdは線要素
                uint quantityId0 = 0;
                LineFE workLineFE = World.GetLineFE(quantityId0, feId);
                uint workMaId = workLineFE.MaterialId;
                if (!World.IsMaterialId(workMaId))
                {
                    return;
                }
                Material workMa0 = World.GetMaterial(workMaId);
                if (!(workMa0 is TimoshenkoBeamMaterial))
                {
                    return;
                }
            }

            uint dQuantityId = 0; // displacement
            uint rQuantityId = 1; // rotation
            System.Diagnostics.Debug.Assert(World.GetDof(dQuantityId) == 1);
            System.Diagnostics.Debug.Assert(World.GetDof(rQuantityId) == 1);
            int dDof = 1; // w
            int rDof = 1; //θ
            int dNodeCnt = (int)World.GetNodeCount(dQuantityId);
            int rNodeCnt = (int)World.GetNodeCount(rQuantityId);
            int offset = dNodeCnt * dDof;

            double dt = TimeStep;
            double beta = NewmarkBeta;
            double gamma = NewmarkGamma;
            System.Diagnostics.Debug.Assert(UValueIds.Count == 2);
            var dFV = World.GetFieldValue(UValueIds[0]);
            var rFV = World.GetFieldValue(UValueIds[1]);

            // Note: feIdは線要素
            LineFE dLineFE = World.GetLineFE(dQuantityId, feId);
            LineFE rLineFE = World.GetLineFE(rQuantityId, feId);
            uint maId = dLineFE.MaterialId;
            if (!World.IsMaterialId(maId))
            {
                return;
            }
            Material ma0 = World.GetMaterial(maId);
            System.Diagnostics.Debug.Assert(ma0 is TimoshenkoBeamMaterial);

            int[] dCoIds = dLineFE.NodeCoordIds;
            uint dElemNodeCnt = dLineFE.NodeCount;
            int[] dNodes = new int[dElemNodeCnt];
            for (int iNode = 0; iNode < dElemNodeCnt; iNode++)
            {
                int coId = dCoIds[iNode];
                int nodeId = World.Coord2Node(dQuantityId, coId);
                dNodes[iNode] = nodeId;
            }
            int[] rCoIds = rLineFE.NodeCoordIds;
            uint rElemNodeCnt = rLineFE.NodeCount;
            int[] rNodes = new int[rElemNodeCnt];
            for (int iNode = 0; iNode < rElemNodeCnt; iNode++)
            {
                int coId = rCoIds[iNode];
                int nodeId = World.Coord2Node(rQuantityId, coId);
                rNodes[iNode] = nodeId;
            }
            int[] vCoIds = dLineFE.VertexCoordIds;
            uint elemVertexCnt = dLineFE.VertexCount;
            double[][] vCoords = new double[elemVertexCnt][];
            for (int iVertex = 0; iVertex < elemVertexCnt; iVertex++)
            {
                int coId = vCoIds[iVertex];
                double[] coord = World.GetCoord(dQuantityId, coId);
                vCoords[iVertex] = coord;
            }

            var ma = ma0 as TimoshenkoBeamMaterial;
            double Ae = ma.Area;
            double Iz = ma.SecondMomentOfArea;
            double rho = ma.MassDensity;
            double E = ma.Young;
            double G = ma.ShearCoefficient;
            double kappa = ma.ShearCorrectionFactor;

            double[] pt1 = vCoords[0];
            double[] pt2 = vCoords[1];
            // y座標は0固定
            System.Diagnostics.Debug.Assert(Math.Abs(pt1[1]) < IvyFEM.Constants.PrecisionLowerLimit);
            System.Diagnostics.Debug.Assert(Math.Abs(pt2[1]) < IvyFEM.Constants.PrecisionLowerLimit);
            double le = Math.Sqrt(
                (pt2[0] - pt1[0]) * (pt2[0] - pt1[0]) +
                (pt2[1] - pt1[1]) * (pt2[1] - pt1[1]));

            var Ke = CalcTimoshenkoBeamKl(
                dElemNodeCnt, rElemNodeCnt, dDof, rDof, dLineFE, rLineFE,
                E, G, kappa, Ae, Iz);
            var Me = CalcTimoshenkoBeamMl(
                dElemNodeCnt, rElemNodeCnt, dDof, rDof, dLineFE, rLineFE,
                rho, Ae, Iz);

            // local dof
            int localDof = 2;
            System.Diagnostics.Debug.Assert(localDof == (dDof + rDof));
            int localOffset = dDof;

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
                    int colDof = 0;

                    double kValue = Ke[row * localDof, col * localDof];
                    double mValue = Me[row * localDof, col * localDof];
                    A[rowNodeId, colNodeId] +=
                        (1.0 / (beta * dt * dt)) * mValue +
                        kValue;

                    B[rowNodeId] +=
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        (1.0 / (2.0 * beta) - 1.0) * acc[colDof]);
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
                    int colDof = 0;

                    double kValue = Ke[row * localDof, col * localDof + localOffset];
                    double mValue = Me[row * localDof, col * localDof + localOffset];
                    A[rowNodeId, colNodeId + offset] +=
                        (1.0 / (beta * dt * dt)) * mValue +
                        kValue;

                    B[rowNodeId] +=
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        (1.0 / (2.0 * beta) - 1.0) * acc[colDof]);
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
                    int colDof = 0;

                    double kValue = Ke[row * localDof + localOffset, col * localDof];
                    double mValue = Me[row * localDof + localOffset, col * localDof];
                    A[rowNodeId + offset, colNodeId] +=
                        (1.0 / (beta * dt * dt)) * mValue +
                        kValue;

                    B[rowNodeId + offset] +=
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        (1.0 / (2.0 * beta) - 1.0) * acc[colDof]);
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
                    int colDof = 0;

                    double kValue = Ke[row * localDof + localOffset, col * localDof + localOffset];
                    double mValue = Me[row * localDof + localOffset, col * localDof + localOffset];
                    A[rowNodeId + offset, colNodeId + offset] +=
                        (1.0 / (beta * dt * dt)) * mValue +
                        kValue;

                    B[rowNodeId + offset] +=
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        (1.0 / (2.0 * beta) - 1.0) * acc[colDof]);
                }
            }
        }
    }
}
