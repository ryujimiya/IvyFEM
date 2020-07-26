using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DTDFEM
    {
        protected void CalcFieldConsistentTLFrameKl(
            double l0, double E, double Ae, double Iz, double[] ul,
            out double[] fl, out IvyFEM.Lapack.DoubleMatrix kl)
        {
            Elastic2DFEMUtils.CalcFieldConsistentTLFrameKl(l0, E, Ae, Iz, ul, out fl, out kl);
        }

        protected void CalcFieldConsistentTLFrameMl(
            double l0, double rho, double Ae, double Iz,
            double[] ul, double[] velUl, double[] accUl,
            out double[] fkl,
            out IvyFEM.Lapack.DoubleMatrix ml,
            out IvyFEM.Lapack.DoubleMatrix ckl, 
            out IvyFEM.Lapack.DoubleMatrix kkl)
        {
            Elastic2DFEMUtils.CalcFieldConsistentTLFrameMl(l0, rho, Ae, Iz, ul, velUl, accUl,
                out fkl, out ml, out ckl, out kkl);
        }

        protected void CalcFieldConsistentTLFrameElementABForLine(
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
                if (!(workMa0 is FieldConsistentTLFrameMaterial))
                {
                    return;
                }
            }

            uint d1QuantityId = 0; // displacement
            uint d2QuantityId = 1; // displacement
            uint rQuantityId = 2; // rotation
            System.Diagnostics.Debug.Assert(World.GetDof(d1QuantityId) == 1);
            System.Diagnostics.Debug.Assert(World.GetDof(d2QuantityId) == 1);
            System.Diagnostics.Debug.Assert(World.GetDof(rQuantityId) == 1);
            int d1Dof = 1; // u
            int d2Dof = 1; // v
            int rDof = 1; //θ
            int d1NodeCnt = (int)World.GetNodeCount(d1QuantityId);
            int d2NodeCnt = (int)World.GetNodeCount(d2QuantityId);
            int rNodeCnt = (int)World.GetNodeCount(rQuantityId);
            int d2Offset = d1NodeCnt * d1Dof;
            int rOffset = d2Offset + d2NodeCnt * d2Dof;

            double dt = TimeStep;
            double beta = NewmarkBeta;
            double gamma = NewmarkGamma;
            System.Diagnostics.Debug.Assert(UValueIds.Count == 3);
            var d1FV = World.GetFieldValue(UValueIds[0]);
            var d2FV = World.GetFieldValue(UValueIds[1]);
            var rFV = World.GetFieldValue(UValueIds[2]);

            // Note: feIdは線要素
            LineFE d1LineFE = World.GetLineFE(d1QuantityId, feId);
            LineFE d2LineFE = World.GetLineFE(d2QuantityId, feId);
            LineFE rLineFE = World.GetLineFE(rQuantityId, feId);
            uint maId = d1LineFE.MaterialId;
            if (!World.IsMaterialId(maId))
            {
                return;
            }
            Material ma0 = World.GetMaterial(maId);
            System.Diagnostics.Debug.Assert(ma0 is FieldConsistentTLFrameMaterial);

            int[] d1CoIds = d1LineFE.NodeCoordIds;
            uint d1ElemNodeCnt = d1LineFE.NodeCount;
            // 暫定で要素マトリクスをこのメソッド内に直接記述
            System.Diagnostics.Debug.Assert(d1ElemNodeCnt == 3); // 2次要素
            int[] d1Nodes = new int[d1ElemNodeCnt];
            for (int iNode = 0; iNode < d1ElemNodeCnt; iNode++)
            {
                int coId = d1CoIds[iNode];
                int nodeId = World.Coord2Node(d1QuantityId, coId);
                d1Nodes[iNode] = nodeId;
            }
            int[] d2CoIds = d2LineFE.NodeCoordIds;
            uint d2ElemNodeCnt = d2LineFE.NodeCount;
            // 暫定で要素マトリクスをこのメソッド内に直接記述
            System.Diagnostics.Debug.Assert(d2ElemNodeCnt == 2); // 1次要素
            int[] d2Nodes = new int[d2ElemNodeCnt];
            for (int iNode = 0; iNode < d2ElemNodeCnt; iNode++)
            {
                int coId = d2CoIds[iNode];
                int nodeId = World.Coord2Node(d2QuantityId, coId);
                d2Nodes[iNode] = nodeId;
            }
            int[] rCoIds = rLineFE.NodeCoordIds;
            uint rElemNodeCnt = rLineFE.NodeCount;
            // 暫定で要素マトリクスをこのメソッド内に直接記述
            System.Diagnostics.Debug.Assert(rElemNodeCnt == 2); // 1次要素
            int[] rNodes = new int[rElemNodeCnt];
            for (int iNode = 0; iNode < rElemNodeCnt; iNode++)
            {
                int coId = rCoIds[iNode];
                int nodeId = World.Coord2Node(rQuantityId, coId);
                rNodes[iNode] = nodeId;
            }
            int[] vCoIds = d1LineFE.VertexCoordIds;
            uint elemVertexCnt = d1LineFE.VertexCount;
            double[][] vCoords = new double[elemVertexCnt][];
            for (int iVertex = 0; iVertex < elemVertexCnt; iVertex++)
            {
                int coId = vCoIds[iVertex];
                double[] coord = World.GetCoord(d1QuantityId, coId);
                vCoords[iVertex] = coord;
            }

            System.Diagnostics.Debug.Assert(d1ElemNodeCnt == 3);
            System.Diagnostics.Debug.Assert(d2ElemNodeCnt == 2);
            System.Diagnostics.Debug.Assert(rElemNodeCnt == d2ElemNodeCnt);
            System.Diagnostics.Debug.Assert(d1ElemNodeCnt + d2ElemNodeCnt + rElemNodeCnt == 7);
            double[] ug = new double[7];
            for (int iNode = 0; iNode < d1ElemNodeCnt; iNode++)
            {
                int d1NodeId = d1Nodes[iNode];
                if (d1NodeId == -1)
                {
                    continue;
                }
                ug[iNode * 3] = U[d1NodeId];
            }
            for (int iNode = 0; iNode < d2ElemNodeCnt; iNode++)
            {
                int d2NodeId = d2Nodes[iNode];
                if (d2NodeId == -1)
                {
                    continue;
                }
                ug[iNode * 3 + 1] = U[d2NodeId + d2Offset];
            }
            for (int iNode = 0; iNode < rElemNodeCnt; iNode++)
            {
                int rNodeId = rNodes[iNode];
                if (rNodeId == -1)
                {
                    continue;
                }
                ug[iNode * 3 + 2] = U[rNodeId + rOffset];
            }

            //---------------------------------------------
            // vel, acc
            double prevU1 = d1FV.GetDoubleValue(d1CoIds[0], FieldDerivativeType.Value)[0];
            double prevVelU1 = d1FV.GetDoubleValue(d1CoIds[0], FieldDerivativeType.Velocity)[0];
            double prevAccU1 = d1FV.GetDoubleValue(d1CoIds[0], FieldDerivativeType.Acceleration)[0];
            double prevU2 = d1FV.GetDoubleValue(d1CoIds[1], FieldDerivativeType.Value)[0];
            double prevVelU2 = d1FV.GetDoubleValue(d1CoIds[1], FieldDerivativeType.Velocity)[0];
            double prevAccU2 = d1FV.GetDoubleValue(d1CoIds[1], FieldDerivativeType.Acceleration)[0];
            double prevV1 = d2FV.GetDoubleValue(d2CoIds[0], FieldDerivativeType.Value)[0];
            double prevVelV1 = d2FV.GetDoubleValue(d2CoIds[0], FieldDerivativeType.Velocity)[0];
            double prevAccV1 = d2FV.GetDoubleValue(d2CoIds[0], FieldDerivativeType.Acceleration)[0];
            double prevV2 = d2FV.GetDoubleValue(d2CoIds[1], FieldDerivativeType.Value)[0];
            double prevVelV2 = d2FV.GetDoubleValue(d2CoIds[1], FieldDerivativeType.Velocity)[0];
            double prevAccV2 = d2FV.GetDoubleValue(d2CoIds[1], FieldDerivativeType.Acceleration)[0];
            double prevT1 = rFV.GetDoubleValue(rCoIds[0], FieldDerivativeType.Value)[0];
            double prevVelT1 = rFV.GetDoubleValue(rCoIds[0], FieldDerivativeType.Velocity)[0];
            double prevAccT1 = rFV.GetDoubleValue(rCoIds[0], FieldDerivativeType.Acceleration)[0];
            double prevT2 = rFV.GetDoubleValue(rCoIds[1], FieldDerivativeType.Value)[0];
            double prevVelT2 = rFV.GetDoubleValue(rCoIds[1], FieldDerivativeType.Velocity)[0];
            double prevAccT2 = rFV.GetDoubleValue(rCoIds[1], FieldDerivativeType.Acceleration)[0];
            double prevU3 = d1FV.GetDoubleValue(d1CoIds[2], FieldDerivativeType.Value)[0];
            double prevVelU3 = d1FV.GetDoubleValue(d1CoIds[2], FieldDerivativeType.Velocity)[0];
            double prevAccU3 = d1FV.GetDoubleValue(d1CoIds[2], FieldDerivativeType.Acceleration)[0];
            double[] prevUg = { prevU1, prevV1, prevT1, prevU2, prevV2, prevT2, prevU3 };
            double[] prevVelUg = { prevVelU1, prevVelV1, prevVelT1, prevVelU2, prevVelV2, prevVelT2, prevVelU3 };
            double[] prevAccUg = { prevAccU1, prevAccV1, prevAccT1, prevAccU2, prevAccV2, prevAccT2, prevAccU3 };
            var velUg = new double[7];
            var accUg = new double[7];
            for (int i = 0; i < 7; i++)
            {
                velUg[i] =
                    (gamma / (beta * dt)) * (ug[i] - prevUg[i]) +
                    (1.0 - gamma / beta) * prevVelUg[i] +
                    dt * (1.0 - gamma / (2.0 * beta)) * prevAccUg[i];
                accUg[i] =
                    (1.0 / (beta * dt * dt)) * (ug[i] - prevUg[i]) -
                    (1.0 / (beta * dt)) * prevVelUg[i] -
                    (1.0 / (2.0 * beta) - 1.0) * prevAccUg[i];
            }

            var ma = ma0 as FieldConsistentTLFrameMaterial;
            double Ae = ma.Area;
            double Iz = ma.SecondMomentOfArea;
            double rho = ma.MassDensity;
            double E = ma.Young;

            double[] pt10 = vCoords[0];
            double[] pt20 = vCoords[1];
            double l0 = Math.Sqrt(
                (pt20[0] - pt10[0]) * (pt20[0] - pt10[0]) +
                (pt20[1] - pt10[1]) * (pt20[1] - pt10[1]));
            double c0 = (pt20[0] - pt10[0]) / l0;
            double s0 = (pt20[1] - pt10[1]) / l0;

            var T = new IvyFEM.Lapack.DoubleMatrix(7, 7);
            T[0, 0] = c0;
            T[0, 1] = s0;
            T[0, 2] = 0.0;
            T[0, 3] = 0.0;
            T[0, 4] = 0.0;
            T[0, 5] = 0.0;
            T[0, 6] = 0.0;
            T[1, 0] = -s0;
            T[1, 1] = c0;
            T[1, 2] = 0.0;
            T[1, 3] = 0.0;
            T[1, 4] = 0.0;
            T[1, 5] = 0.0;
            T[1, 6] = 0.0;
            T[2, 0] = 0.0;
            T[2, 1] = 0.0;
            T[2, 2] = 1.0;
            T[2, 3] = 0.0;
            T[2, 4] = 0.0;
            T[2, 5] = 0.0;
            T[2, 6] = 0.0;
            T[3, 0] = 0.0;
            T[3, 1] = 0.0;
            T[3, 2] = 0.0;
            T[3, 3] = c0;
            T[3, 4] = s0;
            T[3, 5] = 0.0;
            T[3, 6] = 0.0;
            T[4, 0] = 0.0;
            T[4, 1] = 0.0;
            T[4, 2] = 0.0;
            T[4, 3] = -s0;
            T[4, 4] = c0;
            T[4, 5] = 0.0;
            T[4, 6] = 0.0;
            T[5, 0] = 0.0;
            T[5, 1] = 0.0;
            T[5, 2] = 0.0;
            T[5, 3] = 0.0;
            T[5, 4] = 0.0;
            T[5, 5] = 1.0;
            T[5, 6] = 0.0;
            T[6, 0] = 0.0;
            T[6, 1] = 0.0;
            T[6, 2] = 0.0;
            T[6, 3] = 0.0;
            T[6, 4] = 0.0;
            T[6, 5] = 0.0;
            T[6, 6] = 1.0;
            var transT = IvyFEM.Lapack.DoubleMatrix.Transpose(T);

            double[] ul = T * ug;
            double[] velUl = T * velUg;
            double[] accUl = T * accUg;

            //---------------------------------
            // local
            double[] fl;
            IvyFEM.Lapack.DoubleMatrix kl;
            
            CalcFieldConsistentTLFrameKl(l0, E, Ae, Iz, ul, out fl, out kl);
            //---------------------------------

            //---------------------------------
            // f & K
            double[] f = transT * fl;
            var Ke = (transT * kl) * T;
            //---------------------------------

            //-----------------------------
            // local
            double[] fkl;
            IvyFEM.Lapack.DoubleMatrix ml;
            IvyFEM.Lapack.DoubleMatrix ckl;
            IvyFEM.Lapack.DoubleMatrix kkl;
            
            CalcFieldConsistentTLFrameMl(l0, rho, Ae, Iz, ul, velUl, accUl,
                out fkl, out ml, out ckl, out kkl);
            //---------------------------------

            //---------------------------------
            // fk, M, etc
            double[] fk = transT * fkl;
            var Me = (transT * ml) * T;
            var Cke = (transT * ckl) * T;
            var Kke = (transT * kkl) * T;
            //---------------------------------

            // local dof
            int localDof = 3;
            System.Diagnostics.Debug.Assert(localDof == (d1Dof + d2Dof + rDof));
            int localD2Offset = d1Dof;
            int localROffset = localD2Offset + d2Dof;
            // displacement 1
            for (int row = 0; row < d1ElemNodeCnt; row++)
            {
                int rowNodeId = d1Nodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                // displacement 1
                for (int col = 0; col < d1ElemNodeCnt; col++)
                {
                    int colNodeId = d1Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    int colCoId = d1CoIds[col];
                    double[] u = d1FV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                    double[] vel = d1FV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                    double[] acc = d1FV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);
                    int colDof = 0;

                    double kValue = Ke[row * localDof, col * localDof];
                    double mValue = Me[row * localDof, col * localDof];
                    double ckValue = Cke[row * localDof, col * localDof];
                    double kkValue = Kke[row * localDof, col * localDof];
                    double curU = ug[col * localDof];
                    double curVel = velUg[col * localDof];
                    double curAcc = accUg[col * localDof];
                    A[rowNodeId, colNodeId] +=
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue +
                        (gamma / (beta * dt)) * ckValue +
                        kkValue;
                    B[rowNodeId] +=
                        kValue * curU +
                        mValue * curAcc +
                        ckValue * curVel +
                        kkValue * curU +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]) +
                        ckValue * (
                        (gamma / (beta * dt)) * u[colDof] +
                        ((gamma / beta) - 1.0) * vel[colDof] +
                        ((gamma / (2.0 * beta)) - 1.0) * acc[colDof]);
                }
                // displacement 2
                for (int col = 0; col < d2ElemNodeCnt; col++)
                {
                    int colNodeId = d2Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    int colCoId = d2CoIds[col];
                    double[] u = d2FV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                    double[] vel = d2FV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                    double[] acc = d2FV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);
                    int colDof = 0;

                    double kValue = Ke[row * localDof, col * localDof + localD2Offset];
                    double mValue = Me[row * localDof, col * localDof + localD2Offset];
                    double ckValue = Cke[row * localDof, col * localDof + localD2Offset];
                    double kkValue = Kke[row * localDof, col * localDof + localD2Offset];
                    double curU = ug[col * localDof + localD2Offset];
                    double curVel = velUg[col * localDof + localD2Offset];
                    double curAcc = accUg[col * localDof + localD2Offset];
                    A[rowNodeId, colNodeId + d2Offset] +=
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue +
                        (gamma / (beta * dt)) * ckValue +
                        kkValue;
                    B[rowNodeId] +=
                        kValue * curU +
                        mValue * curAcc +
                        ckValue * curVel +
                        kkValue * curU +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]) +
                        ckValue * (
                        (gamma / (beta * dt)) * u[colDof] +
                        ((gamma / beta) - 1.0) * vel[colDof] +
                        ((gamma / (2.0 * beta)) - 1.0) * acc[colDof]);
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

                    double kValue = Ke[row * localDof, col * localDof + localROffset];
                    double mValue = Me[row * localDof, col * localDof + localROffset];
                    double ckValue = Cke[row * localDof, col * localDof + localROffset];
                    double kkValue = Kke[row * localDof, col * localDof + localROffset];
                    double curU = ug[col * localDof + localROffset];
                    double curVel = velUg[col * localDof + localROffset];
                    double curAcc = accUg[col * localDof + localROffset];
                    A[rowNodeId, colNodeId + rOffset] +=
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue +
                        (gamma / (beta * dt)) * ckValue +
                        kkValue;
                    B[rowNodeId] +=
                        kValue * curU +
                        mValue * curAcc +
                        ckValue * curVel +
                        kkValue * curU +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]) +
                        ckValue * (
                        (gamma / (beta * dt)) * u[colDof] +
                        ((gamma / beta) - 1.0) * vel[colDof] +
                        ((gamma / (2.0 * beta)) - 1.0) * acc[colDof]);
                }
            }
            // displacement 2
            for (int row = 0; row < d2ElemNodeCnt; row++)
            {
                int rowNodeId = d2Nodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                // displacement 1
                for (int col = 0; col < d1ElemNodeCnt; col++)
                {
                    int colNodeId = d1Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    int colCoId = d1CoIds[col];
                    double[] u = d1FV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                    double[] vel = d1FV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                    double[] acc = d1FV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);
                    int colDof = 0;

                    double kValue = Ke[row * localDof + localD2Offset, col * localDof];
                    double mValue = Me[row * localDof + localD2Offset, col * localDof];
                    double ckValue = Cke[row * localDof + localD2Offset, col * localDof];
                    double kkValue = Kke[row * localDof + localD2Offset, col * localDof];
                    double curU = ug[col * localDof];
                    double curVel = velUg[col * localDof];
                    double curAcc = accUg[col * localDof];
                    A[rowNodeId + d2Offset, colNodeId] +=
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue +
                        (gamma / (beta * dt)) * ckValue +
                        kkValue;
                    B[rowNodeId + d2Offset] +=
                        kValue * curU +
                        mValue * curAcc +
                        ckValue * curVel +
                        kkValue * curU +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]) +
                        ckValue * (
                        (gamma / (beta * dt)) * u[colDof] +
                        ((gamma / beta) - 1.0) * vel[colDof] +
                        ((gamma / (2.0 * beta)) - 1.0) * acc[colDof]);
                }
                // displacement 2
                for (int col = 0; col < d2ElemNodeCnt; col++)
                {
                    int colNodeId = d2Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    int colCoId = d2CoIds[col];
                    double[] u = d2FV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                    double[] vel = d2FV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                    double[] acc = d2FV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);
                    int colDof = 0;

                    double kValue = Ke[row * localDof + localD2Offset, col * localDof + localD2Offset];
                    double mValue = Me[row * localDof + localD2Offset, col * localDof + localD2Offset];
                    double ckValue = Cke[row * localDof + localD2Offset, col * localDof + localD2Offset];
                    double kkValue = Kke[row * localDof + localD2Offset, col * localDof + localD2Offset];
                    double curU = ug[col * localDof + localD2Offset];
                    double curVel = velUg[col * localDof + localD2Offset];
                    double curAcc = accUg[col * localDof + localD2Offset];
                    A[rowNodeId + d2Offset, colNodeId + d2Offset] +=
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue +
                        (gamma / (beta * dt)) * ckValue +
                        kkValue;
                    B[rowNodeId + d2Offset] +=
                        kValue * curU +
                        mValue * curAcc +
                        ckValue * curVel +
                        kkValue * curU +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]) +
                        ckValue * (
                        (gamma / (beta * dt)) * u[colDof] +
                        ((gamma / beta) - 1.0) * vel[colDof] +
                        ((gamma / (2.0 * beta)) - 1.0) * acc[colDof]);
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

                    double kValue = Ke[row * localDof + localD2Offset, col * localDof + localROffset];
                    double mValue = Me[row * localDof + localD2Offset, col * localDof + localROffset];
                    double ckValue = Cke[row * localDof + localD2Offset, col * localDof + localROffset];
                    double kkValue = Kke[row * localDof + localD2Offset, col * localDof + localROffset];
                    double curU = ug[col * localDof + localROffset];
                    double curVel = velUg[col * localDof + localROffset];
                    double curAcc = accUg[col * localDof + localROffset];
                    A[rowNodeId + d2Offset, colNodeId + rOffset] +=
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue +
                        (gamma / (beta * dt)) * ckValue +
                        kkValue;
                    B[rowNodeId + d2Offset] +=
                        kValue * curU +
                        mValue * curAcc +
                        ckValue * curVel +
                        kkValue * curU +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]) +
                        ckValue * (
                        (gamma / (beta * dt)) * u[colDof] +
                        ((gamma / beta) - 1.0) * vel[colDof] +
                        ((gamma / (2.0 * beta)) - 1.0) * acc[colDof]);
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
                // displacement 1
                for (int col = 0; col < d1ElemNodeCnt; col++)
                {
                    int colNodeId = d1Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    int colCoId = d1CoIds[col];
                    double[] u = d1FV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                    double[] vel = d1FV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                    double[] acc = d1FV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);
                    int colDof = 0;

                    double kValue = Ke[row * localDof + localROffset, col * localDof];
                    double mValue = Me[row * localDof + localROffset, col * localDof];
                    double ckValue = Cke[row * localDof + localROffset, col * localDof];
                    double kkValue = Kke[row * localDof + localROffset, col * localDof];
                    double curU = ug[col * localDof];
                    double curVel = velUg[col * localDof];
                    double curAcc = accUg[col * localDof];
                    A[rowNodeId + rOffset, colNodeId] +=
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue +
                        (gamma / (beta * dt)) * ckValue +
                        kkValue;
                    B[rowNodeId + rOffset] +=
                        kValue * curU +
                        mValue * curAcc +
                        ckValue * curVel +
                        kkValue * curU +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]) +
                        ckValue * (
                        (gamma / (beta * dt)) * u[colDof] +
                        ((gamma / beta) - 1.0) * vel[colDof] +
                        ((gamma / (2.0 * beta)) - 1.0) * acc[colDof]);
                }
                // displacement 2
                for (int col = 0; col < d2ElemNodeCnt; col++)
                {
                    int colNodeId = d2Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    int colCoId = d2CoIds[col];
                    double[] u = d2FV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                    double[] vel = d2FV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                    double[] acc = d2FV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);
                    int colDof = 0;

                    double kValue = Ke[row * localDof + localROffset, col * localDof + localD2Offset];
                    double mValue = Me[row * localDof + localROffset, col * localDof + localD2Offset];
                    double ckValue = Cke[row * localDof + localROffset, col * localDof + localD2Offset];
                    double kkValue = Kke[row * localDof + localROffset, col * localDof + localD2Offset];
                    double curU = ug[col * localDof + localD2Offset];
                    double curVel = velUg[col * localDof + localD2Offset];
                    double curAcc = accUg[col * localDof + localD2Offset];
                    A[rowNodeId + rOffset, colNodeId + d2Offset] +=
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue +
                        (gamma / (beta * dt)) * ckValue +
                        kkValue;
                    B[rowNodeId + rOffset] +=
                        kValue * curU +
                        mValue * curAcc +
                        ckValue * curVel +
                        kkValue * curU +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]) +
                        ckValue * (
                        (gamma / (beta * dt)) * u[colDof] +
                        ((gamma / beta) - 1.0) * vel[colDof] +
                        ((gamma / (2.0 * beta)) - 1.0) * acc[colDof]);
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

                    double kValue = Ke[row * localDof + localROffset, col * localDof + localROffset];
                    double mValue = Me[row * localDof + localROffset, col * localDof + localROffset];
                    double ckValue = Cke[row * localDof + localROffset, col * localDof + localROffset];
                    double kkValue = Kke[row * localDof + localROffset, col * localDof + localROffset];
                    double curU = ug[col * localDof + localROffset];
                    double curVel = velUg[col * localDof + localROffset];
                    double curAcc = accUg[col * localDof + localROffset];
                    A[rowNodeId + rOffset, colNodeId + rOffset] +=
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue +
                        (gamma / (beta * dt)) * ckValue +
                        kkValue;
                    B[rowNodeId + rOffset] +=
                        kValue * curU +
                        mValue * curAcc +
                        ckValue * curVel +
                        kkValue * curU +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]) +
                        ckValue * (
                        (gamma / (beta * dt)) * u[colDof] +
                        ((gamma / beta) - 1.0) * vel[colDof] +
                        ((gamma / (2.0 * beta)) - 1.0) * acc[colDof]);
                }
            }
            // internal force
            for (int row = 0; row < d1ElemNodeCnt; row++)
            {
                int rowNodeId = d1Nodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                B[rowNodeId] +=
                    -f[row * localDof] - fk[row * localDof];
            }
            for (int row = 0; row < d2ElemNodeCnt; row++)
            {
                int rowNodeId = d2Nodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                B[rowNodeId + d2Offset] +=
                    -f[row * localDof + localD2Offset] - fk[row * localDof + localD2Offset];
            }
            for (int row = 0; row < rElemNodeCnt; row++)
            {
                int rowNodeId = rNodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                B[rowNodeId + rOffset] +=
                    -f[row * localDof + localROffset] - fk[row * localDof + localROffset];
            }
        }
    }
}
