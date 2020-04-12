using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DTDFEM
    {
        protected void CalcTimoshenkoFrameElementABForLine(
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
                if (!(workMa0 is TimoshenkoFrameMaterial))
                {
                    return;
                }
            }

            System.Diagnostics.Debug.Assert(DisplacementQuantityIds.Count == 2);
            System.Diagnostics.Debug.Assert(DisplacementQuantityIds[0] == 0);
            System.Diagnostics.Debug.Assert(DisplacementQuantityIds[1] == 1);
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
            System.Diagnostics.Debug.Assert(ma0 is TimoshenkoFrameMaterial);

            int[] d1CoIds = d1LineFE.NodeCoordIds;
            uint d1ElemNodeCnt = d1LineFE.NodeCount;
            // 暫定で要素マトリクスをこのメソッド内に直接記述
            System.Diagnostics.Debug.Assert(d1ElemNodeCnt == 2); // 1次要素
            int[] d1Nodes = new int[d1ElemNodeCnt];
            for (int iNode = 0; iNode < d1ElemNodeCnt; iNode++)
            {
                int coId = d1CoIds[iNode];
                int nodeId = World.Coord2Node(d1QuantityId, coId);
                d1Nodes[iNode] = nodeId;
            }
            int[] d2CoIds = d2LineFE.NodeCoordIds;
            uint d2ElemNodeCnt = d2LineFE.NodeCount;
            int[] d2Nodes = new int[d2ElemNodeCnt];
            for (int iNode = 0; iNode < d2ElemNodeCnt; iNode++)
            {
                int coId = d2CoIds[iNode];
                int nodeId = World.Coord2Node(d2QuantityId, coId);
                d2Nodes[iNode] = nodeId;
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
            int[] vCoIds = d1LineFE.VertexCoordIds;
            uint elemVertexCnt = d1LineFE.VertexCount;
            double[][] vCoords = new double[elemVertexCnt][];
            for (int iVertex = 0; iVertex < d1ElemNodeCnt; iVertex++)
            {
                int coId = vCoIds[iVertex];
                double[] coord = World.GetCoord(d1QuantityId, coId);
                vCoords[iVertex] = coord;
            }

            var ma = ma0 as TimoshenkoFrameMaterial;
            double Ae = ma.Area;
            double Iz = ma.SecondMomentOfArea;
            double Ix = ma.PolarSecondMomentOfArea;
            double rho = ma.MassDensity;
            double E = ma.Young;
            double G = ma.ShearCoefficient;
            double kappa = ma.TimoshenkoShearCoefficient;

            // local dof
            double[] pt1 = vCoords[0];
            double[] pt2 = vCoords[1];
            double le = Math.Sqrt(
                (pt2[0] - pt1[0]) * (pt2[0] - pt1[0]) +
                (pt2[1] - pt1[1]) * (pt2[1] - pt1[1]));
            double cosXX = (pt2[0] - pt1[0]) / le;
            double cosXY = (pt2[1] - pt1[1]) / le;
            double cosYX = -1.0 * (pt2[1] - pt1[1]) / le;
            double cosYY = (pt2[0] - pt1[0]) / le;

            // local dof
            int localDof = 3;
            System.Diagnostics.Debug.Assert(localDof == (d1Dof + d2Dof + rDof));
            int localD2Offset = d1Dof;
            int localROffset = localD2Offset + d2Dof;
            // 次数の高い要素に合わせる
            int localMaxNodeCnt = (int)Math.Max(Math.Max(d1ElemNodeCnt, d2ElemNodeCnt), rElemNodeCnt);
            var T = new IvyFEM.Lapack.DoubleMatrix(localMaxNodeCnt * localDof, localMaxNodeCnt * localDof);
            // d1
            for (int row = 0; row < d1ElemNodeCnt; row++)
            {
                // d1
                {
                    int col = row;
                    T[row * localDof, col * localDof] = cosXX;
                }
                // d2
                {
                    int col = row;
                    T[row * localDof, col * localDof + localD2Offset] = cosXY;
                }
                // r
                {
                    int col = row;
                    T[row * localDof, col * localDof + localROffset] = 0.0;
                }
            }
            // d2
            for (int row = 0; row < d2ElemNodeCnt; row++)
            {
                // d1
                {
                    int col = row;
                    T[row * localDof + localD2Offset, col * localDof] = cosYX;
                }
                // d2
                {
                    int col = row;
                    T[row * localDof + localD2Offset, col * localDof + localD2Offset] = cosYY;
                }
                // r
                {
                    int col = row;
                    T[row * localDof + localD2Offset, col * localDof + localROffset] = 0.0;
                }
            }
            // r
            for (int row = 0; row < rElemNodeCnt; row++)
            {
                // d1
                {
                    int col = row;
                    T[row * localDof + localROffset, col * localDof] = 0.0;
                }
                // d2
                {
                    int col = row;
                    T[row * localDof + localROffset, col * localDof + localD2Offset] = 0.0;
                }
                // r
                {
                    int col = row;
                    T[row * localDof + localROffset, col * localDof + localROffset] = 1.0;
                }
            }
            var transT = IvyFEM.Lapack.DoubleMatrix.Transpose(T);

            // Truss
            var localKeTruss = CalcTrussLocalKe(le, E, Ae);
            var localMeTruss = CalcTrussLocalMe(le, rho, Ae);

            // Timoshenko Beam
            int localBeamMaxNodeCnt =(int) Math.Max(d2ElemNodeCnt, rElemNodeCnt);
            var localKeBeam = new IvyFEM.Lapack.DoubleMatrix(
                localBeamMaxNodeCnt * (d2Dof + rDof), localBeamMaxNodeCnt * (d2Dof + rDof));
            var localMeBeam = new IvyFEM.Lapack.DoubleMatrix(
                localBeamMaxNodeCnt * (d2Dof + rDof), localBeamMaxNodeCnt * (d2Dof + rDof));
            int beamDof = 2;
            int localBeamROffset = 1;
            System.Diagnostics.Debug.Assert(beamDof == (d2Dof + rDof));
            System.Diagnostics.Debug.Assert(localBeamROffset == d2Dof);
            IntegrationPoints ipK;
            if (d2LineFE.Order == 1 && rLineFE.Order == 1)
            {
                // 低減積分
                ipK = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point1);
                System.Diagnostics.Debug.Assert(ipK.Ls.Length == 1);
            }
            else
            {
                ipK = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);
                System.Diagnostics.Debug.Assert(ipK.Ls.Length == 5);
            }
            for (int ipPt = 0; ipPt < ipK.PointCount; ipPt++)
            {
                double[] L = ipK.Ls[ipPt];
                double[] d2N = d2LineFE.CalcN(L);
                double[][] d2Nu = d2LineFE.CalcNu(L);
                double[] d2Nx = d2Nu[0];
                double[] rN = rLineFE.CalcN(L);
                double[][] rNu = rLineFE.CalcNu(L);
                double[] rNx = rNu[0];
                double lineLen = d2LineFE.GetLineLength();
                double weight = ipK.Weights[ipPt];
                double detJWeight = (lineLen / 2.0) * weight;

                // displacement
                for (int row = 0; row < d2ElemNodeCnt; row++)
                {
                    // displacement
                    for (int col = 0; col < d2ElemNodeCnt; col++)
                    {
                        double kValue = detJWeight * kappa * G * Ae * d2Nx[row] * d2Nx[col];
                        localKeBeam[row * beamDof, col * beamDof] += kValue;
                    }
                    // rotation
                    for (int col = 0; col < rElemNodeCnt; col++)
                    {
                        double kValue = -1.0 * detJWeight * kappa * G * Ae * d2Nx[row] * rN[col];
                        localKeBeam[row * beamDof, col * beamDof + localBeamROffset] += kValue;
                    }
                }
                // rotation
                for (int row = 0; row < rElemNodeCnt; row++)
                {
                    // displacement
                    for (int col = 0; col < d2ElemNodeCnt; col++)
                    {
                        double kValue = -1.0 * detJWeight * kappa * G * Ae * rN[row] * d2Nx[col];
                        localKeBeam[row * beamDof + localBeamROffset, col * beamDof] += kValue;
                    }
                    // rotation
                    for (int col = 0; col < rElemNodeCnt; col++)
                    {
                        double kValue1 = detJWeight * kappa * G * Ae * rN[row] * rN[col];
                        double kValue2 = detJWeight * E * Iz * rNx[row] * rNx[col];
                        localKeBeam[row * beamDof  + localBeamROffset, col * beamDof + localBeamROffset] +=
                            kValue1 + kValue2;
                    }
                }
            }
            IntegrationPoints ipM = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);
            System.Diagnostics.Debug.Assert(ipM.Ls.Length == 5);
            for (int ipPt = 0; ipPt < ipM.PointCount; ipPt++)
            {
                double[] L = ipM.Ls[ipPt];
                double[] d2N = d2LineFE.CalcN(L);
                double[][] d2Nu = d2LineFE.CalcNu(L);
                double[] d2Nx = d2Nu[0];
                double[] rN = rLineFE.CalcN(L);
                double[][] rNu = rLineFE.CalcNu(L);
                double[] rNx = rNu[0];
                double lineLen = d2LineFE.GetLineLength();
                double weight = ipM.Weights[ipPt];
                double detJWeight = (lineLen / 2.0) * weight;

                // displacement
                for (int row = 0; row < d2ElemNodeCnt; row++)
                {
                    // displacement
                    for (int col = 0; col < d2ElemNodeCnt; col++)
                    {
                        double mValue = detJWeight * rho * Ae * d2N[row] * d2N[col];
                        localMeBeam[row * beamDof, col * beamDof] += mValue;
                    }
                    // rotation
                    for (int col = 0; col < rElemNodeCnt; col++)
                    {
                        double mValue = 0.0;
                        localMeBeam[row * beamDof, col * beamDof + localBeamROffset] += mValue;
                    }
                }
                // rotation
                for (int row = 0; row < rElemNodeCnt; row++)
                {
                    // displacement
                    for (int col = 0; col < d2ElemNodeCnt; col++)
                    {
                        double mValue = 0.0;
                        localMeBeam[row * beamDof + localBeamROffset, col * beamDof] += mValue;
                    }
                    // rotation
                    for (int col = 0; col < rElemNodeCnt; col++)
                    {
                        double mValue = detJWeight * rho * Ix * rN[row] * rN[col];
                        localMeBeam[row * beamDof + localBeamROffset, col * beamDof + localBeamROffset] +=
                            mValue;
                    }
                }
            }

            var localKe = new IvyFEM.Lapack.DoubleMatrix(localMaxNodeCnt * localDof, localMaxNodeCnt * localDof);
            var localMe = new IvyFEM.Lapack.DoubleMatrix(localMaxNodeCnt * localDof, localMaxNodeCnt * localDof);
            // d1
            for (int row = 0; row < d1ElemNodeCnt; row++)
            {
                // d1
                for (int col = 0; col < d1ElemNodeCnt; col++)
                {
                    localKe[row * localDof, col * localDof] = localKeTruss[row * d1Dof, col * d1Dof];
                    localMe[row * localDof, col * localDof] = localMeTruss[row * d1Dof, col * d1Dof];
                }
                // d2
                for (int col = 0; col < d2ElemNodeCnt; col++)
                {
                    localKe[row * localDof, col * localDof + localD2Offset] = 0.0;
                    localMe[row * localDof, col * localDof + localD2Offset] = 0.0;
                }
                // r
                for (int col = 0; col < rElemNodeCnt; col++)
                {
                    localKe[row * localDof, col * localDof + localROffset] = 0.0;
                    localMe[row * localDof, col * localDof + localROffset] = 0.0;
                }
            }
            // d2
            for (int row = 0; row < d2ElemNodeCnt; row++)
            {
                // d1
                for (int col = 0; col < d1ElemNodeCnt; col++)
                {
                    localKe[row * localDof + localD2Offset, col * localDof] = 0.0;
                    localMe[row * localDof + localD2Offset, col * localDof] = 0.0;
                }
                // d2
                for (int col = 0; col < d2ElemNodeCnt; col++)
                {
                    localKe[row * localDof + localD2Offset, col * localDof + localD2Offset] =
                        localKeBeam[row * (d2Dof + rDof), col * (d2Dof + rDof)];
                    localMe[row * localDof + localD2Offset, col * localDof + localD2Offset] =
                        localMeBeam[row * (d2Dof + rDof), col * (d2Dof + rDof)];
                }
                // r
                for (int col = 0; col < rElemNodeCnt; col++)
                {
                    localKe[row * localDof + localD2Offset, col * localDof + localROffset] =
                        localKeBeam[row * (d2Dof + rDof), col * (d2Dof + rDof) + localBeamROffset];
                    localMe[row * localDof + localD2Offset, col * localDof + localROffset] =
                        localMeBeam[row * (d2Dof + rDof), col * (d2Dof + rDof) + localBeamROffset];
                }
            }
            // r
            for (int row = 0; row < rElemNodeCnt; row++)
            {
                // d1
                for (int col = 0; col < d1ElemNodeCnt; col++)
                {
                    localKe[row * localDof + localROffset, col * localDof] = 0.0;
                    localMe[row * localDof + localROffset, col * localDof] = 0.0;
                }
                // d2
                for (int col = 0; col < d2ElemNodeCnt; col++)
                {
                    localKe[row * localDof + localROffset, col * localDof + localD2Offset] =
                        localKeBeam[row * (d2Dof + rDof) + localBeamROffset, col * (d2Dof + rDof)];
                    localMe[row * localDof + localROffset, col * localDof + localD2Offset] =
                        localMeBeam[row * (d2Dof + rDof) + localBeamROffset, col * (d2Dof + rDof)];
                }
                // r
                for (int col = 0; col < rElemNodeCnt; col++)
                {
                    localKe[row * localDof + localROffset, col * localDof + localROffset] =
                        localKeBeam[row * (d2Dof + rDof) + localBeamROffset, col * (d2Dof + rDof) + localBeamROffset];
                    localMe[row * localDof + localROffset, col * localDof + localROffset] =
                        localMeBeam[row * (d2Dof + rDof) + localBeamROffset, col * (d2Dof + rDof) + localBeamROffset];
                }
            }

            var Ke = (transT * localKe) * T;
            var Me = (transT * localMe) * T;

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
                    A[rowNodeId, colNodeId] +=
                        (1.0 / (beta * dt * dt)) * mValue +
                        kValue;

                    B[rowNodeId] +=
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        (1.0 / (2.0 * beta) - 1.0) * acc[colDof]);
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
                    A[rowNodeId, colNodeId + d2Offset] +=
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

                    double kValue = Ke[row * localDof, col * localDof + localROffset];
                    double mValue = Me[row * localDof, col * localDof + localROffset];
                    A[rowNodeId, colNodeId + rOffset] +=
                        (1.0 / (beta * dt * dt)) * mValue +
                        kValue;

                    B[rowNodeId] +=
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        (1.0 / (2.0 * beta) - 1.0) * acc[colDof]);
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
                    A[rowNodeId + d2Offset, colNodeId] +=
                        (1.0 / (beta * dt * dt)) * mValue +
                        kValue;

                    B[rowNodeId + d2Offset] +=
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        (1.0 / (2.0 * beta) - 1.0) * acc[colDof]);
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
                    A[rowNodeId + d2Offset, colNodeId + d2Offset] +=
                        (1.0 / (beta * dt * dt)) * mValue +
                        kValue;

                    B[rowNodeId + d2Offset] +=
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

                    double kValue = Ke[row * localDof + localD2Offset, col * localDof + localROffset];
                    double mValue = Me[row * localDof + localD2Offset, col * localDof + localROffset];
                    A[rowNodeId + d2Offset, colNodeId + rOffset] +=
                        (1.0 / (beta * dt * dt)) * mValue +
                        kValue;

                    B[rowNodeId + d2Offset] +=
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
                    A[rowNodeId + rOffset, colNodeId] +=
                        (1.0 / (beta * dt * dt)) * mValue +
                        kValue;

                    B[rowNodeId + rOffset] +=
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        (1.0 / (2.0 * beta) - 1.0) * acc[colDof]);
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
                    A[rowNodeId + rOffset, colNodeId + d2Offset] +=
                        (1.0 / (beta * dt * dt)) * mValue +
                        kValue;

                    B[rowNodeId + rOffset] +=
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

                    double kValue = Ke[row * localDof + localROffset, col * localDof + localROffset];
                    double mValue = Me[row * localDof + localROffset, col * localDof + localROffset];
                    A[rowNodeId + rOffset, colNodeId + rOffset] +=
                        (1.0 / (beta * dt * dt)) * mValue +
                        kValue;

                    B[rowNodeId + rOffset] +=
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        (1.0 / (2.0 * beta) - 1.0) * acc[colDof]);
                }
            }
        }
    }
}
