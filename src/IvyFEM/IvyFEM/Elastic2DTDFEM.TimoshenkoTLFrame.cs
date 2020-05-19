using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DTDFEM
    {
        protected void CalcTimoshenkoTLFrameElementABForLine(
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
                if (!(workMa0 is TimoshenkoTLFrameMaterial))
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
            System.Diagnostics.Debug.Assert(ma0 is TimoshenkoTLFrameMaterial);

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
            for (int iVertex = 0; iVertex < elemVertexCnt; iVertex++)
            {
                int coId = vCoIds[iVertex];
                double[] coord = World.GetCoord(d1QuantityId, coId);
                vCoords[iVertex] = coord;
            }

            var ma = ma0 as TimoshenkoTLFrameMaterial;
            double Ae = ma.Area;
            double Iz = ma.SecondMomentOfArea;
            double rho = ma.MassDensity;
            double E = ma.Young;
            double G = ma.ShearCoefficient;
            double kappa = ma.ShearCorrectionFactor;

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
            System.Diagnostics.Debug.Assert(d1ElemNodeCnt == d2ElemNodeCnt);
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

            double[] ug = new double[localMaxNodeCnt * localDof];
            // d1
            for (int row = 0; row < d1ElemNodeCnt; row++)
            {
                int rowNodeId = d1Nodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                ug[row * localDof] = U[rowNodeId];
            }
            // d2
            for (int row = 0; row < d2ElemNodeCnt; row++)
            {
                int rowNodeId = d2Nodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                ug[row * localDof + localD2Offset] = U[rowNodeId + d2Offset];
            }
            // r
            for (int row = 0; row < rElemNodeCnt; row++)
            {
                int rowNodeId = rNodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                ug[row * localDof + localROffset] = U[rowNodeId + rOffset];
            }

            double[] ul = T * ug;

            //--------------------------
            // local
            double[] fl = new double[localMaxNodeCnt * localDof];
            var localKe = new IvyFEM.Lapack.DoubleMatrix(
                localMaxNodeCnt * localDof, localMaxNodeCnt * localDof);
            var localMe = new IvyFEM.Lapack.DoubleMatrix(
                localMaxNodeCnt * localDof, localMaxNodeCnt * localDof);

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
                double[] d1N = d1LineFE.CalcN(L);
                double[][] d1Nu = d1LineFE.CalcNu(L);
                double[] d1Np = d1Nu[0];
                double[] d2N = d2LineFE.CalcN(L);
                double[][] d2Nu = d2LineFE.CalcNu(L);
                double[] d2Np = d2Nu[0];
                double[] rN = rLineFE.CalcN(L);
                double[][] rNu = rLineFE.CalcNu(L);
                double[] rNp = rNu[0];
                double lineLen = d2LineFE.GetLineLength();
                double weight = ipK.Weights[ipPt];
                double detJWeight = (lineLen / 2.0) * weight;

                double[] Nx = new double[localMaxNodeCnt * localDof];
                double[] Ny = new double[localMaxNodeCnt * localDof];
                double[] Nt = new double[localMaxNodeCnt * localDof];
                double[] Npx = new double[localMaxNodeCnt * localDof];
                double[] Npy = new double[localMaxNodeCnt * localDof];
                double[] Npt = new double[localMaxNodeCnt * localDof];
                for (int i = 0; i < d1ElemNodeCnt; i++)
                {
                    Nx[i * localDof] = d1N[i];
                    Npx[i * localDof] = d1Np[i];
                }
                for (int i = 0; i < d2ElemNodeCnt; i++)
                {
                    Ny[i * localDof + localD2Offset] = d2N[i];
                    Npy[i * localDof + localD2Offset] = d2Np[i];
                }
                for (int i = 0; i < rElemNodeCnt; i++)
                {
                    Nt[i * localDof + localROffset] = rN[i];
                    Npt[i * localDof + localROffset] = rNp[i];
                }
                double[] upxu = Npx;
                double[] upyu = Npy;
                double[] tu = Nt;
                double[] tpu = Npt;

                double upx = 0.0;
                for (int i = 0; i < d1ElemNodeCnt; i++)
                {
                    upx += d1Np[i] * ul[i * localDof];
                }
                double upy = 0.0;
                for (int i = 0; i < d2ElemNodeCnt; i++)
                {
                    upy += d2Np[i] * ul[i * localDof + localD2Offset];
                }
                double t = 0.0;
                double tp = 0.0;
                for (int i = 0; i < rElemNodeCnt; i++)
                {
                    t += rN[i] * ul[i * localDof + localROffset];
                    tp += rNp[i] * ul[i * localDof + localROffset];
                }
                double cost = Math.Cos(t);
                double sint = Math.Sin(t);

                double c1 = (1.0 + upx) * cost * cost + upy * sint * cost - cost;
                double c2 = (1.0 + upx) * cost * sint + upy * sint * sint - sint;
                double c3 = -(1.0 + upx) * (1.0 + upx) * cost * sint + (1.0 + upx) * upy * cost * cost;
                double c4 = -(1.0 + upx) * upy * sint * sint + upy * upy * sint * cost;
                double c5 = (1.0 + upx) * sint - upy * cost;
                double c6 = tp;
                double c7 = (1.0 + upx) * sint * sint - upy * cost * sint;
                double c8 = -(1.0 + upx) * sint * cost + upy * cost * cost;
                double c9 = (1.0 + upx) * (1.0 + upx) * sint * cost + (1.0 + upx) * upy * sint * sint;
                double c10 = -(1.0 + upx) * upy * cost * cost - upy * upy * cost * sint;

                for (int row = 0; row < localMaxNodeCnt * localDof; row++)
                {
                    fl[row] += detJWeight * E * Ae * c1 * Npx[row] +
                        detJWeight * E * Ae * c2 * Npy[row] +
                        detJWeight * E * Ae * c3 * Nt[row] +
                        detJWeight * E * Ae * c4 * Nt[row] +
                        detJWeight * E * Ae * c5 * Nt[row] +
                        detJWeight * E * Iz * c6 * Npt[row] +
                        detJWeight * kappa * G * Ae * c7 * Npx[row] +
                        detJWeight * kappa * G * Ae * c8 * Npy[row] +
                        detJWeight * kappa * G * Ae * c9 * Nt[row] +
                        detJWeight * kappa * G * Ae * c10 * Nt[row];
                }
                for (int row = 0; row < localMaxNodeCnt * localDof; row++)
                {
                    for (int col = 0; col < localMaxNodeCnt * localDof; col++)
                    {
                        double c1u = upxu[col] * cost * cost -
                            (1.0 + upx) * 2.0 * cost * sint * tu[col] +
                            upyu[col] * sint * cost +
                            upy * (cost * cost - sint * sint) * tu[col] +
                            sint * tu[col];
                        double c2u = upxu[col] * cost * sint +
                            (1.0 + upx) * (cost * cost - sint * sint) * tu[col] +
                            upyu[col] * sint * sint +
                            upy * 2.0 * sint * cost * tu[col] -
                            cost * tu[col];
                        double c3u = -2.0 * (1.0 + upx) * upxu[col] * cost * sint -
                            (1.0 + upx) * (1.0 + upx) * (cost * cost - sint * sint) * tu[col] +
                            upxu[col] * upy * cost * cost +
                            (1.0 + upx) * upyu[col] * cost * cost -
                            (1.0 + upx) * upy * 2.0 * cost * sint * tu[col];
                        double c4u = -upxu[col] * upy * sint * sint -
                            (1.0 + upx) * upyu[col] * sint * sint -
                            (1.0 + upx) * upy * 2.0 * sint * cost * tu[col] +
                            2.0 * upy * upyu[col] * sint * cost +
                            upy * upy * (cost * cost - sint * sint) * tu[col];
                        double c5u = upxu[col] * sint +
                            (1.0 + upx) * cost * tu[col] -
                            upyu[col] * cost +
                            upy * sint * tu[col];
                        double c6u = tpu[col];
                        double c7u = upxu[col] * sint * sint +
                            (1.0 + upx) * 2.0 * sint * cost * tu[col] -
                            upyu[col] * cost * sint -
                            upy * (cost * cost - sint * sint) * tu[col];
                        double c8u = -upxu[col] * sint * cost -
                            (1.0 + upx) * (cost * cost - sint * sint) * tu[col] +
                            upyu[col] * cost * cost -
                            upy * 2.0 * cost * sint * tu[col];
                        double c9u = 2.0 * (1.0 + upx) * upxu[col] * sint * cost +
                            (1.0 + upx) * (1.0 + upx) * (cost * cost - sint * sint) * tu[col] +
                            upxu[col] * upy * sint * sint +
                            (1.0 + upx) * upyu[col] * sint * sint +
                            (1.0 + upx) * upy * 2.0 * sint * cost * tu[col];
                        double c10u = -upxu[col] * upy * cost * cost -
                            (1.0 + upx) * upyu[col] * cost * cost +
                            (1.0 + upx) * upy * 2.0 * cost * sint * tu[col] -
                            2.0 * upy * upyu[col] * cost * sint -
                            upy * upy * (cost * cost - sint * sint) * tu[col];

                        localKe[row, col] += detJWeight * E * Ae * c1u * Npx[row] +
                            detJWeight * E * Ae * c2u * Npy[row] +
                            detJWeight * E * Ae * c3u * Nt[row] +
                            detJWeight * E * Ae * c4u * Nt[row] +
                            detJWeight * E * Ae * c5u * Nt[row] +
                            detJWeight * E * Iz * c6u * Npt[row] +
                            detJWeight * kappa * G * Ae * c7u * Npx[row] +
                            detJWeight * kappa * G * Ae * c8u * Npy[row] +
                            detJWeight * kappa * G * Ae * c9u * Nt[row] +
                            detJWeight * kappa * G * Ae * c10u * Nt[row];
                    }
                }
            }
            IntegrationPoints ipM;
            ipM = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);
            System.Diagnostics.Debug.Assert(ipM.Ls.Length == 5);
            for (int ipPt = 0; ipPt < ipK.PointCount; ipPt++)
            {
                double[] L = ipK.Ls[ipPt];
                double[] d1N = d1LineFE.CalcN(L);
                double[][] d1Nu = d1LineFE.CalcNu(L);
                double[] d1Np = d1Nu[0];
                double[] d2N = d2LineFE.CalcN(L);
                double[][] d2Nu = d2LineFE.CalcNu(L);
                double[] d2Np = d2Nu[0];
                double[] rN = rLineFE.CalcN(L);
                double[][] rNu = rLineFE.CalcNu(L);
                double[] rNp = rNu[0];
                double lineLen = d2LineFE.GetLineLength();
                double weight = ipK.Weights[ipPt];
                double detJWeight = (lineLen / 2.0) * weight;

                double[] Nx = new double[localMaxNodeCnt * localDof];
                double[] Ny = new double[localMaxNodeCnt * localDof];
                double[] Nt = new double[localMaxNodeCnt * localDof];
                for (int i = 0; i < d1ElemNodeCnt; i++)
                {
                    Nx[i * localDof] = d1N[i];
                }
                for (int i = 0; i < d2ElemNodeCnt; i++)
                {
                    Ny[i * localDof + localD2Offset] = d2N[i];
                }
                for (int i = 0; i < rElemNodeCnt; i++)
                {
                    Nt[i * localDof + localROffset] = rN[i];
                }
                for (int row = 0; row < localMaxNodeCnt * localDof; row++)
                {
                    for (int col = 0; col < localMaxNodeCnt * localDof; col++)
                    {
                        double mex = rho * Ae * Nx[row] * Nx[col];
                        double mey = rho * Ae * Ny[row] * Ny[col];
                        double met = rho * Iz * Nt[row] * Nt[col];
                        localMe[row, col] += mex + mey + met;
                    }
                }
            }
            //--------------------------

            double[] f = transT * fl;
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
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue;
                    B[rowNodeId] +=
                        kValue * U[colNodeId] +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]);
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
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue;
                    B[rowNodeId] +=
                        kValue * U[colNodeId + d2Offset] +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]);
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
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue;
                    B[rowNodeId] +=
                        kValue * U[colNodeId + rOffset] +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]);
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
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue;
                    B[rowNodeId + d2Offset] +=
                        kValue * U[colNodeId] +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]);
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
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue;
                    B[rowNodeId + d2Offset] +=
                        kValue * U[colNodeId + d2Offset] +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]);
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
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue;
                    B[rowNodeId + d2Offset] +=
                        kValue * U[colNodeId + rOffset] +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]);
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
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue;
                    B[rowNodeId + rOffset] +=
                        kValue * U[colNodeId] +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]);
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
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue;
                    B[rowNodeId + rOffset] +=
                        kValue * U[colNodeId + d2Offset] +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]);
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
                        kValue +
                        (1.0 / (beta * dt * dt)) * mValue;
                    B[rowNodeId + rOffset] +=
                        kValue * U[colNodeId + rOffset] +
                        mValue * (
                        (1.0 / (beta * dt * dt)) * u[colDof] +
                        (1.0 / (beta * dt)) * vel[colDof] +
                        ((1.0 / (2.0 * beta)) - 1.0) * acc[colDof]);
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
                B[rowNodeId] += -f[row * localDof];
            }
            for (int row = 0; row < d2ElemNodeCnt; row++)
            {
                int rowNodeId = d2Nodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                B[rowNodeId + d2Offset] += -f[row * localDof + localD2Offset];
            }
            for (int row = 0; row < rElemNodeCnt; row++)
            {
                int rowNodeId = rNodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                B[rowNodeId + rOffset] += -f[row * localDof + localROffset];
            }
        }
    }
}
