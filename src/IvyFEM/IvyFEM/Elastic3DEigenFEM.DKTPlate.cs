using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Animation;

namespace IvyFEM
{
    partial class Elastic3DEigenFEM
    {
        protected IvyFEM.Lapack.DoubleMatrix CalcDKTPlateKl(
            TriangleFE d1TriFE, double h, IvyFEM.Lapack.DoubleMatrix C)
        {
            return Elastic3DFEMUtils.CalcDKTPlateKl(d1TriFE, h, C);
        }

        protected IvyFEM.Lapack.DoubleMatrix CalcDKTPlateMl(
            TriangleFE d1TriFE, double h, double rho)
        {
            return Elastic3DFEMUtils.CalcDKTPlateMl(d1TriFE, h, rho);
        }

        protected void CalcDKTPlateElementKM(
            uint feId, IvyFEM.Lapack.DoubleMatrix K, IvyFEM.Lapack.DoubleMatrix M)
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
                if (!(workMa0 is DKTPlateMaterial))
                {
                    return;
                }
            }

            System.Diagnostics.Debug.Assert(DisplacementQuantityIds.Count == 2);
            System.Diagnostics.Debug.Assert(DisplacementQuantityIds[0] == 0);
            System.Diagnostics.Debug.Assert(DisplacementQuantityIds[1] == 1);
            uint d1QuantityId = 0; // displacement
            uint d2QuantityId = 1; // displacement
            uint r1QuantityId = 2; // rotation
            uint r2QuantityId = 1; // rotation
            System.Diagnostics.Debug.Assert(World.GetDof(d1QuantityId) == 2);
            System.Diagnostics.Debug.Assert(World.GetDof(d2QuantityId) == 1);
            System.Diagnostics.Debug.Assert(World.GetDof(r1QuantityId) == 2);
            System.Diagnostics.Debug.Assert(World.GetDof(r2QuantityId) == 1);
            int d1Dof = 2; // u v
            int d2Dof = 1; // w
            int r1Dof = 2; //θx θy
            int r2Dof = 1; //θz
            int d1NodeCnt = (int)World.GetNodeCount(d1QuantityId);
            int d2NodeCnt = (int)World.GetNodeCount(d2QuantityId);
            int r1NodeCnt = (int)World.GetNodeCount(r1QuantityId);
            int r2NodeCnt = (int)World.GetNodeCount(r2QuantityId);
            int d2Offset = d1NodeCnt * d1Dof;
            int r1Offset = d2Offset + d2NodeCnt * d2Dof;
            int r2Offset = r1Offset + r1NodeCnt * r1Dof;

            TriangleFE d1TriFE = World.GetTriangleFE(d1QuantityId, feId);
            TriangleFE d2TriFE = World.GetTriangleFE(d2QuantityId, feId);
            TriangleFE r1TriFE = World.GetTriangleFE(r1QuantityId, feId);
            TriangleFE r2TriFE = World.GetTriangleFE(r2QuantityId, feId);
            uint maId = d1TriFE.MaterialId;
            if (!World.IsMaterialId(maId))
            {
                return;
            }
            Material ma0 = World.GetMaterial(maId);
            System.Diagnostics.Debug.Assert(ma0 is DKTPlateMaterial);

            int[] d1CoIds = d1TriFE.NodeCoordIds;
            uint d1ElemNodeCnt = d1TriFE.NodeCount;
            // 暫定で要素マトリクスをこのメソッド内に直接記述
            System.Diagnostics.Debug.Assert(d1ElemNodeCnt == 3); // 1次要素
            int[] d1Nodes = new int[d1ElemNodeCnt];
            for (int iNode = 0; iNode < d1ElemNodeCnt; iNode++)
            {
                int coId = d1CoIds[iNode];
                int nodeId = World.Coord2Node(d1QuantityId, coId);
                d1Nodes[iNode] = nodeId;
            }
            int[] d2CoIds = d2TriFE.NodeCoordIds;
            uint d2ElemNodeCnt = d2TriFE.NodeCount;
            // 暫定で要素マトリクスをこのメソッド内に直接記述
            System.Diagnostics.Debug.Assert(d2ElemNodeCnt == 3); // 1次要素
            int[] d2Nodes = new int[d2ElemNodeCnt];
            for (int iNode = 0; iNode < d2ElemNodeCnt; iNode++)
            {
                int coId = d2CoIds[iNode];
                int nodeId = World.Coord2Node(d2QuantityId, coId);
                d2Nodes[iNode] = nodeId;
            }
            int[] r1CoIds = r1TriFE.NodeCoordIds;
            uint r1ElemNodeCnt = r1TriFE.NodeCount;
            // 暫定で要素マトリクスをこのメソッド内に直接記述
            System.Diagnostics.Debug.Assert(r1ElemNodeCnt == 3); // 1次要素
            int[] r1Nodes = new int[r1ElemNodeCnt];
            for (int iNode = 0; iNode < r1ElemNodeCnt; iNode++)
            {
                int coId = r1CoIds[iNode];
                int nodeId = World.Coord2Node(r1QuantityId, coId);
                r1Nodes[iNode] = nodeId;
            }
            int[] r2CoIds = r2TriFE.NodeCoordIds;
            uint r2ElemNodeCnt = r2TriFE.NodeCount;
            // 暫定で要素マトリクスをこのメソッド内に直接記述
            System.Diagnostics.Debug.Assert(r1ElemNodeCnt == 3); // 1次要素
            int[] r2Nodes = new int[r2ElemNodeCnt];
            for (int iNode = 0; iNode < r2ElemNodeCnt; iNode++)
            {
                int coId = r2CoIds[iNode];
                int nodeId = World.Coord2Node(r2QuantityId, coId);
                r2Nodes[iNode] = nodeId;
            }

            var ma = ma0 as DKTPlateMaterial;
            double h = ma.Thickness;
            double rho = ma.MassDensity;
            double E = ma.Young;
            double nu = ma.Poisson;
            double Ks = ma.ShearCorrectionFactor;
            IvyFEM.Lapack.DoubleMatrix C = new IvyFEM.Lapack.DoubleMatrix(3, 3);
            {
                double k1 = E / (1.0 - nu * nu);
                C[0, 0] = k1 * 1.0;
                C[0, 1] = k1 * nu;
                C[0, 2] = k1 * 0.0;
                C[1, 1] = k1 * 1.0;
                C[1, 2] = k1 * 0.0;
                C[2, 2] = k1 * (1.0 / 2.0) * (1.0 - nu);
                for (int i = 0; i < 3; i++)
                {
                    for (int j = i + 1; j < 3; j++)
                    {
                        C[j, i] = C[i, j];
                    }
                }
            }

            //------------------------------------------
            // local
            var Kl = CalcDKTPlateKl(d1TriFE, h, C);
            var Ml = CalcDKTPlateMl(d1TriFE, h, rho);

            //-------------------------------------------
            // global
            IvyFEM.Lapack.DoubleMatrix T = new Lapack.DoubleMatrix(18, 18);
            {
                double[] co1 = World.GetCoord(d1QuantityId, d1CoIds[0]);
                double[] co2 = World.GetCoord(d1QuantityId, d1CoIds[1]);
                double[] co3 = World.GetCoord(d1QuantityId, d1CoIds[2]);
                OpenTK.Vector3d v1 = new OpenTK.Vector3d(co1[0], co1[1], co1[2]);
                OpenTK.Vector3d v2 = new OpenTK.Vector3d(co2[0], co2[1], co2[2]);
                OpenTK.Vector3d v3 = new OpenTK.Vector3d(co3[0], co3[1], co3[2]);
                var vx = v2 - v1;
                var vr = v3 - v1;
                var vz = OpenTK.Vector3d.Cross(vx, vr);
                var ex = OpenTK.Vector3d.Normalize(vx);
                var ez = OpenTK.Vector3d.Normalize(vz);
                var ey = OpenTK.Vector3d.Cross(ez, ex);
                double[] doubleEx = { ex.X, ex.Y, ex.Z };
                double[] doubleEy = { ey.X, ey.Y, ey.Z };
                double[] doubleEz = { ez.X, ez.Y, ez.Z };
                IvyFEM.Lapack.DoubleMatrix T3 = new IvyFEM.Lapack.DoubleMatrix(3, 3);
                for (int j = 0; j < 3; j++)
                {
                    T3[0, j] = doubleEx[j];
                    T3[1, j] = doubleEy[j];
                    T3[2, j] = doubleEz[j];
                }
                for (int iT3 = 0; iT3 < 6; iT3++)
                {
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            T[i + 3 * iT3, j + 3 * iT3] = T3[i, j];
                        }
                    }
                }
            }
            var transT = IvyFEM.Lapack.DoubleMatrix.Transpose(T);
            var Ke = transT * Kl * T;
            var Me = transT * Ml * T;
            //--------------------------------

            // local dof
            int localDof = 6;
            System.Diagnostics.Debug.Assert(localDof == (d1Dof + d2Dof + r1Dof + r2Dof));
            int localD2Offset = d1Dof;
            int localR1Offset = localD2Offset + d2Dof;
            int localR2Offset = localR1Offset + r1Dof;
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
                    for (int rowDof = 0; rowDof < d1Dof; rowDof++)
                    {
                        for (int colDof = 0; colDof < d1Dof; colDof++)
                        {
                            double kValue = Ke[row * localDof + rowDof, col * localDof + colDof];
                            double mValue = Me[row * localDof + rowDof, col * localDof + colDof];
                            K[rowNodeId * d1Dof + rowDof, colNodeId * d1Dof + colDof] += kValue;
                            M[rowNodeId * d1Dof + rowDof, colNodeId * d1Dof + colDof] += mValue;
                        }
                    }
                }
                // displacement 2
                for (int col = 0; col < d2ElemNodeCnt; col++)
                {
                    int colNodeId = d2Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    for (int rowDof = 0; rowDof < d1Dof; rowDof++)
                    {
                        double kValue = Ke[row * localDof + rowDof, col * localDof + localD2Offset];
                        double mValue = Me[row * localDof + rowDof, col * localDof + localD2Offset];
                        K[rowNodeId * d1Dof + rowDof, colNodeId + d2Offset] += kValue;
                        M[rowNodeId * d1Dof + rowDof, colNodeId + d2Offset] += mValue;
                    }
                }
                // rotation 1
                for (int col = 0; col < r1ElemNodeCnt; col++)
                {
                    int colNodeId = r1Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    for (int rowDof = 0; rowDof < d1Dof; rowDof++)
                    {
                        for (int colDof = 0; colDof < r1Dof; colDof++)
                        {
                            double kValue = Ke[row * localDof + rowDof, col * localDof + colDof + localR1Offset];
                            double mValue = Me[row * localDof + rowDof, col * localDof + colDof + localR1Offset];
                            K[rowNodeId * d1Dof + rowDof, colNodeId * r1Dof + colDof + r1Offset] += kValue;
                            M[rowNodeId * d1Dof + rowDof, colNodeId * r1Dof + colDof + r1Offset] += mValue;
                        }
                    }
                }
                // rotation 2
                for (int col = 0; col < r2ElemNodeCnt; col++)
                {
                    int colNodeId = r2Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    for (int rowDof = 0; rowDof < d1Dof; rowDof++)
                    {
                        double kValue = Ke[row * localDof + rowDof, col * localDof + localR2Offset];
                        double mValue = Me[row * localDof + rowDof, col * localDof + localR2Offset];
                        K[rowNodeId * d1Dof + rowDof, colNodeId + r2Offset] += kValue;
                        M[rowNodeId * d1Dof + rowDof, colNodeId + r2Offset] += mValue;
                    }
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
                    for (int colDof = 0; colDof < d1Dof; colDof++)
                    {
                        double kValue = Ke[row * localDof + localD2Offset, col * localDof + colDof];
                        double mValue = Me[row * localDof + localD2Offset, col * localDof + colDof];
                        K[rowNodeId + d2Offset, colNodeId * d1Dof + colDof] += kValue;
                        M[rowNodeId + d2Offset, colNodeId * d1Dof + colDof] += mValue;
                    }
                }
                // displacement 2
                for (int col = 0; col < d2ElemNodeCnt; col++)
                {
                    int colNodeId = d2Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    double kValue = Ke[row * localDof + localD2Offset, col * localDof + localD2Offset];
                    double mValue = Me[row * localDof + localD2Offset, col * localDof + localD2Offset];
                    K[rowNodeId + d2Offset, colNodeId + d2Offset] += kValue;
                    M[rowNodeId + d2Offset, colNodeId + d2Offset] += mValue;
                }
                // rotation 1
                for (int col = 0; col < r1ElemNodeCnt; col++)
                {
                    int colNodeId = r1Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    for (int colDof = 0; colDof < r1Dof; colDof++)
                    {
                        double kValue = Ke[row * localDof + localD2Offset, col * localDof + colDof + localR1Offset];
                        double mValue = Me[row * localDof + localD2Offset, col * localDof + colDof + localR1Offset];
                        K[rowNodeId + d2Offset, colNodeId * r1Dof + colDof + r1Offset] += kValue;
                        M[rowNodeId + d2Offset, colNodeId * r1Dof + colDof + r1Offset] += mValue;
                    }
                }
                // rotation 2
                for (int col = 0; col < r2ElemNodeCnt; col++)
                {
                    int colNodeId = r2Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    double kValue = Ke[row * localDof + localD2Offset, col * localDof + localR2Offset];
                    double mValue = Me[row * localDof + localD2Offset, col * localDof + localR2Offset];
                    K[rowNodeId + d2Offset, colNodeId + r2Offset] += kValue;
                    M[rowNodeId + d2Offset, colNodeId + r2Offset] += mValue;
                }
            }
            // rotation 1
            for (int row = 0; row < r1ElemNodeCnt; row++)
            {
                int rowNodeId = r1Nodes[row];
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
                    for (int rowDof = 0; rowDof < r1Dof; rowDof++)
                    {
                        for (int colDof = 0; colDof < d1Dof; colDof++)
                        {
                            double kValue = Ke[row * localDof + rowDof + localR1Offset, col * localDof + colDof];
                            double mValue = Me[row * localDof + rowDof + localR1Offset, col * localDof + colDof];
                            K[rowNodeId * r1Dof + rowDof + r1Offset, colNodeId * d1Dof + colDof] += kValue;
                            M[rowNodeId * r1Dof + rowDof + r1Offset, colNodeId * d1Dof + colDof] += mValue;
                        }
                    }
                }
                // displacement 2
                for (int col = 0; col < d2ElemNodeCnt; col++)
                {
                    int colNodeId = d2Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    for (int rowDof = 0; rowDof < r1Dof; rowDof++)
                    {
                        double kValue = Ke[row * localDof + rowDof + localR1Offset, col * localDof + localD2Offset];
                        double mValue = Me[row * localDof + rowDof + localR1Offset, col * localDof + localD2Offset];
                        K[rowNodeId * r1Dof + rowDof + r1Offset, colNodeId + d2Offset] += kValue;
                        M[rowNodeId * r1Dof + rowDof + r1Offset, colNodeId + d2Offset] += mValue;
                    }
                }
                // rotation 1
                for (int col = 0; col < r1ElemNodeCnt; col++)
                {
                    int colNodeId = r1Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    for (int rowDof = 0; rowDof < r1Dof; rowDof++)
                    {
                        for (int colDof = 0; colDof < r1Dof; colDof++)
                        {
                            double kValue = 
                                Ke[row * localDof + rowDof + localR1Offset, col * localDof + colDof + localR1Offset];
                            double mValue =
                                Me[row * localDof + rowDof + localR1Offset, col * localDof + colDof + localR1Offset];
                            K[rowNodeId * r1Dof + rowDof + r1Offset, colNodeId * r1Dof + colDof + r1Offset] += kValue;
                            M[rowNodeId * r1Dof + rowDof + r1Offset, colNodeId * r1Dof + colDof + r1Offset] += mValue;
                        }
                    }
                }
                // rotation 2
                for (int col = 0; col < r2ElemNodeCnt; col++)
                {
                    int colNodeId = r2Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    for (int rowDof = 0; rowDof < r1Dof; rowDof++)
                    {
                        double kValue = Ke[row * localDof + rowDof + localR1Offset, col * localDof + localR2Offset];
                        double mValue = Me[row * localDof + rowDof + localR1Offset, col * localDof + localR2Offset];
                        K[rowNodeId * r1Dof + rowDof + r1Offset, colNodeId + r2Offset] += kValue;
                        M[rowNodeId * r1Dof + rowDof + r1Offset, colNodeId + r2Offset] += mValue;
                    }
                }
            }
            // rotation 2
            for (int row = 0; row < r2ElemNodeCnt; row++)
            {
                int rowNodeId = r2Nodes[row];
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
                    for (int colDof = 0; colDof < d1Dof; colDof++)
                    {
                        double kValue = Ke[row * localDof + localR2Offset, col * localDof + colDof];
                        double mValue = Me[row * localDof + localR2Offset, col * localDof + colDof];
                        K[rowNodeId + r2Offset, colNodeId * d1Dof + colDof] += kValue;
                        M[rowNodeId + r2Offset, colNodeId * d1Dof + colDof] += mValue;
                    }
                }
                // displacement 2
                for (int col = 0; col < d2ElemNodeCnt; col++)
                {
                    int colNodeId = d2Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    double kValue = Ke[row * localDof + localR2Offset, col * localDof + localD2Offset];
                    double mValue = Me[row * localDof + localR2Offset, col * localDof + localD2Offset];
                    K[rowNodeId + r2Offset, colNodeId + d2Offset] += kValue;
                    M[rowNodeId + r2Offset, colNodeId + d2Offset] += mValue;
                }
                // rotation 1
                for (int col = 0; col < r1ElemNodeCnt; col++)
                {
                    int colNodeId = r1Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    for (int colDof = 0; colDof < r1Dof; colDof++)
                    {
                        double kValue = Ke[row * localDof + localR2Offset, col * localDof + colDof + localR1Offset];
                        double mValue = Me[row * localDof + localR2Offset, col * localDof + colDof + localR1Offset];
                        K[rowNodeId + r2Offset, colNodeId * r1Dof + colDof + r1Offset] += kValue;
                        M[rowNodeId + r2Offset, colNodeId * r1Dof + colDof + r1Offset] += mValue;
                    }
                }
                // rotation 2
                for (int col = 0; col < r2ElemNodeCnt; col++)
                {
                    int colNodeId = r2Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    double kValue = Ke[row * localDof + localR2Offset, col * localDof + localR2Offset];
                    double mValue = Me[row * localDof + localR2Offset, col * localDof + localR2Offset];
                    K[rowNodeId + r2Offset, colNodeId + r2Offset] += kValue;
                    M[rowNodeId + r2Offset, colNodeId + r2Offset] += mValue;
                }
            }
        }
    }
}
