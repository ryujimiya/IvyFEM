using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DFEM
    {
        protected void CalcSimpleCorotationalFrameKl(
            double l0, double E, double Ae, double Iz,
            double barU, double barT1, double barT2,
            out double[] fl, out IvyFEM.Lapack.DoubleMatrix kl)
        {
            Elastic2DFEMUtils.CalcSimpleCorotationalFrameKl(
                l0, E, Ae, Iz,
                barU, barT1, barT2,
                out fl, out kl);
        }

        protected void CalcBernoulliCorotationalFrameKl(
            double l0, double E, double Ae, double Iz,
            double barU, double barT1, double barT2,
            out double[] fl, out IvyFEM.Lapack.DoubleMatrix kl)
        {
            Elastic2DFEMUtils.CalcBernoulliCorotationalFrameKl(
                l0, E, Ae, Iz,
                barU, barT1, barT2,
                out fl, out kl);
        }

        protected void CalcShallowArchCorotationalFrameKl(
            double l0, double E, double Ae, double Iz,
            double barU, double barT1, double barT2,
            out double[] fl, out IvyFEM.Lapack.DoubleMatrix kl)
        {
            Elastic2DFEMUtils.CalcShallowArchCorotationalFrameKl(
                l0, E, Ae, Iz,
                barU, barT1, barT2,
                out fl, out kl);
        }

        protected void CalcCorotationalFrameElementABForLine(
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
                if (!(workMa0 is CorotationalFrameMaterial))
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
            System.Diagnostics.Debug.Assert(ma0 is CorotationalFrameMaterial);

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

            System.Diagnostics.Debug.Assert(d1ElemNodeCnt == 2);
            System.Diagnostics.Debug.Assert(d2ElemNodeCnt == d1ElemNodeCnt);
            System.Diagnostics.Debug.Assert(rElemNodeCnt == d1ElemNodeCnt);
            uint elemNodeCnt = d1ElemNodeCnt;
            double[] pg = new double[elemNodeCnt * 3];
            for (int iNode = 0; iNode < d1ElemNodeCnt; iNode++)
            {
                int d1NodeId = d1Nodes[iNode];
                if (d1NodeId == -1)
                {
                    continue;
                }
                pg[iNode * 3] = U[d1NodeId];
            }
            for (int iNode = 0; iNode < d2ElemNodeCnt; iNode++)
            {
                int d2NodeId = d2Nodes[iNode];
                if (d2NodeId == -1)
                {
                    continue;
                }
                pg[iNode * 3 + 1] = U[d2NodeId + d2Offset];
            }
            for (int iNode = 0; iNode < rElemNodeCnt; iNode++)
            {
                int rNodeId = rNodes[iNode];
                if (rNodeId == -1)
                {
                    continue;
                }
                pg[iNode * 3 + 2] = U[rNodeId + rOffset];
            }
            double u1 = pg[0];
            double v1 = pg[1];
            double t1 = pg[2];
            double u2 = pg[3];
            double v2 = pg[4];
            double t2 = pg[5];

            var ma = ma0 as CorotationalFrameMaterial;
            double Ae = ma.Area;
            double Iz = ma.SecondMomentOfArea;
            double rho = ma.MassDensity;
            double E = ma.Young;

            double[] pt10 = vCoords[0];
            double[] pt20 = vCoords[1];
            double[] pt1n = { pt10[0] + u1, pt10[1] + v1 };
            double[] pt2n = { pt20[0] + u2, pt20[1] + v2 };
            double l0 = Math.Sqrt(
                (pt20[0] - pt10[0]) * (pt20[0] - pt10[0]) +
                (pt20[1] - pt10[1]) * (pt20[1] - pt10[1]));
            double ln = Math.Sqrt(
                (pt2n[0] - pt1n[0]) * (pt2n[0] - pt1n[0]) +
                (pt2n[1] - pt1n[1]) * (pt2n[1] - pt1n[1]));
            double c0 = (pt20[0] - pt10[0]) / l0;
            double s0 = (pt20[1] - pt10[1]) / l0;
            double c = (pt2n[0] - pt1n[0]) / ln;
            double s = (pt2n[1] - pt1n[1]) / ln;
            double sa = c0 * s - s0 * c;
            double ca = c0 * c + s0 * s;
            double alpha = 0.0;
            if (sa >= 0.0)
            {
                if (ca >= 0.0)
                {
                    alpha = Math.Asin(sa);
                }
                else
                {
                    alpha = Math.Acos(ca);
                }
            }
            else
            {
                if (ca >= 0.0)
                {
                    alpha = Math.Asin(sa);
                }
                else
                {
                    alpha = -Math.Acos(ca);
                }
            }
            System.Diagnostics.Debug.Assert(Math.Abs(alpha) <= Math.PI);

            double barU = ln - l0;
            double barT1 = t1 - alpha;
            double barT2 = t2 - alpha;

            double[][] bVec = new double[3][];
            bVec[0] = new double[] { -c, -s, 0.0, c, s, 0.0 };
            bVec[1] = new double[] { -s / ln, c / ln, 1.0, s / ln, -c / ln, 0.0 };
            bVec[2] = new double[] { -s / ln, c / ln, 0.0, s / ln, -c / ln, 1.0 };
            var bb = new IvyFEM.Lapack.DoubleMatrix(3, 6);
            for (int row = 0; row < 3; row++)
            {
                for (int col = 0; col < 6; col++)
                {
                    bb[row, col] = bVec[row][col];
                }
            }
            var transbb = IvyFEM.Lapack.DoubleMatrix.Transpose(bb);
            double[] rVec = { -c, -s, 0.0, c, s, 0.0 };
            double[] zVec = { s, -c, 0.0, -s, c, 0.0 };
            var r = new IvyFEM.Lapack.DoubleMatrix(6, 1);
            var z = new IvyFEM.Lapack.DoubleMatrix(6, 1);
            for (int i = 0; i < 6; i++)
            {
                r[i, 0] = rVec[i];
                z[i, 0] = zVec[i];
            }
            var transr = IvyFEM.Lapack.DoubleMatrix.Transpose(r);
            var transz = IvyFEM.Lapack.DoubleMatrix.Transpose(z);

            //---------------------------------
            // local
            double[] fl;
            IvyFEM.Lapack.DoubleMatrix kl;

            // 
            //CalcSimpleCorotationalFrameKl(l0, E, Ae, Iz, barU, barT1, barT2, out fl, out kl);
            // Bernoulli
            //CalcBernoulliCorotationalFrameKl(l0, E, Ae, Iz, barU, barT1, barT2, out fl, out kl);
            // shallow arch beam
            CalcShallowArchCorotationalFrameKl(l0, E, Ae, Iz, barU, barT1, barT2, out fl, out kl);
            
            double n = fl[0];
            double m1 = fl[1];
            double m2 = fl[2];
            //---------------------------------

            //---------------------------------
            // f & K
            double[] f = transbb * fl;

            var Ke = new IvyFEM.Lapack.DoubleMatrix(6, 6);
            {
                var work1 = (transbb * kl) * bb;
                var work2 = z * transz;
                work2 = IvyFEM.Lapack.DoubleMatrix.Scal(work2, n / ln);
                var work31 = r * transz;
                var work32 = z * transr;
                var work3 = new IvyFEM.Lapack.DoubleMatrix(6, 6);
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        work3[i, j] = (1.0 / (ln * ln)) * (work31[i, j] + work32[i, j])  * (m1 + m2);
                    }
                }
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        Ke[i, j] = work1[i, j] + work2[i, j] + work3[i, j];
                    }
                }
            }
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
                    double kValue = Ke[row * localDof, col * localDof];
                    A[rowNodeId, colNodeId] += kValue;
                    B[rowNodeId] += kValue * U[colNodeId];
                }
                // displacement 2
                for (int col = 0; col < d2ElemNodeCnt; col++)
                {
                    int colNodeId = d2Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    double kValue = Ke[row * localDof, col * localDof + localD2Offset];
                    A[rowNodeId, colNodeId + d2Offset] += kValue;
                    B[rowNodeId] += kValue * U[colNodeId + d2Offset];
                }
                // rotation
                for (int col = 0; col < rElemNodeCnt; col++)
                {
                    int colNodeId = rNodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    double kValue = Ke[row * localDof, col * localDof + localROffset];
                    A[rowNodeId, colNodeId + rOffset] += kValue;
                    B[rowNodeId] += kValue * U[colNodeId + rOffset];
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
                    double kValue = Ke[row * localDof + localD2Offset, col * localDof];
                    A[rowNodeId + d2Offset, colNodeId] += kValue;
                    B[rowNodeId + d2Offset] += kValue * U[colNodeId];
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
                    A[rowNodeId + d2Offset, colNodeId + d2Offset] += kValue;
                    B[rowNodeId + d2Offset] += kValue * U[colNodeId + d2Offset];
                }
                // rotation
                for (int col = 0; col < rElemNodeCnt; col++)
                {
                    int colNodeId = rNodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    double kValue = Ke[row * localDof + localD2Offset, col * localDof + localROffset];
                    A[rowNodeId + d2Offset, colNodeId + rOffset] += kValue;
                    B[rowNodeId + d2Offset] += kValue * U[colNodeId + rOffset];
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
                    double kValue = Ke[row * localDof + localROffset, col * localDof];
                    A[rowNodeId + rOffset, colNodeId] += kValue;
                    B[rowNodeId + rOffset] += kValue * U[colNodeId];
                }
                // displacement 2
                for (int col = 0; col < d2ElemNodeCnt; col++)
                {
                    int colNodeId = d2Nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    double kValue = Ke[row * localDof + localROffset, col * localDof + localD2Offset];
                    A[rowNodeId + rOffset, colNodeId + d2Offset] += kValue;
                    B[rowNodeId + rOffset] += kValue * U[colNodeId + d2Offset];
                }
                // rotation
                for (int col = 0; col < rElemNodeCnt; col++)
                {
                    int colNodeId = rNodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    double kValue = Ke[row * localDof + localROffset, col * localDof + localROffset];
                    A[rowNodeId + rOffset, colNodeId + rOffset] += kValue;
                    B[rowNodeId + rOffset] += kValue * U[colNodeId + rOffset];
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
