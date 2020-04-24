using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DTDFEM
    {
        protected void CalcSimpleCorotationalFrameKl(
            double E, double Ae, double Iz,
            double l0, double barU, double barT1, double barT2,
            out double[] fl, out IvyFEM.Lapack.DoubleMatrix kl)
        {
            Elastic2DFEMUtils.CalcSimpleCorotationalFrameKl(
                E, Ae, Iz,
                l0, barU, barT1, barT2,
                out fl, out kl);
        }

        protected void CalcBernoulliCorotationalFrameKl(
            double E, double Ae, double Iz,
            double l0, double barU, double barT1, double barT2,
            out double[] fl, out IvyFEM.Lapack.DoubleMatrix kl)
        {
            Elastic2DFEMUtils.CalcBernoulliCorotationalFrameKl(
                E, Ae, Iz,
                l0, barU, barT1, barT2,
                out fl, out kl);
        }

        protected void CalcShallowArchCorotationalFrameKl(
            double E, double Ae, double Iz,
            double l0, double barU, double barT1, double barT2,
            out double[] fl, out IvyFEM.Lapack.DoubleMatrix kl)
        {
            Elastic2DFEMUtils.CalcShallowArchCorotationalFrameKl(
                E, Ae, Iz,
                l0, barU, barT1, barT2,
                out fl, out kl);
        }

        protected IvyFEM.Lapack.DoubleMatrix CalcCorotationalFrameMl(
            double rho, double Ae, double Ix, double l0)
        {
            return Elastic2DFEMUtils.CalcCorotationalFrameMl(
                rho, Ae, Ix, l0);
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
            double Ix = ma.PolarSecondMomentOfArea;
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
            //CalcSimpleCorotationalFrameKl(E, Ae, Iz, l0, barU, barT1, barT2, out fl, out kl);
            // Bernoulli
            //CalcBernoulliCorotationalFrameKl(E, Ae, Iz, l0, barU, barT1, barT2, out fl, out kl);
            // shallow arch beam
            CalcShallowArchCorotationalFrameKl(E, Ae, Iz, l0, barU, barT1, barT2, out fl, out kl);
            
            double n = fl[0];
            double m1 = fl[1];
            double m2 = fl[2];
            //---------------------------------

            //---------------------------------
            // f & K
            double[] f = transbb * fl;

            var K = new IvyFEM.Lapack.DoubleMatrix(6, 6);
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
                        K[i, j] = work1[i, j] + work2[i, j] + work3[i, j];
                    }
                }
            }
            //---------------------------------

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
            double[] prevPg = { prevU1, prevV1, prevT1, prevU2, prevV2, prevT2 };
            double[] prevVelPg = { prevVelU1, prevVelV1, prevVelT1, prevVelU2, prevVelV2, prevVelT2 };
            double[] prevAccPg = { prevAccU1, prevAccV1, prevAccT1, prevAccU2, prevAccV2, prevAccT2 }; 
            var velPg = new double[6];
            var accPg = new double[6];
            for (int i = 0; i < 6; i++)
            {
                velPg[i] =
                    (gamma / (beta * dt)) * (pg[i] - prevPg[i]) +
                    (1.0 - gamma / beta) * prevVelPg[i] +
                    dt * (1.0 - gamma / (2.0 * beta)) * prevAccPg[i];
                accPg[i] =
                    (1.0 / (beta * dt * dt)) * (pg[i] - prevPg[i]) -
                    (1.0 / (beta * dt)) * prevVelPg[i] -
                    (1.0 / (2.0 * beta) - 1.0) * prevAccPg[i];
            }

            var pgM = new IvyFEM.Lapack.DoubleMatrix(6, 1);
            pg.CopyTo(pgM.Buffer, 0);
            var velPgM = new IvyFEM.Lapack.DoubleMatrix(6, 1);
            velPg.CopyTo(velPgM.Buffer, 0);
            var accPgM = new IvyFEM.Lapack.DoubleMatrix(6, 1);
            accPg.CopyTo(accPgM.Buffer, 0);
            var transvelPgM = IvyFEM.Lapack.DoubleMatrix.Transpose(velPgM);

            //-----------------------------
            // local
            var ml = CalcCorotationalFrameMl(rho, Ae, Ix, l0);
            //---------------------------------

            //---------------------------------
            // T, M
            double cosb = (pt2n[0] - pt1n[0]) / ln;
            double sinb = (pt2n[1] - pt1n[1]) / ln;
            var T = new IvyFEM.Lapack.DoubleMatrix(6, 6);
            T[0, 0] = cosb;
            T[0, 1] = sinb;
            T[0, 2] = 0.0;
            T[0, 3] = 0.0;
            T[0, 4] = 0.0;
            T[0, 5] = 0.0;
            T[1, 0] = -sinb;
            T[1, 1] = cosb;
            T[1, 2] = 0.0;
            T[1, 3] = 0.0;
            T[1, 4] = 0.0;
            T[1, 5] = 0.0;
            T[2, 0] = 0.0;
            T[2, 1] = 0.0;
            T[2, 2] = 1.0;
            T[2, 3] = 0.0;
            T[2, 4] = 0.0;
            T[2, 5] = 0.0;
            T[3, 1] = 0.0;
            T[3, 2] = 0.0;
            T[3, 3] = cosb;
            T[3, 4] = sinb;
            T[3, 5] = 0.0;
            T[4, 0] = 0.0;
            T[4, 1] = 0.0;
            T[4, 2] = 0.0;
            T[4, 3] = -sinb;
            T[4, 4] = cosb;
            T[4, 5] = 0.0;
            T[5, 0] = 0.0;
            T[5, 1] = 0.0;
            T[5, 2] = 0.0;
            T[5, 3] = 0.0;
            T[5, 4] = 0.0;
            T[5, 5] = 1.0;
            var transT = IvyFEM.Lapack.DoubleMatrix.Transpose(T);
            var M = (transT * ml) * T;

            var I0 = new IvyFEM.Lapack.DoubleMatrix(6, 6);
            I0[0, 0] = 0.0;
            I0[0, 1] = 1.0;
            I0[0, 2] = 0.0;
            I0[0, 3] = 0.0;
            I0[0, 4] = 0.0;
            I0[0, 5] = 0.0;
            I0[1, 0] = -1.0;
            I0[1, 1] = 0.0;
            I0[1, 2] = 0.0;
            I0[1, 3] = 0.0;
            I0[1, 4] = 0.0;
            I0[1, 5] = 0.0;
            I0[2, 0] = 0.0;
            I0[2, 1] = 0.0;
            I0[2, 2] = 0.0;
            I0[2, 3] = 0.0;
            I0[2, 4] = 0.0;
            I0[2, 5] = 0.0;
            I0[3, 0] = 0.0;
            I0[3, 1] = 0.0;
            I0[3, 2] = 0.0;
            I0[3, 3] = 0.0;
            I0[3, 4] = 1.0;
            I0[3, 5] = 0.0;
            I0[4, 0] = 0.0;
            I0[4, 1] = 0.0;
            I0[4, 2] = 0.0;
            I0[4, 3] = -1.0;
            I0[4, 4] = 0.0;
            I0[4, 5] = 0.0;
            I0[5, 0] = 0.0;
            I0[5, 1] = 0.0;
            I0[5, 2] = 0.0;
            I0[5, 3] = 0.0;
            I0[5, 4] = 0.0;
            I0[5, 5] = 0.0;
            var transI0 = IvyFEM.Lapack.DoubleMatrix.Transpose(I0);
            IvyFEM.Lapack.DoubleMatrix Mb;
            {
                var work1 = (transI0 * ml);
                var work2 = ml * I0;
                var work3 = new IvyFEM.Lapack.DoubleMatrix(6, 6);
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        work3[i, j] = work1[i, j] + work2[i, j];
                    }
                }
                Mb = (transT * work3) * T;
            }
            IvyFEM.Lapack.DoubleMatrix dotM;
            {
                double d = IvyFEM.Lapack.Functions.ddot(z.Buffer, velPg);
                dotM = IvyFEM.Lapack.DoubleMatrix.Scal(Mb, d / ln);
            }
            double[] fk = new double[6];
            {
                var work1 = M * accPgM;
                var work2 = dotM * velPgM;
                double d3 = IvyFEM.Lapack.Functions.ddot((transvelPgM * Mb).Buffer, velPg);
                var work3 = new IvyFEM.Lapack.DoubleMatrix(6, 1);
                for (int i = 0; i < 6; i++)
                {
                    work3[i, 0] = (1.0 / 2.0) * d3 * (1.0 / ln) * z[i, 0];
                }
                for (int i = 0; i < 6; i++)
                {
                    fk[i] = work1[i, 0] + work2[i, 0] - work3[i, 0];
                }
            }

            IvyFEM.Lapack.DoubleMatrix Ck1;
            {
                var work = velPgM * transz;
                work = IvyFEM.Lapack.DoubleMatrix.Scal(work, (1.0 / ln));
                Ck1 = Mb * work;
            }
            var transCk1 = IvyFEM.Lapack.DoubleMatrix.Transpose(Ck1);
            var Ck = new IvyFEM.Lapack.DoubleMatrix(6, 6);
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    Ck[i, j] = dotM[i, j] + Ck1[i, j] - transCk1[i, j];
                }
            }

            IvyFEM.Lapack.DoubleMatrix Mbb = new IvyFEM.Lapack.DoubleMatrix(6, 6);
            {
                var work11 = transI0 * ml;
                var work12 = ml * I0;
                var work1 = new IvyFEM.Lapack.DoubleMatrix(6, 6);
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        work1[i, j] = work11[i, j] + work12[i, j];
                    }
                }
                var work2 = ((transT * transI0) * work1) * T;
                var work3 = ((transT * work1) * I0) * T;
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        Mbb[i, j] = work2[i, j] + work3[i, j];
                    }
                }
            }
            IvyFEM.Lapack.DoubleMatrix Kk1;
            {
                var work = accPgM * transz;
                work = IvyFEM.Lapack.DoubleMatrix.Scal(work, (1.0 / ln));
                Kk1 = Mb * work;
            }
            var Kk2 = new IvyFEM.Lapack.DoubleMatrix(6, 6);
            {
                double d1 = IvyFEM.Lapack.Functions.ddot(z.Buffer, velPg);
                d1 *= (1.0 / ln);
                var work1 = (Mbb * velPgM) * transz;
                work1 = IvyFEM.Lapack.DoubleMatrix.Scal(work1, d1 * (1.0 / ln));

                var work21 = (Mb * pgM) * transvelPgM;
                var work221 = r * transz;
                var work222 = z * transr;
                var work22 = new IvyFEM.Lapack.DoubleMatrix(6, 6);
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        work22[i, j] = (1.0 / (ln * ln)) * (work221[i, j] + work222[i, j]);
                    }
                }
                var work2 = work21 * work22;

                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        Kk2[i, j] = work1[i, j] - work2[i, j];
                    }
                }
            }
            var Kk3 = new IvyFEM.Lapack.DoubleMatrix(6, 6);
            {
                double d1 = IvyFEM.Lapack.Functions.ddot((transvelPgM * Mbb).Buffer, velPg);
                var work1 = z * transz;
                work1 = IvyFEM.Lapack.DoubleMatrix.Scal(work1, d1 * (1.0 / (ln * ln)));

                double d2 = IvyFEM.Lapack.Functions.ddot((transvelPgM * Mb).Buffer, velPg);
                var work221 = r * transz;
                var work222 = z * transr;
                var work2 = new IvyFEM.Lapack.DoubleMatrix(6, 6);
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        work2[i, j] = d2 * (1.0 / (ln * ln)) * (work221[i, j] + work222[i, j]);
                    }
                }

                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        Kk3[i, j] = (1.0 / 2.0) * (work1[i, j] - work2[i, j]);
                    }
                }
            }
            var Kk = new IvyFEM.Lapack.DoubleMatrix(6, 6);
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    Kk[i, j] = Kk1[i, j] + Kk2[i, j] - Kk3[i, j];
                }
            }
            //-----------------------

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

                    double kValue = K[row * localDof, col * localDof];
                    double mValue = M[row * localDof, col * localDof];
                    double ckValue = Ck[row * localDof, col * localDof];
                    double kkValue = Kk[row * localDof, col * localDof];
                    double curU = pg[col * localDof];
                    double curVel = velPg[col * localDof];
                    double curAcc = accPg[col * localDof];
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

                    double kValue = K[row * localDof, col * localDof + localD2Offset];
                    double mValue = M[row * localDof, col * localDof + localD2Offset];
                    double ckValue = Ck[row * localDof, col * localDof + localD2Offset];
                    double kkValue = Kk[row * localDof, col * localDof + localD2Offset];
                    double curU = pg[col * localDof + localD2Offset];
                    double curVel = velPg[col * localDof + localD2Offset];
                    double curAcc = accPg[col * localDof + localD2Offset];
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

                    double kValue = K[row * localDof, col * localDof + localROffset];
                    double mValue = M[row * localDof, col * localDof + localROffset];
                    double ckValue = Ck[row * localDof, col * localDof + localROffset];
                    double kkValue = Kk[row * localDof, col * localDof + localROffset];
                    double curU = pg[col * localDof + localROffset];
                    double curVel = velPg[col * localDof + localROffset];
                    double curAcc = accPg[col * localDof + localROffset];
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

                    double kValue = K[row * localDof + localD2Offset, col * localDof];
                    double mValue = M[row * localDof + localD2Offset, col * localDof];
                    double ckValue = Ck[row * localDof + localD2Offset, col * localDof];
                    double kkValue = Kk[row * localDof + localD2Offset, col * localDof];
                    double curU = pg[col * localDof];
                    double curVel = velPg[col * localDof];
                    double curAcc = accPg[col * localDof];
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

                    double kValue = K[row * localDof + localD2Offset, col * localDof + localD2Offset];
                    double mValue = M[row * localDof + localD2Offset, col * localDof + localD2Offset];
                    double ckValue = Ck[row * localDof + localD2Offset, col * localDof + localD2Offset];
                    double kkValue = Kk[row * localDof + localD2Offset, col * localDof + localD2Offset];
                    double curU = pg[col * localDof + localD2Offset];
                    double curVel = velPg[col * localDof + localD2Offset];
                    double curAcc = accPg[col * localDof + localD2Offset];
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

                    double kValue = K[row * localDof + localD2Offset, col * localDof + localROffset];
                    double mValue = M[row * localDof + localD2Offset, col * localDof + localROffset];
                    double ckValue = Ck[row * localDof + localD2Offset, col * localDof + localROffset];
                    double kkValue = Kk[row * localDof + localD2Offset, col * localDof + localROffset];
                    double curU = pg[col * localDof + localROffset];
                    double curVel = velPg[col * localDof + localROffset];
                    double curAcc = accPg[col * localDof + localROffset];
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

                    double kValue = K[row * localDof + localROffset, col * localDof];
                    double mValue = M[row * localDof + localROffset, col * localDof];
                    double ckValue = Ck[row * localDof + localROffset, col * localDof];
                    double kkValue = Kk[row * localDof + localROffset, col * localDof];
                    double curU = pg[col * localDof];
                    double curVel = velPg[col * localDof];
                    double curAcc = accPg[col * localDof];
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

                    double kValue = K[row * localDof + localROffset, col * localDof + localD2Offset];
                    double mValue = M[row * localDof + localROffset, col * localDof + localD2Offset];
                    double ckValue = Ck[row * localDof + localROffset, col * localDof + localD2Offset];
                    double kkValue = Kk[row * localDof + localROffset, col * localDof + localD2Offset];
                    double curU = pg[col * localDof + localD2Offset];
                    double curVel = velPg[col * localDof + localD2Offset];
                    double curAcc = accPg[col * localDof + localD2Offset];
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

                    double kValue = K[row * localDof + localROffset, col * localDof + localROffset];
                    double mValue = M[row * localDof + localROffset, col * localDof + localROffset];
                    double ckValue = Ck[row * localDof + localROffset, col * localDof + localROffset];
                    double kkValue = Kk[row * localDof + localROffset, col * localDof + localROffset];
                    double curU = pg[col * localDof + localROffset];
                    double curVel = velPg[col * localDof + localROffset];
                    double curAcc = accPg[col * localDof + localROffset];
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
