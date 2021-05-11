using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class FluidNSDroplet2DRKTDFEM : Fluid2DTDFEM
    {
        //
        public double[] Gravity { get; set; } = { 0.0, -9.81 }; // 外から指定
        // viscosity
        public double Muf { get; set; } = 0.0; // 外から指定
        // 密度
        public double Rhof { get; set; } = 0.0; // 外から指定
        public double Rhop { get; set; } = 0.0; // 外から指定
        // particle体積、半径、質量
        public double Vp { get; set; } = 0.0;  // 外から指定
        public double Rp { get; set; } = 0.0;  // 外から指定
        public double Mp { get; set; } = 0.0;  // 外から指定

        // droplet位置(in/out)
        public double[] PosXd { get; set; } = new double[2];

        // output
        public double[] particleU { get; private set; } = new double[2];

        public FluidNSDroplet2DRKTDFEM(FEWorld world,
            double timeStep,
            double newmarkBeta, double newmarkGamma,
            uint valueId, uint prevValueId)
            : base(world,
                timeStep,
                newmarkBeta, newmarkGamma,
                valueId, prevValueId)
        {
            //EquationType = FluidEquationType.StdGNavierStokes;
            EquationType = FluidEquationType.SUPGNavierStokes;
        }

        private void SolveDropletEquation()
        {
            double dt = TimeStep;

            uint vQuantityId = 0;
            uint pQuantityId = 1;
            int vDof = 2;
            int pDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int pOffset = vNodeCnt * vDof;

            var FV = World.GetFieldValue(ValueId);

            //----------------------------------------------------
            // 速度
            double[] uf = new double[] { 0.0, 0.0 }; // RANSの結果からとってくる
            double[] dotUf = new double[] { 0.0, 0.0 }; // RANSの結果からとってくる
            double[] up = new double[] { 0.0, 0.0 }; // 前回のparticleの速度を入れる
            double[] gradP = new double[] { 0.0, 0.0 }; // RANSの結果からとってくる
            //
            double[] xd = PosXd;
            {
                uint feId = World.GetTriangleFEWithPointInside(vQuantityId, xd);
                System.Diagnostics.Debug.Assert(feId != 0);
                TriangleFE vTriFE = World.GetTriangleFE(vQuantityId, feId);
                TriangleFE pTriFE = World.GetTriangleFE(pQuantityId, feId);
                int[] vCoIds = vTriFE.NodeCoordIds;
                uint vElemNodeCnt = vTriFE.NodeCount;
                int[] vNodes = new int[vElemNodeCnt];
                for (int iNode = 0; iNode < vElemNodeCnt; iNode++)
                {
                    int coId = vCoIds[iNode];
                    int nodeId = World.Coord2Node(vQuantityId, coId);
                    vNodes[iNode] = nodeId;
                }
                int[] pCoIds = pTriFE.NodeCoordIds;
                uint pElemNodeCnt = pTriFE.NodeCount;
                int[] pNodes = new int[pElemNodeCnt];
                for (int iNode = 0; iNode < pElemNodeCnt; iNode++)
                {
                    int coId = pCoIds[iNode];
                    int nodeId = World.Coord2Node(pQuantityId, coId);
                    pNodes[iNode] = nodeId;
                }

                double[] L = vTriFE.Coord2L(xd);
                double[] vN = vTriFE.CalcN(L);
                double[][] vNu = vTriFE.CalcNu(L);
                double[] vNx = vNu[0];
                double[] vNy = vNu[1];
                double[] pN = pTriFE.CalcN(L);
                double[][] pNu = pTriFE.CalcNu(L);
                double[] pNx = pNu[0];
                double[] pNy = pNu[1];

                double[] v = new double[vDof];
                double[] vx = new double[vDof];
                double[] vy = new double[vDof];
                for (int iNode = 0; iNode < vElemNodeCnt; iNode++)
                {
                    int nodeId = vNodes[iNode];
                    if (nodeId == -1)
                    {
                        continue;
                    }

                    for (int iDof = 0; iDof < vDof; iDof++)
                    {
                        double vValue = U[nodeId * vDof + iDof];
                        v[iDof] += vValue * vN[iNode];
                        vx[iDof] += vValue * vNx[iNode];
                        vy[iDof] += vValue * vNy[iNode];
                    }
                }
                double[][] vu = { vx, vy };

                double press = 0.0;
                double pressX = 0.0;
                double pressY = 0.0;
                for (int iNode = 0; iNode < pElemNodeCnt; iNode++)
                {
                    int nodeId = pNodes[iNode];
                    if (nodeId == -1)
                    {
                        continue;
                    }

                    double pValue = U[pOffset + nodeId];
                    press += pValue * pN[iNode];
                    pressX += pValue * pNx[iNode];
                    pressY += pValue * pNy[iNode];
                }

                double[] velov = new double[vDof];
                for (int iNode = 0; iNode < vElemNodeCnt; iNode++)
                {
                    int coId = vCoIds[iNode];
                    double[] u = FV.GetDoubleValue(coId, FieldDerivativeType.Value);
                    double[] vel = FV.GetDoubleValue(coId, FieldDerivativeType.Velocity);
                    double[] acc = FV.GetDoubleValue(coId, FieldDerivativeType.Acceleration);

                    for (int iDof = 0; iDof < vDof; iDof++)
                    {
                        velov[iDof] += vel[iDof] * vN[iNode];
                    }
                }

                ///////////////////
                uf[0] = v[0];
                uf[1] = v[1];
                dotUf[0] = velov[0];
                dotUf[1] = velov[1];
                gradP[0] = pressX;
                gradP[1] = pressY;
            }

            up[0] = particleU[0];
            up[1] = particleU[1];
            //----------------------------------------------------

            ///////////////////////////////////////////////////////////////
            // RK4
            int rkStepCnt = 4;
            double[][] rkkn = new double[rkStepCnt][];
            for (int iStep = 0; iStep < rkStepCnt; iStep++)
            {
                double[] tmpUp = new double[2];
                up.CopyTo(tmpUp, 0);

                if (iStep == 0)
                {
                    rkkn[0] = CalcDropletRKFunc(dt, tmpUp, uf, dotUf, Muf, Rhof, Rhop, Rp, Gravity, gradP);
                }
                else if (iStep == 1)
                {
                    for (int i = 0; i < 2; i++)
                    {
                        tmpUp[i] += (1.0 / 2.0) * rkkn[0][i];
                    }
                    rkkn[1] = CalcDropletRKFunc(dt, tmpUp, uf, dotUf, Muf, Rhof, Rhop, Rp, Gravity, gradP);
                }
                else if (iStep == 2)
                {
                    for (int i = 0; i < 2; i++)
                    {
                        tmpUp[i] += (1.0 / 2.0) * rkkn[1][i];
                    }
                    rkkn[2] = CalcDropletRKFunc(dt, tmpUp, uf, dotUf, Muf, Rhof, Rhop, Rp, Gravity, gradP);
                }
                else if (iStep == 3)
                {
                    for (int i = 0; i < 2; i++)
                    {
                        tmpUp[i] += rkkn[2][i];
                    }
                    rkkn[3] = CalcDropletRKFunc(dt, tmpUp, uf, dotUf, Muf, Rhof, Rhop, Rp, Gravity, gradP);
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
            }
            double[] newUp = new double[vDof];
            for (int i = 0; i < 2; i++)
            {
                newUp[i] = up[i] + (1.0 / 6.0) * (rkkn[0][i] + 2.0 * rkkn[1][i] + 2.0 * rkkn[2][i] + rkkn[3][i]);
            }

            newUp.CopyTo(particleU, 0);

            xd[0] += particleU[0] * dt;
            xd[1] += particleU[1] * dt;
            xd.CopyTo(PosXd, 0);
        }

        private double[] CalcDropletRKFunc(
            double dt,
            double[] up,
            double[] uf,
            double[] dotUf,
            double muf,
            double rhof, double rhop,
            double Rp,
            double[] g,
            double[] gradP)
        {
            double Rep = (2.0 * Rp *
                Math.Sqrt((uf[0] - up[0]) * (uf[0] - up[0]) + (uf[1] - up[1]) * (uf[1] - up[1])) *
                rhof) / muf;
            if (Rep < 1.0e-20)
            {
                Rep = 1.0e-20;
            }
            double cd = 0.0;
            if (Rep < 1.0)
            {
                cd = 24.0 / Rep;
            }
            else if (Rep >= 1.0 && Rep <= 1000.0)
            {
                cd = (24.0 / Rep) * (1.0 + 0.5 * Math.Pow(Rep, 0.687));
            }
            else
            {
                // Rep > 1000.0
                cd = 0.44;
            }

            double[] ret = new double[] { 0.0, 0.0 };

            double m = 1.0 + rhof / (2.0 * rhop);
            double k1 = (3.0 / 4.0) * (cd / (2.0 * Rp)) * (rhof / rhop) *
                Math.Sqrt((uf[0] - up[0]) * (uf[0] - up[0]) + (uf[1] - up[1]) * (uf[1] - up[1]));
            for (int iDim = 0; iDim < 2; iDim++)
            {
                ret[iDim] += k1 * (uf[iDim] - up[iDim]);
            }
            double k2 = (1.0 - rhof / rhop);
            for (int iDim = 0; iDim < 2; iDim++)
            {
                ret[iDim] += k2 * g[iDim];
            }
            double k3 = rhof / (2.0 * rhop);
            for (int iDim = 0; iDim < 2; iDim++)
            {
                ret[iDim] += k3 * dotUf[iDim];
            }
            double k4 = (1.0 / rhop);
            for (int iDim = 0; iDim < 2; iDim++)
            {
                ret[iDim] += k4 * gradP[iDim];
            }

            for (int iDim = 0; iDim < 2; iDim++)
            {
                ret[iDim] *= dt * (1.0 / m);
            }

            return ret;
        }

        public override void Solve()
        {
            // Fluid2DTD
            base.Solve();

            SolveDropletEquation();
        }

    }
}
