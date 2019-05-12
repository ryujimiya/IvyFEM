using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Fluid2DTDFEM
    {
        private void CalcSUPGVorticityAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint wQuantityId = 0;
            uint pQuantityId = 1;
            int wNodeCnt = (int)World.GetNodeCount(wQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int offset = wNodeCnt;
            int vDof = 2; // 速度

            double dt = TimeStep;
            double beta = NewmarkBeta;
            double gamma = NewmarkGamma;
            var FV = World.GetFieldValue(ValueId);

            IList<uint> feIds = World.GetTriangleFEIds(wQuantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE wTriFE = World.GetTriangleFE(wQuantityId, feId);
                TriangleFE pTriFE = World.GetTriangleFE(pQuantityId, feId);
                uint vertexCnt = wTriFE.VertexCount;
                for (int iVertex = 0; iVertex < vertexCnt; iVertex++)
                {
                    System.Diagnostics.Debug.Assert(wTriFE.VertexCoordIds[iVertex] == pTriFE.VertexCoordIds[iVertex]);
                }

                int[] wCoIds = wTriFE.NodeCoordIds;
                uint wElemNodeCnt = wTriFE.NodeCount;
                int[] wNodes = new int[wElemNodeCnt];
                for (int iNode = 0; iNode < wElemNodeCnt; iNode++)
                {
                    int coId = wCoIds[iNode];
                    int nodeId = World.Coord2Node(wQuantityId, coId);
                    wNodes[iNode] = nodeId;
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

                Material ma0 = World.GetMaterial(wTriFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is NewtonFluidMaterial);
                var ma = ma0 as NewtonFluidMaterial;
                double rho = ma.MassDensity;
                double mu = ma.Mu;
                double nu = mu / rho;
                double[] g = { ma.GravityX, ma.GravityY };

                double[][] velos = new double[pElemNodeCnt][];
                for (int iNode = 0; iNode < pElemNodeCnt; iNode++)
                {
                    int coId = pCoIds[iNode];
                    double[] velo = new double[vDof];
                    for (int iDof = 0; iDof < vDof; iDof++)
                    {
                        velo[iDof] = CoordV[coId * vDof + iDof];
                    }
                    velos[iNode] = velo;
                }

                double taum = 0;
                double tauc = 0;
                {
                    double[] aveVelo = {
                        (velos[0][0] + velos[1][0] + velos[2][0]) / 3.0,
                        (velos[0][1] + velos[1][1] + velos[2][1]) / 3.0
                    };
                    double veloNorm = Math.Sqrt(aveVelo[0] * aveVelo[0] + aveVelo[1] * aveVelo[1]);
                    double[][] Lu = new double[vDof][];
                    {
                        double[] a;
                        double[] b;
                        double[] c;
                        wTriFE.CalcTransMatrix(out a, out b, out c);
                        // Lx
                        Lu[0] = b;
                        // Ly
                        Lu[1] = c;
                    }
                    IvyFEM.Lapack.DoubleMatrix GMat = new IvyFEM.Lapack.DoubleMatrix(vDof, vDof);
                    double[] gVec = new double[vDof];
                    for (int iDof = 0; iDof < vDof; iDof++)
                    {
                        for (int jDof = 0; jDof < vDof; jDof++)
                        {
                            for (int kDof = 0; kDof < vDof; kDof++)
                            {
                                GMat[iDof, jDof] += Lu[iDof][kDof] * Lu[jDof][kDof];
                            }
                        }
                    }
                    for (int iDof = 0; iDof < vDof; iDof++)
                    {
                        for (int kDof = 0; kDof < vDof; kDof++)
                        {
                            gVec[iDof] += Lu[iDof][kDof];
                        }
                    }

                    double sqinvtaum1 = 0;
                    double sqinvtaum2 = 0;
                    {
                        double[] tmpVec = GMat * aveVelo;
                        sqinvtaum2 = IvyFEM.Lapack.Functions.ddot(aveVelo, tmpVec);
                    }
                    double sqinvtaum3 = 0;
                    {
                        IvyFEM.Lapack.DoubleMatrix GMatT = new Lapack.DoubleMatrix(GMat);
                        GMatT.Transpose();
                        double GMatDoubleDot = IvyFEM.Lapack.DoubleMatrix.DoubleDot(GMat, GMatT);
                        sqinvtaum3 = nu * nu * GMatDoubleDot;
                    }
                    double sqinvtaum = sqinvtaum1 + sqinvtaum2 + sqinvtaum3;
                    taum = 1.0 / Math.Sqrt(sqinvtaum);

                    double gDot = IvyFEM.Lapack.Functions.ddot(gVec, gVec);
                    tauc = 1.0 / (taum * gDot);
                }

                IntegrationPoints ip = TriangleFE.GetIntegrationPoints(World.TriIntegrationPointCount);//Point7
                for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                {
                    double[] L = ip.Ls[ipPt];
                    double[] wN = wTriFE.CalcN(L);
                    double[][] wNu = wTriFE.CalcNu(L);
                    double[] wNx = wNu[0];
                    double[] wNy = wNu[1];
                    double[,][] wNuv = wTriFE.CalcNuv(L);
                    double[] pN = pTriFE.CalcN(L);
                    double[][] pNu = pTriFE.CalcNu(L);
                    double[] pNx = pNu[0];
                    double[] pNy = pNu[1];

                    double detJ = wTriFE.GetDetJacobian(L);
                    double weight = ip.Weights[ipPt];
                    double detJWeight = (1.0 / 2.0) * weight * detJ;

                    double w = 0;
                    double wx = 0;
                    double wy = 0;
                    double wxx = 0;
                    double wxy = 0;
                    double wyx = 0;
                    double wyy = 0;
                    double velow = 0;
                    for (int iNode = 0; iNode < wElemNodeCnt; iNode++)
                    {
                        int nodeId = wNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }
                        double wValue = U[nodeId];
                        w += wValue * wN[iNode];
                        wx += wValue * wNx[iNode];
                        wy += wValue * wNy[iNode];
                        wxx += wValue * wNuv[0, 0][iNode];
                        wxy += wValue * wNuv[0, 1][iNode];
                        wyx += wValue * wNuv[1, 0][iNode];
                        wyy += wValue * wNuv[1, 1][iNode];

                        int coId = wCoIds[iNode];
                        double[] prevU = FV.GetDoubleValue(coId, FieldDerivativeType.Value);
                        double[] prevVel = FV.GetDoubleValue(coId, FieldDerivativeType.Velocity);
                        double[] prevAcc = FV.GetDoubleValue(coId, FieldDerivativeType.Acceleration);
                        System.Diagnostics.Debug.Assert(prevU.Length == 1);
                        System.Diagnostics.Debug.Assert(prevVel.Length == 1);
                        System.Diagnostics.Debug.Assert(prevAcc.Length == 1);
                        double vel = 0;
                        vel = (gamma / (beta * dt)) * (wValue - prevU[0]) +
                            (1.0 - gamma / beta) * prevVel[0] +
                            dt * (1.0 - gamma / (2.0 * beta)) * prevAcc[0];
                        velow += vel * wN[iNode];
                    }
                    double p = 0;
                    double px = 0;
                    double py = 0;
                    for (int iNode = 0; iNode < pElemNodeCnt; iNode++)
                    {
                        int nodeId = pNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }
                        double pValue = U[offset + nodeId];
                        p += pValue * pN[iNode];
                        px += pValue * pNx[iNode];
                        py += pValue * pNy[iNode];
                    }
                    // dg/du
                    double[] gx = new double[2];
                    double[] gy = new double[2];
                    for (int iNode = 0; iNode < wElemNodeCnt; iNode++)
                    {
                        int nodeId = wNodes[iNode];
                        if (nodeId == -1)
                        {
                            continue;
                        }
                        for (int iDof = 0; iDof < 2; iDof++)
                        {
                            double gValue = g[iDof];
                            gx[iDof] += gValue * wNx[iNode];
                            gy[iDof] += gValue * wNy[iNode];
                        }
                    }
                    double[] v = new double[2];
                    v[0] = py;
                    v[1] = -px;

                    for (int row = 0; row < wElemNodeCnt; row++)
                    {
                        int rowNodeId = wNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < wElemNodeCnt; col++)
                        {
                            int colNodeId = wNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }
                            int colCoId = wCoIds[col];
                            // ω、dω/dt、d2ω/dt2
                            double[] u = FV.GetDoubleValue(colCoId, FieldDerivativeType.Value);
                            double[] vel = FV.GetDoubleValue(colCoId, FieldDerivativeType.Velocity);
                            double[] acc = FV.GetDoubleValue(colCoId, FieldDerivativeType.Acceleration);
                            System.Diagnostics.Debug.Assert(u.Length == 1);
                            System.Diagnostics.Debug.Assert(vel.Length == 1);
                            System.Diagnostics.Debug.Assert(acc.Length == 1);

                            double kww1 = detJWeight * mu * (wNx[row] * wNx[col] + wNy[row] * wNy[col]);
                            double kww2 = detJWeight * rho * wN[row] * (v[0] * wNx[col] + v[1] * wNy[col]);

                            double m = detJWeight * rho * wN[row] * wN[col];

                            A[rowNodeId, colNodeId] +=
                                kww1 + kww2 +
                                (gamma / (beta * dt)) * m;
                            // v = f(ψ) : kww2 Newton-Raphson
                            B[rowNodeId] += kww2 * U[colNodeId] +
                                m * (
                                    (gamma / (beta * dt)) * u[0] -
                                    (1.0 - gamma / beta) * vel[0] -
                                    dt * (1.0 - gamma / (2.0 * beta)) * acc[0]
                                );
                        }
                    }

                    // v = f(ψ): kwp Newton-Raphson
                    for (int row = 0; row < wElemNodeCnt; row++)
                    {
                        int rowNodeId = wNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < pElemNodeCnt; col++)
                        {
                            int colNodeId = pNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double kwp = detJWeight * rho * wN[row] * (pNy[col] * wx - pNx[col] * wy);
                            A[rowNodeId, offset + colNodeId] += kwp;
                            B[rowNodeId] += kwp * U[offset + colNodeId];
                        }
                    }

                    for (int row = 0; row < pElemNodeCnt; row++)
                    {
                        int rowNodeId = pNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < wElemNodeCnt; col++)
                        {
                            int colNodeId = wNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double kpw = -detJWeight * pN[row] * wN[col];
                            A[offset + rowNodeId, colNodeId] += kpw;
                        }
                    }

                    for (int row = 0; row < pElemNodeCnt; row++)
                    {
                        int rowNodeId = pNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < pElemNodeCnt; col++)
                        {
                            int colNodeId = pNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double kpp = detJWeight * (pNx[row] * pNx[col] + pNy[row] * pNy[col]);
                            A[offset + rowNodeId, offset + colNodeId] += kpp;
                        }
                    }

                    // v = f(ψ): qw Newton-Raphson
                    for (int row = 0; row < wElemNodeCnt; row++)
                    {
                        int rowNodeId = wNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        double qw = detJWeight * rho * wN[row] * (v[0] * wx + v[1] * wy);
                        B[rowNodeId] += -qw;
                    }

                    for (int row = 0; row < wElemNodeCnt; row++)
                    {
                        int rowNodeId = wNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        double f = 0;
                        for (int kNode = 0; kNode < wElemNodeCnt; kNode++)
                        {
                            int kNodeId = wNodes[kNode];
                            if (kNodeId == -1)
                            {
                                continue;
                            }
                            f += detJWeight * rho * (wNx[kNode] * g[1] - wNy[kNode] * g[0]);
                        }
                        B[rowNodeId] += f;
                    }

                    //////////////////////////////////////////////////////////////
                    // SUPG
                    double rm = 0;
                    double[] rmw = new double[wElemNodeCnt];
                    double[] rmp = new double[pElemNodeCnt];
                    {
                        rm =
                            rho * velow +
                            -mu * (wxx + wyy) +
                            rho * (v[0] * wx + v[1] * wy) +
                            -rho * (gx[1] - gy[0]);

                        for (int jNode = 0; jNode < wElemNodeCnt; jNode++)
                        {
                            int jNodeId = wNodes[jNode];
                            if (jNodeId == -1)
                            {
                                continue;
                            }

                            rmw[jNode] =
                                rho * (gamma / (beta * dt)) * wN[jNode] +
                                -mu * (wNuv[0, 0][jNode] + wNuv[1, 1][jNode]) +
                                rho * (v[0] * wNu[0][jNode] + v[1] * wNu[1][jNode]);
                        }

                        for (int jNode = 0; jNode < pElemNodeCnt; jNode++)
                        {
                            int jNodeId = pNodes[jNode];
                            if (jNodeId == -1)
                            {
                                continue;
                            }
                            rmp[jNode] = rho * (pNu[1][jNode] * wx - pNu[0][jNode] * wy);
                        }
                    }
                    // kww
                    for (int row = 0; row < wElemNodeCnt; row++)
                    {
                        int rowNodeId = wNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < wElemNodeCnt; col++)
                        {
                            int colNodeId = wNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double kww =
                                detJWeight * taum * (v[0] * wNu[0][row] + v[1] * wNu[1][row]) * rmw[col];

                            A[rowNodeId, colNodeId] += kww;

                            B[rowNodeId] +=
                                kww * U[colNodeId];
                        }
                    }
                    // kwp
                    for (int row = 0; row < wElemNodeCnt; row++)
                    {
                        int rowNodeId = wNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < pElemNodeCnt; col++)
                        {
                            int colNodeId = pNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double kwp =
                                detJWeight * taum * (pNu[1][col] * wNu[0][row] - pNu[0][col] * wNu[1][row]) * rm +
                                detJWeight * taum * (v[0] * wNu[0][row] + v[1] * wNu[1][row]) * rmp[col];

                            A[rowNodeId, offset + colNodeId] += kwp;

                            B[rowNodeId] +=
                                kwp * U[offset + colNodeId];
                        }
                    }

                    for (int row = 0; row < wElemNodeCnt; row++)
                    {
                        int rowNodeId = wNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }

                        double qw = detJWeight * taum *
                            (v[0] * wNu[0][row] + v[1] * wNu[1][row]) * rm;

                        B[rowNodeId] += -qw;
                    }
                }

            }
        }
    }
}
