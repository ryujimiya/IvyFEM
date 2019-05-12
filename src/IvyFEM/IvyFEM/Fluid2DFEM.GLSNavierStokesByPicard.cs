using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Fluid2DFEM
    {
        private void CalcGLSNavierStokesByPicardAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint vQuantityId = 0;
            uint pQuantityId = 1;
            int vDof = 2;
            int pDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int offset = vNodeCnt * vDof;

            IList<uint> feIds = World.GetTriangleFEIds(vQuantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE vTriFE = World.GetTriangleFE(vQuantityId, feId);
                TriangleFE pTriFE = World.GetTriangleFE(pQuantityId, feId);
                uint vertexCnt = vTriFE.VertexCount;
                for (int iVertex = 0; iVertex < vertexCnt; iVertex++)
                {
                    System.Diagnostics.Debug.Assert(vTriFE.VertexCoordIds[iVertex] == pTriFE.VertexCoordIds[iVertex]);
                }

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

                Material ma0 = World.GetMaterial(vTriFE.MaterialId);
                System.Diagnostics.Debug.Assert(ma0 is NewtonFluidMaterial);
                var ma = ma0 as NewtonFluidMaterial;
                double rho = ma.MassDensity;
                double mu = ma.Mu;
                double nu = mu / rho;
                double[] g = { ma.GravityX, ma.GravityY };

                double[][] velos = new double[vElemNodeCnt][];
                for (int iNode = 0; iNode < vElemNodeCnt; iNode++)
                {
                    double[] velo = new double[vDof];
                    int nodeId = vNodes[iNode];
                    if (nodeId == -1)
                    {
                        // 0
                    }
                    else
                    {
                        for (int iDof = 0; iDof < vDof; iDof++)
                        {
                            velo[iDof] = U[nodeId * vDof + iDof];
                        }
                    }
                    velos[iNode] = velo;
                }

                /*
                double taum = 0;
                double tauc = 0;
                {
                    double[] aveVelo = {
                        (velos[0][0] + velos[1][0] + velos[2][0]) / 3.0,
                        (velos[0][1] + velos[1][1] + velos[2][1]) / 3.0
                    };
                    double veloNorm = Math.Sqrt(aveVelo[0] * aveVelo[0] + aveVelo[1] * aveVelo[1]);
                    double Ae = vTriFE.GetArea();
                    double h = 2.0 * Math.Sqrt(Ae / Math.PI);
                    double sqinvtaum1 = 0;
                    double sqinvtaum2 = (2.0 * 2.0 * veloNorm * veloNorm) / (h * h);
                    double sqinvtaum3 = (4.0 * 4.0 * nu * nu) / (h * h * h * h);
                    double sqinvtaum = sqinvtaum1 + sqinvtaum2 + sqinvtaum3;
                    taum = 1.0 / Math.Sqrt(sqinvtaum);

                    double re = veloNorm * h / (2.0 * nu);
                    if (re < 3.0)
                    {
                        tauc = (1.0 / 2.0) * h * veloNorm * re / 3.0;
                    }
                    else
                    {
                        tauc = (1.0 / 2.0) * h * veloNorm;
                    }
                }
                */
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
                        vTriFE.CalcTransMatrix(out a, out b, out c);
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
                        sqinvtaum3 = 30.0 * nu * nu * GMatDoubleDot;
                    }
                    double sqinvtaum = sqinvtaum1 + sqinvtaum2 + sqinvtaum3;
                    taum = 1.0 / Math.Sqrt(sqinvtaum);

                    double gDot = IvyFEM.Lapack.Functions.ddot(gVec, gVec);
                    tauc = 1.0 / (taum * gDot);
                }

                double[] vSN = vTriFE.CalcSN();
                IntegrationPoints ip = TriangleFE.GetIntegrationPoints(World.TriIntegrationPointCount);//Point7
                for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                {
                    double[] L = ip.Ls[ipPt];
                    double[] vN = vTriFE.CalcN(L);
                    double[][] vNu = vTriFE.CalcNu(L);
                    double[] vNx = vNu[0];
                    double[] vNy = vNu[1];
                    double[,][] vNuv = vTriFE.CalcNuv(L);
                    double[] pN = pTriFE.CalcN(L);
                    double[][] pNu = pTriFE.CalcNu(L);
                    double[] pNx = pNu[0];
                    double[] pNy = pNu[1];

                    double detJ = vTriFE.GetDetJacobian(L);
                    double weight = ip.Weights[ipPt];
                    double detJWeight = (1.0 / 2.0) * weight * detJ;

                    double[] v = new double[vDof];
                    double[] vx = new double[vDof];
                    double[] vy = new double[vDof];
                    double[] vxx = new double[vDof];
                    double[] vxy = new double[vDof];
                    double[] vyx = new double[vDof];
                    double[] vyy = new double[vDof];
                    double p = 0;
                    double px = 0;
                    double py = 0;
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
                            vxx[iDof] += vValue * vNuv[0, 0][iNode];
                            vxy[iDof] += vValue * vNuv[0, 1][iNode];
                            vyx[iDof] += vValue * vNuv[1, 0][iNode];
                            vyy[iDof] += vValue * vNuv[1, 1][iNode];
                        }
                    }
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
                    double[][] vu = new double[vDof][];
                    vu[0] = vx;
                    vu[1] = vy;
                    double[,][] vuv = new double[vDof, vDof][];
                    vuv[0, 0] = vxx;
                    vuv[0, 1] = vxy;
                    vuv[1, 0] = vyx;
                    vuv[1, 1] = vyy;
                    double[] pu = new double[vDof];
                    pu[0] = px;
                    pu[1] = py;

                    for (int row = 0; row < vElemNodeCnt; row++)
                    {
                        int rowNodeId = vNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < vElemNodeCnt; col++)
                        {
                            int colNodeId = vNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double[,] kvv1 = new double[vDof, vDof];
                            kvv1[0, 0] = detJWeight * mu * (vNx[row] * vNx[col] +
                                vNx[row] * vNx[col] + vNy[row] * vNy[col]);
                            kvv1[0, 1] = detJWeight * mu * vNy[row] * vNx[col];
                            kvv1[1, 0] = detJWeight * mu * vNx[row] * vNy[col];
                            kvv1[1, 1] = detJWeight * mu * (vNy[row] * vNy[col] +
                                vNx[row] * vNx[col] + vNy[row] * vNy[col]);

                            double[,] kvv2 = new double[vDof, vDof];
                            //kvv2[0, 0] = detJWeight * rho * vN[row] * (
                            //    vN[col] * vx[0] + v[0] * vNx[col] + v[1] * vNy[col]);
                            //kvv2[0, 1] = detJWeight * rho * vN[row] * vN[col] * vy[0];
                            //kvv2[1, 0] = detJWeight * rho * vN[row] * vN[col] * vx[1];
                            //kvv2[1, 1] = detJWeight * rho * vN[row] * (
                            //    vN[col] * vy[1] + v[0] * vNx[col] + v[1] * vNy[col]);
                            // Picard
                            kvv2[0, 0] = detJWeight * rho * vN[row] * (v[0] * vNx[col] + v[1] * vNy[col]);
                            kvv2[0, 1] = 0;
                            kvv2[1, 0] = 0;
                            kvv2[1, 1] = detJWeight * rho * vN[row] * (v[0] * vNx[col] + v[1] * vNy[col]);

                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                for (int colDof = 0; colDof < vDof; colDof++)
                                {
                                    A[rowNodeId * vDof + rowDof, colNodeId * vDof + colDof] +=
                                        kvv1[rowDof, colDof] + kvv2[rowDof, colDof];

                                    //B[rowNodeId * vDof + rowDof] +=
                                    //    kvv2[rowDof, colDof] * U[colNodeId * vDof + colDof];
                                    // Picard
                                    // nothing
                                }
                            }
                        }
                    }

                    for (int row = 0; row < vElemNodeCnt; row++)
                    {
                        int rowNodeId = vNodes[row];
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

                            double[,] kvp = new double[vDof, pDof];
                            kvp[0, 0] = -detJWeight * vNx[row] * pN[col];
                            kvp[1, 0] = -detJWeight * vNy[row] * pN[col];

                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                A[rowNodeId * vDof + rowDof, offset + colNodeId] += kvp[rowDof, 0];
                                A[offset + colNodeId, rowNodeId * vDof + rowDof] += kvp[rowDof, 0];
                            }
                        }
                    }

                    for (int row = 0; row < vElemNodeCnt; row++)
                    {
                        int rowNodeId = vNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        double[] q2 = new double[vDof];
                        for (int rowDof = 0; rowDof < vDof; rowDof++)
                        {
                            //q2[rowDof] = detJWeight * rho * vN[row] * (v[0] * vx[rowDof] + v[1] * vy[rowDof]);
                            // Picard
                            // nothing
                        }
                        for (int rowDof = 0; rowDof < vDof; rowDof++)
                        {
                            B[rowNodeId * vDof + rowDof] += -q2[rowDof];
                        }
                    }

                    //////////////////////////////////////////////////////////////
                    // SUPG
                    double[] rmi = new double[vDof];
                    double[,,] rmivj = new double[vDof, vDof, vElemNodeCnt];
                    double[,] rmip = new double[vDof, pElemNodeCnt];
                    double rc = 0;
                    double[,] rcvj = new double[vDof, vElemNodeCnt];
                    for (int iDof = 0; iDof < vDof; iDof++)
                    {
                        rmi[iDof] =
                            -mu * (vuv[0, 0][iDof] + vuv[1, 1][iDof]) +
                            rho * (v[0] * vx[iDof] + v[1] * vy[iDof]) +
                            pu[iDof] - rho * g[iDof];

                        for (int jDof = 0; jDof < vDof; jDof++)
                        {
                            for (int jNode = 0; jNode < vElemNodeCnt; jNode++)
                            {
                                int jNodeId = vNodes[jNode];
                                if (jNodeId == -1)
                                {
                                    continue;
                                }

                                rmivj[iDof, jDof, jNode] = 0;
                                if (iDof == jDof)
                                {
                                    rmivj[iDof, jDof, jNode] +=
                                        -mu * (vNuv[0, 0][jNode] + vNuv[1, 1][jNode]);
                                }
                                //rmivj[iDof, jDof, jNode] +=
                                //    rho * vN[jNode] * vu[jDof][iDof];
                                // Picard
                                // nothing
                                if (iDof == jDof)
                                {
                                    rmivj[iDof, jDof, jNode] +=
                                        rho * (v[0] * vNu[0][jNode] + v[1] * vNu[1][jNode]);
                                }
                            }
                        }

                        for (int jNode = 0; jNode < pElemNodeCnt; jNode++)
                        {
                            int jNodeId = pNodes[jNode];
                            if (jNodeId == -1)
                            {
                                continue;
                            }
                            rmip[iDof, jNode] = pNu[iDof][jNode];
                        }
                    }
                    {
                        rc = vx[0] + vy[1];
                        for (int jDof = 0; jDof < vDof; jDof++)
                        {
                            for (int jNode = 0; jNode < vElemNodeCnt; jNode++)
                            {
                                int jNodeId = vNodes[jNode];
                                if (jNodeId == -1)
                                {
                                    continue;
                                }
                                rcvj[jDof, jNode] = vNu[jDof][jNode];
                            }
                        }
                    }
                    // kvv
                    for (int row = 0; row < vElemNodeCnt; row++)
                    {
                        int rowNodeId = vNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < vElemNodeCnt; col++)
                        {
                            int colNodeId = vNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double[,] kvv1 = new double[vDof, vDof];
                            double[,] kvv1GLSAdd = new double[vDof, vDof];
                            double[,] kvv2 = new double[vDof, vDof];
                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                for (int colDof = 0; colDof < vDof; colDof++)
                                {
                                    // Picard
                                    kvv1[rowDof, colDof] =
                                        detJWeight * taum * (v[0] * vNu[0][row] + v[1] * vNu[1][row]) *
                                        rmivj[rowDof, colDof, col];
                                    kvv2[rowDof, colDof] =
                                        detJWeight * tauc * rho * vNu[rowDof][row] * rcvj[colDof, col];
                                }
                            }
                            //  GLS追加項
                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                for (int colDof = 0; colDof < vDof; colDof++)
                                {
                                    kvv1GLSAdd[rowDof, colDof] =
                                        -detJWeight * (1.0 / rho) * taum * mu * (vNuv[0, 0][row] + vNuv[1, 1][row]) *
                                        rmivj[rowDof, colDof, col];
                                }
                            }

                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                for (int colDof = 0; colDof < vDof; colDof++)
                                {
                                    A[rowNodeId * vDof + rowDof, colNodeId * vDof + colDof] +=
                                        kvv1[rowDof, colDof] + kvv2[rowDof, colDof] + kvv1GLSAdd[rowDof, colDof];

                                    // Picard
                                    // nothing
                                }
                            }
                        }
                    }
                    // kvp
                    for (int row = 0; row < vElemNodeCnt; row++)
                    {
                        int rowNodeId = vNodes[row];
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

                            double[,] kvp = new double[vDof, pDof];
                            double[,] kvpGLSAdd = new double[vDof, pDof];
                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                kvp[rowDof, 0] =
                                    detJWeight * taum *
                                    (v[0] * vNu[0][row] + v[1] * vNu[1][row]) * rmip[rowDof, col];
                            }
                            // GLS追加項
                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                kvpGLSAdd[rowDof, 0] =
                                    -detJWeight * (1.0 / rho) * taum * mu *
                                    (vNuv[0, 0][row] + vNuv[1, 1][row]) * rmip[rowDof, col];
                            }

                            for (int rowDof = 0; rowDof < vDof; rowDof++)
                            {
                                A[rowNodeId * vDof + rowDof, offset + colNodeId] += 
                                    kvp[rowDof, 0] + kvpGLSAdd[rowDof, 0];

                                // Picard
                                // nothing
                            }
                        }
                    }
                    // kpv
                    for (int row = 0; row < pElemNodeCnt; row++)
                    {
                        int rowNodeId = pNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        for (int col = 0; col < vElemNodeCnt; col++)
                        {
                            int colNodeId = vNodes[col];
                            if (colNodeId == -1)
                            {
                                continue;
                            }

                            double[,] kpv = new double[pDof, vDof];
                            for (int colDof = 0; colDof < vDof; colDof++)
                            {
                                kpv[0, colDof] =
                                    -detJWeight * (1.0 / rho) * taum *
                                    (pNu[0][row] * rmivj[0, colDof, col] + pNu[1][row] * rmivj[1, colDof, col]);
                            }

                            for (int colDof = 0; colDof < vDof; colDof++)
                            {
                                A[offset + rowNodeId, colNodeId * vDof + colDof] += kpv[0, colDof];

                                // Picard
                                // nothing
                            }
                        }
                    }
                    // kpp
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

                            double[,] kpp = new double[pDof, pDof];
                            kpp[0, 0] =
                                -detJWeight * (1.0 / rho) * taum *
                                (pNu[0][row] * rmip[0, col] + pNu[1][row] * rmip[1, col]);

                            A[offset + rowNodeId, offset + colNodeId] += kpp[0, 0];

                            // Picard
                            // nothing
                        }
                    }

                    for (int row = 0; row < vElemNodeCnt; row++)
                    {
                        int rowNodeId = vNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        double[] qv1 = new double[vDof];
                        double[] qv2 = new double[vDof];
                        for (int rowDof = 0; rowDof < vDof; rowDof++)
                        {
                            // Picard
                            // nothing
                        }
                        for (int rowDof = 0; rowDof < vDof; rowDof++)
                        {
                            B[rowNodeId * vDof + rowDof] += -(qv1[rowDof] + qv2[rowDof]);
                        }
                    }

                    for (int row = 0; row < pElemNodeCnt; row++)
                    {
                        int rowNodeId = pNodes[row];
                        if (rowNodeId == -1)
                        {
                            continue;
                        }
                        double[] qp = new double[pDof];
                        // Picard
                        // nothing
                        B[offset + rowNodeId] += -qp[0];
                    }
                }

                for (int row = 0; row < vElemNodeCnt; row++)
                {
                    int rowNodeId = vNodes[row];
                    if (rowNodeId == -1)
                    {
                        continue;
                    }
                    for (int rowDof = 0; rowDof < vDof; rowDof++)
                    {
                        B[rowNodeId * vDof + rowDof] += rho * g[rowDof] * vSN[row];
                    }
                }
            }
        }
    }
}
