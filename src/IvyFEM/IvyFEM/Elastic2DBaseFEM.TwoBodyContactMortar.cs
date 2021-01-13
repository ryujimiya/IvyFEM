using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class Elastic2DBaseFEM
    {
        private void SetTwoBodyContactMortarSegmentationQuantitySpecialBC(
            uint cQuantityId, IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint uQuantityId = 0;
            System.Diagnostics.Debug.Assert(World.GetFEOrder(uQuantityId) == 1);
            System.Diagnostics.Debug.Assert(World.GetCoordCount(uQuantityId) ==
                World.GetCoordCount(cQuantityId));
            System.Diagnostics.Debug.Assert(World.GetDof(uQuantityId) == 2);
            System.Diagnostics.Debug.Assert(World.GetDof(cQuantityId) == 2);
            int uDof = 2;
            int cDof = 2;
            int uNodeCnt = (int)World.GetNodeCount(uQuantityId);
            int cNodeCnt = (int)World.GetNodeCount(cQuantityId);
            int offset = World.GetOffset(cQuantityId);

            // 線要素の変位を更新
            UpdateLineFEDisplacements(uQuantityId, uDof, cQuantityId);

            // 節点法線ベクトルの計算
            Dictionary<int, double[]> co2Normal = GetSlaveLineFECo2Normal(uQuantityId, uDof, cQuantityId);

            // MasterからSlaveに垂らした点を計算する
            Dictionary<uint, IList<double>> slaveFEL2s = GetSlavePointFromMasterNodes(
                co2Normal,
                uQuantityId, uDof, cQuantityId);

            bool[] lConstraintNodeIds = new bool[cNodeCnt];
            IList<uint> slaveFEIds = World.GetContactSlaveLineFEIds(cQuantityId);
            foreach (uint slaveFEId in slaveFEIds)
            {
                LineFE lineFE = World.GetLineFE(uQuantityId, slaveFEId);
                uint elemNodeCnt = lineFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = lineFE.NodeCoordIds[iNode];
                    int nodeId = World.Coord2Node(uQuantityId, coId);
                    nodes[iNode] = nodeId;
                }

                LineFE lLineFE = World.GetLineFE(cQuantityId, slaveFEId);
                int[] lNodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = lLineFE.NodeCoordIds[iNode];
                    int lNodeId = World.Coord2Node(cQuantityId, coId);
                    lNodes[iNode] = lNodeId;
                }

                double[][] curNodeCoord = new double[elemNodeCnt][];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = lineFE.NodeCoordIds[iNode];
                    double[] coord = World.GetCoord(uQuantityId, coId);
                    int iNodeId = nodes[iNode];
                    curNodeCoord[iNode] = new double[uDof];
                    if (iNodeId == -1)
                    {
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            curNodeCoord[iNode][iDof] = coord[iDof];
                        }
                    }
                    else
                    {
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            curNodeCoord[iNode][iDof] = coord[iDof] + U[iNodeId * uDof + iDof];
                        }
                    }
                }

                //IntegrationPoints ip = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point3);
                //System.Diagnostics.Debug.Assert(ip.Ls.Length == 3);
                IntegrationPoints ip = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);
                System.Diagnostics.Debug.Assert(ip.Ls.Length == 5);

                IList<double> L2s = new List<double>();
                L2s.Add(0.0);
                if (slaveFEL2s.ContainsKey(slaveFEId))
                {
                    foreach (double L2 in slaveFEL2s[slaveFEId])
                    {
                        L2s.Add(L2);
                    }
                }
                L2s.Add(1.0);
                int segCnt = L2s.Count - 1;
                for (int iSeg = 0; iSeg < segCnt; iSeg++)
                {
                    double sL2 = L2s[iSeg];
                    double eL2 = L2s[iSeg + 1];
                    OpenTK.Vector2d sPt = new OpenTK.Vector2d(
                        (1.0 - sL2) * curNodeCoord[0][0] + sL2 * curNodeCoord[1][0],
                        (1.0 - sL2) * curNodeCoord[0][1] + sL2 * curNodeCoord[1][1]
                    );
                    OpenTK.Vector2d ePt = new OpenTK.Vector2d(
                        (1.0 - eL2) * curNodeCoord[0][0] + eL2 * curNodeCoord[1][0],
                        (1.0 - eL2) * curNodeCoord[0][1] + eL2 * curNodeCoord[1][1]
                    );
                    double segLen = (ePt - sPt).Length;
                    for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                    {
                        double[] segL = ip.Ls[ipPt];
                        double weight = ip.Weights[ipPt];
                        double L2 = segL[0] * sL2 + segL[1] * eL2;
                        double[] L = { 1.0 - L2, L2 };
                        double[] N = lineFE.CalcN(L);
                        double[] lN = lLineFE.CalcN(L);
                        double detJWeight = (segLen / 2.0) * weight;
                        // 現在の位置
                        double[] curCoord = new double[uDof];
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            for (int iDof = 0; iDof < uDof; iDof++)
                            {
                                curCoord[iDof] += curNodeCoord[iNode][iDof] * N[iNode];
                            }
                        }
                        // 連続な近似法線ベクトルを計算する
                        double[] normal = new double[uDof];
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            int coId = lineFE.NodeCoordIds[iNode];
                            double[] nodeNormal = co2Normal[coId];
                            for (int iDof = 0; iDof < uDof; iDof++)
                            {
                                normal[iDof] += nodeNormal[iDof] * N[iNode];
                            }
                        }
                        normal = IvyFEM.Lapack.Utils.NormalizeDoubleVector(normal);
                        double[] tan = { normal[1], -normal[0] };

                        // 対応するMasterの点を取得する
                        uint masterFEId;
                        double[] masterL;
                        GetMasterLineFEPoint(
                            curCoord, normal,
                            uQuantityId, uDof, cQuantityId,
                            out masterFEId, out masterL);
                        if (masterFEId == 0)
                        {
                            // 対応するMasterの点がない
                            continue;
                        }
                        LineFE masterLineFE = World.GetLineFE(uQuantityId, masterFEId);
                        System.Diagnostics.Debug.Assert(masterLineFE.NodeCount == elemNodeCnt);
                        double[] masterN = masterLineFE.CalcN(masterL);
                        int[] masterNodes = new int[elemNodeCnt];
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            int coId = masterLineFE.NodeCoordIds[iNode];
                            int nodeId = World.Coord2Node(uQuantityId, coId);
                            masterNodes[iNode] = nodeId;
                        }
                        // 現在の位置
                        double[] masterCurCoord = new double[uDof];
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            int coId = masterLineFE.NodeCoordIds[iNode];
                            double[] coord = World.GetCoord(uQuantityId, coId);
                            int iNodeId = masterNodes[iNode];
                            if (iNodeId == -1)
                            {
                                for (int iDof = 0; iDof < uDof; iDof++)
                                {
                                    masterCurCoord[iDof] += coord[iDof] * masterN[iNode];
                                }
                            }
                            else
                            {
                                for (int iDof = 0; iDof < uDof; iDof++)
                                {
                                    masterCurCoord[iDof] +=
                                        (coord[iDof] + U[iNodeId * uDof + iDof]) * masterN[iNode];
                                }
                            }
                        }

                        // ギャップの計算
                        double gap = 0;
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            gap += -normal[iDof] * (curCoord[iDof] - masterCurCoord[iDof]);
                        }

                        // ラグランジュの未定乗数
                        double[] l = new double[cDof];
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            int iNodeId = lNodes[iNode];
                            if (iNodeId == -1)
                            {
                                continue;
                            }
                            for (int iDof = 0; iDof < cDof; iDof++)
                            {
                                l[iDof] += U[offset + iNodeId * cDof + iDof] * lN[iNode];
                            }
                        }
                        // Karush-Kuhn-Tucker条件
                        double tolerance = IvyFEM.Linear.Constants.ConvRatioTolerance;
                        if (l[0] <= tolerance &&
                            Math.Abs(l[1]) <= tolerance &&
                            gap >= -tolerance)
                        {
                            // 拘束しない
                            continue;
                        }

                        ////////////////////////////////////////
                        // これ以降、条件を付加する処理
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            int iNodeId = lNodes[iNode];
                            if (iNodeId == -1)
                            {
                                continue;
                            }
                            lConstraintNodeIds[iNodeId] = true;
                        }

                        // Slave
                        for (int row = 0; row < elemNodeCnt; row++)
                        {
                            int rowNodeId = nodes[row];
                            if (rowNodeId == -1)
                            {
                                continue;
                            }
                            for (int col = 0; col < elemNodeCnt; col++)
                            {
                                int colNodeId = lNodes[col];
                                if (colNodeId == -1)
                                {
                                    continue;
                                }

                                double[,] kul = new double[uDof, cDof];
                                double[,] klu = new double[cDof, uDof];
                                for (int rowDof = 0; rowDof < uDof; rowDof++)
                                {
                                    kul[rowDof, 0] +=
                                        detJWeight * normal[rowDof] * N[row] * lN[col];
                                    klu[0, rowDof] +=
                                        detJWeight * normal[rowDof] * N[row] * lN[col];
                                    kul[rowDof, 1] +=
                                        detJWeight * tan[rowDof] * N[row] * lN[col];
                                    klu[1, rowDof] +=
                                        detJWeight * tan[rowDof] * N[row] * lN[col];
                                }

                                for (int rowDof = 0; rowDof < uDof; rowDof++)
                                {
                                    for (int colDof = 0; colDof < cDof; colDof++)
                                    {
                                        A[rowNodeId * uDof + rowDof, offset + colNodeId * cDof + colDof] +=
                                            kul[rowDof, colDof];
                                        A[offset + colNodeId * cDof + colDof, rowNodeId * uDof + rowDof] +=
                                            klu[colDof, rowDof];
                                        B[rowNodeId * uDof + rowDof] +=
                                            kul[rowDof, colDof] * U[offset + colNodeId * cDof + colDof];
                                        B[offset + colNodeId * cDof + colDof] +=
                                            klu[colDof, rowDof] * U[rowNodeId * uDof + rowDof];
                                    }
                                }
                            }
                        }

                        // Master
                        for (int row = 0; row < elemNodeCnt; row++)
                        {
                            int rowNodeId = masterNodes[row];
                            if (rowNodeId == -1)
                            {
                                continue;
                            }
                            for (int col = 0; col < elemNodeCnt; col++)
                            {
                                int colNodeId = lNodes[col];
                                if (colNodeId == -1)
                                {
                                    continue;
                                }

                                double[,] kul = new double[uDof, cDof];
                                double[,] klu = new double[cDof, uDof];
                                for (int rowDof = 0; rowDof < uDof; rowDof++)
                                {
                                    kul[rowDof, 0] +=
                                        -detJWeight * normal[rowDof] * masterN[row] * lN[col];
                                    klu[0, rowDof] +=
                                        -detJWeight * normal[rowDof] * masterN[row] * lN[col];
                                    kul[rowDof, 1] +=
                                        -detJWeight * tan[rowDof] * masterN[row] * lN[col];
                                    klu[1, rowDof] +=
                                        -detJWeight * tan[rowDof] * masterN[row] * lN[col];
                                }

                                for (int rowDof = 0; rowDof < uDof; rowDof++)
                                {
                                    for (int colDof = 0; colDof < cDof; colDof++)
                                    {
                                        A[rowNodeId * uDof + rowDof, offset + colNodeId * cDof + colDof] +=
                                            kul[rowDof, colDof];
                                        A[offset + colNodeId * cDof + colDof, rowNodeId * uDof + rowDof] +=
                                            klu[colDof, rowDof];
                                        B[rowNodeId * uDof + rowDof] +=
                                            kul[rowDof, colDof] * U[offset + colNodeId * cDof + colDof];
                                        B[offset + colNodeId * cDof + colDof] +=
                                            klu[colDof, rowDof] * U[rowNodeId * uDof + rowDof];
                                    }
                                }
                            }
                        }

                        // Slave
                        double[,] qu = new double[elemNodeCnt, uDof];
                        double[,] ql = new double[elemNodeCnt, cDof];
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            for (int iDof = 0; iDof < uDof; iDof++)
                            {
                                qu[iNode, iDof] += 
                                    detJWeight * (l[0] * normal[iDof] + l[1] * tan[iDof]) * N[iNode];
                            }
                        }
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            for (int iDof = 0; iDof < cDof; iDof++)
                            {
                                ql[iNode, 0] +=
                                    detJWeight * lN[iNode] * normal[iDof] * curCoord[iDof];
                                ql[iNode, 1] +=
                                    detJWeight * lN[iNode] * tan[iDof] * curCoord[iDof];
                            }
                        }

                        for (int row = 0; row < elemNodeCnt; row++)
                        {
                            int rowNodeId = nodes[row];
                            if (rowNodeId == -1)
                            {
                                continue;
                            }
                            for (int rowDof = 0; rowDof < uDof; rowDof++)
                            {
                                B[rowNodeId * uDof + rowDof] += -qu[row, rowDof];
                            }
                        }
                        for (int row = 0; row < elemNodeCnt; row++)
                        {
                            int rowNodeId = lNodes[row];
                            if (rowNodeId == -1)
                            {
                                continue;
                            }
                            for (int rowDof = 0; rowDof < cDof; rowDof++)
                            {
                                B[offset + rowNodeId * cDof + rowDof] += -ql[row, rowDof];
                            }
                        }

                        // Master
                        double[,] masterQu = new double[elemNodeCnt, uDof];
                        double[,] masterQl = new double[elemNodeCnt, cDof];
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            for (int iDof = 0; iDof < uDof; iDof++)
                            {
                                masterQu[iNode, iDof] +=
                                    -detJWeight * (l[0] * normal[iDof] + l[1] * tan[iDof]) * masterN[iNode];
                            }
                        }
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            for (int iDof = 0; iDof < uDof; iDof++)
                            {
                                masterQl[iNode, 0] += 
                                    -detJWeight * lN[iNode] * normal[iDof] * masterCurCoord[iDof];
                                masterQl[iNode, 1] +=
                                    -detJWeight * lN[iNode] * tan[iDof] * masterCurCoord[iDof];
                            }
                        }

                        for (int row = 0; row < elemNodeCnt; row++)
                        {
                            int rowNodeId = masterNodes[row];
                            if (rowNodeId == -1)
                            {
                                continue;
                            }
                            for (int rowDof = 0; rowDof < uDof; rowDof++)
                            {
                                B[rowNodeId * uDof + rowDof] += -masterQu[row, rowDof];
                            }
                        }
                        for (int row = 0; row < elemNodeCnt; row++)
                        {
                            int rowNodeId = lNodes[row];
                            if (rowNodeId == -1)
                            {
                                continue;
                            }
                            for (int rowDof = 0; rowDof < cDof; rowDof++)
                            {
                                B[offset + rowNodeId * cDof + rowDof] += -masterQl[row, rowDof];
                            }
                        }
                    }
                }
            }

            // 条件をセットしなかった節点
            for (int iNodeId = 0; iNodeId < cNodeCnt; iNodeId++)
            {
                if (lConstraintNodeIds[iNodeId])
                {
                    continue;
                }
                for (int iDof = 0; iDof < cDof; iDof++)
                {
                    A[offset + iNodeId * cDof + iDof, offset + iNodeId * cDof + iDof] = 1.0;
                    B[offset + iNodeId * cDof + iDof] = 0;
                }
            }
        }
    }
}
