using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class Elastic2DBaseFEM
    {
        protected void CalcTwoBodyContactAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            for (uint quantityId = 0; quantityId < World.GetQuantityCount(); quantityId++)
            {
                int slaveCnt = World.GetContactSlaveEIds(quantityId).Count;
                int masterCnt = World.GetContactMasterEIds(quantityId).Count;
                if (slaveCnt > 0 && masterCnt > 0)
                {
                    //CalcTwoBodyContactSimpleQuantityAB(quantityId, A, B);
                    //CalcTwoBodyContactSegmentationQuantityAB(quantityId, A, B);
                    CalcTwoBodyContactMortarSegmentationQuantityAB(quantityId, A, B);
                }
            }
        }

        private void CalcTwoBodyContactSimpleQuantityAB(
            uint cQuantityId, IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint uQuantityId = 0;
            System.Diagnostics.Debug.Assert(World.GetFEOrder(uQuantityId) == 1);
            System.Diagnostics.Debug.Assert(World.GetCoordCount(uQuantityId) ==
                World.GetCoordCount(cQuantityId));
            System.Diagnostics.Debug.Assert(World.GetDof(uQuantityId) == 2);
            System.Diagnostics.Debug.Assert(World.GetDof(cQuantityId) == 1);
            int uDof = 2;
            int cDof = 1;
            int uNodeCnt = (int)World.GetNodeCount(uQuantityId);
            int cNodeCnt = (int)World.GetNodeCount(cQuantityId);
            int offset = GetOffset(cQuantityId);

            // 線要素の変位を更新
            UpdateLineFEDisplacements(uQuantityId, uDof, cQuantityId);

            // 節点法線ベクトルの計算
            Dictionary<int, double[]> co2Normal = GetSlaveLineFECo2Normal(uQuantityId, uDof, cQuantityId);

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
                //IntegrationPoints ip = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point3);
                //System.Diagnostics.Debug.Assert(ip.Ls.Length == 3);
                IntegrationPoints ip = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);
                System.Diagnostics.Debug.Assert(ip.Ls.Length == 5);

                for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
                {
                    double[] L = ip.Ls[ipPt];
                    double[] N = lineFE.CalcN(L);
                    double[] lN = lLineFE.CalcN(L);
                    double lineLen = lineFE.GetLineLength();
                    //double[] normal = lineFE.GetNormal();
                    double weight = ip.Weights[ipPt];
                    double detJWeight = (lineLen / 2.0) * weight;

                    // 現在の位置
                    double[] curCoord = new double[uDof];
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int coId = lineFE.NodeCoordIds[iNode];
                        double[] coord = World.GetCoord(uQuantityId, coId);
                        int iNodeId = nodes[iNode];
                        if (iNodeId == -1)
                        {
                            for (int iDof = 0; iDof < uDof; iDof++)
                            {
                                curCoord[iDof] += coord[iDof] * N[iNode];
                            }
                        }
                        else
                        {
                            for (int iDof = 0; iDof < uDof; iDof++)
                            {
                                curCoord[iDof] += (coord[iDof] + U[iNodeId * uDof + iDof]) * N[iNode];
                            }
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
                    double l = 0;
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int iNodeId = lNodes[iNode];
                        if (iNodeId == -1)
                        {
                            continue;
                        }
                        l += U[offset + iNodeId] * lN[iNode];
                    }

                    // Karush-Kuhn-Tucker条件
                    double tolerance = IvyFEM.Linear.Constants.ConvRatioTolerance;
                    if (l <= tolerance && gap >= -tolerance)
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
                            }

                            for (int rowDof = 0; rowDof < uDof; rowDof++)
                            {
                                A[rowNodeId * uDof + rowDof, offset + colNodeId] += kul[rowDof, 0];
                                A[offset + colNodeId, rowNodeId * uDof + rowDof] += klu[0, rowDof];
                                B[rowNodeId * uDof + rowDof] +=
                                    kul[rowDof, 0] * U[offset + colNodeId];
                                B[offset + colNodeId] +=
                                    klu[0, rowDof] * U[rowNodeId * uDof + rowDof];
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
                            }

                            for (int rowDof = 0; rowDof < uDof; rowDof++)
                            {
                                A[rowNodeId * uDof + rowDof, offset + colNodeId] += kul[rowDof, 0];
                                A[offset + colNodeId, rowNodeId * uDof + rowDof] += klu[0, rowDof];
                                B[rowNodeId * uDof + rowDof] +=
                                    kul[rowDof, 0] * U[offset + colNodeId];
                                B[offset + colNodeId] +=
                                    klu[0, rowDof] * U[rowNodeId * uDof + rowDof];
                            }
                        }
                    }

                    // Slave
                    double[,] qu = new double[elemNodeCnt, uDof];
                    double[] ql = new double[elemNodeCnt];
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            qu[iNode, iDof] += detJWeight * l * normal[iDof] * N[iNode];
                        }
                    }
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            ql[iNode] += detJWeight * lN[iNode] * normal[iDof] * curCoord[iDof];
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
                        B[offset + rowNodeId] += -ql[row];
                    }

                    // Master
                    double[,] masterQu = new double[elemNodeCnt, uDof];
                    double[] masterQl = new double[elemNodeCnt];
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            masterQu[iNode, iDof] += -detJWeight * l * normal[iDof] * masterN[iNode];
                        }
                    }
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            masterQl[iNode] += -detJWeight * lN[iNode] * normal[iDof] * masterCurCoord[iDof];
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
                        B[offset + rowNodeId] += -masterQl[row];
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
                A[offset + iNodeId, offset + iNodeId] = 1.0;
                B[offset + iNodeId] = 0;
            }
        }

        private void CalcTwoBodyContactSegmentationQuantityAB(
            uint cQuantityId, IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint uQuantityId = 0;
            System.Diagnostics.Debug.Assert(World.GetFEOrder(uQuantityId) == 1);
            System.Diagnostics.Debug.Assert(World.GetCoordCount(uQuantityId) ==
                World.GetCoordCount(cQuantityId));
            System.Diagnostics.Debug.Assert(World.GetDof(uQuantityId) == 2);
            System.Diagnostics.Debug.Assert(World.GetDof(cQuantityId) == 1);
            int uDof = 2;
            int cDof = 1;
            int uNodeCnt = (int)World.GetNodeCount(uQuantityId);
            int cNodeCnt = (int)World.GetNodeCount(cQuantityId);
            int offset = GetOffset(cQuantityId);

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
                        double L2 = segL[0] * sL2  + segL[1] * eL2;
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
                        double l = 0;
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            int iNodeId = lNodes[iNode];
                            if (iNodeId == -1)
                            {
                                continue;
                            }
                            l += U[offset + iNodeId] * lN[iNode];
                        }

                        // Karush-Kuhn-Tucker条件
                        double tolerance = IvyFEM.Linear.Constants.ConvRatioTolerance;
                        if (l <= tolerance && gap >= -tolerance)
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
                                }

                                for (int rowDof = 0; rowDof < uDof; rowDof++)
                                {
                                    A[rowNodeId * uDof + rowDof, offset + colNodeId] += kul[rowDof, 0];
                                    A[offset + colNodeId, rowNodeId * uDof + rowDof] += klu[0, rowDof];
                                    B[rowNodeId * uDof + rowDof] +=
                                        kul[rowDof, 0] * U[offset + colNodeId];
                                    B[offset + colNodeId] +=
                                        klu[0, rowDof] * U[rowNodeId * uDof + rowDof];
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
                                }

                                for (int rowDof = 0; rowDof < uDof; rowDof++)
                                {
                                    A[rowNodeId * uDof + rowDof, offset + colNodeId] += kul[rowDof, 0];
                                    A[offset + colNodeId, rowNodeId * uDof + rowDof] += klu[0, rowDof];
                                    B[rowNodeId * uDof + rowDof] +=
                                        kul[rowDof, 0] * U[offset + colNodeId];
                                    B[offset + colNodeId] +=
                                        klu[0, rowDof] * U[rowNodeId * uDof + rowDof];
                                }
                            }
                        }

                        // Slave
                        double[,] qu = new double[elemNodeCnt, uDof];
                        double[] ql = new double[elemNodeCnt];
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            for (int iDof = 0; iDof < uDof; iDof++)
                            {
                                qu[iNode, iDof] += detJWeight * l * normal[iDof] * N[iNode];
                            }
                        }
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            for (int iDof = 0; iDof < uDof; iDof++)
                            {
                                ql[iNode] += detJWeight * lN[iNode] * normal[iDof] * curCoord[iDof];
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
                            B[offset + rowNodeId] += -ql[row];
                        }

                        // Master
                        double[,] masterQu = new double[elemNodeCnt, uDof];
                        double[] masterQl = new double[elemNodeCnt];
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            for (int iDof = 0; iDof < uDof; iDof++)
                            {
                                masterQu[iNode, iDof] += -detJWeight * l * normal[iDof] * masterN[iNode];
                            }
                        }
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            for (int iDof = 0; iDof < uDof; iDof++)
                            {
                                masterQl[iNode] += -detJWeight * lN[iNode] * normal[iDof] * masterCurCoord[iDof];
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
                            B[offset + rowNodeId] += -masterQl[row];
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
                A[offset + iNodeId, offset + iNodeId] = 1.0;
                B[offset + iNodeId] = 0;
            }
        }

        private void UpdateLineFEDisplacements(uint uQuantityId, int uDof, uint cQuantityId)
        {
            IList<uint> slaveFEIds = World.GetContactSlaveLineFEIds(cQuantityId);
            IList<uint> masterFEIds = World.GetContactMasterLineFEIds(cQuantityId);
            IList<uint>[] feIdss = { slaveFEIds, masterFEIds };

            foreach (IList<uint> feIds in feIdss)
            {
                foreach (uint feId in feIds)
                {
                    LineFE lineFE = World.GetLineFE(uQuantityId, feId);
                    uint elemNodeCnt = lineFE.NodeCount;
                    int[] nodes = new int[elemNodeCnt];
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        int coId = lineFE.NodeCoordIds[iNode];
                        int nodeId = World.Coord2Node(uQuantityId, coId);
                        nodes[iNode] = nodeId;
                    }
                    double[][] displacements = new double[elemNodeCnt][];
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        double[] u = new double[uDof];
                        int nodeId = nodes[iNode];
                        if (nodeId == -1)
                        {
                            for (int iDof = 0; iDof < uDof; iDof++)
                            {
                                u[iDof] = 0;
                            }
                        }
                        else
                        {
                            for (int iDof = 0; iDof < uDof; iDof++)
                            {
                                u[iDof] = U[nodeId * uDof + iDof];
                            }
                        }
                        displacements[iNode] = u;
                    }
                    lineFE.SetDisplacements(displacements);

                    LineFE lLineFE = World.GetLineFE(cQuantityId, feId);
                    lLineFE.SetDisplacements(displacements);
                }
            }
        }

        private Dictionary<int, double[]> GetSlaveLineFECo2Normal(uint uQuantityId, int uDof, uint cQuantityId)
        {
            Dictionary<int, double[]> co2Normal = new Dictionary<int, double[]>();
            IList<uint> slaveFEIds = World.GetContactSlaveLineFEIds(cQuantityId);
            Dictionary<int, IList<double[]>> co2NormalList = new Dictionary<int, IList<double[]>>();
            foreach (uint slaveFEId in slaveFEIds)
            {
                LineFE lineFE = World.GetLineFE(uQuantityId, slaveFEId);
                uint elemNodeCnt = lineFE.NodeCount;
                double lineLen = lineFE.GetLineLength();
                double[] normal = lineFE.GetNormal();
                // 法線ベクトルに重みを付ける
                int dim = normal.Length;
                for (int iDim = 0; iDim < dim; iDim++)
                {
                    normal[iDim] /= lineLen;
                }

                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = lineFE.NodeCoordIds[iNode];
                    if (!co2NormalList.ContainsKey(coId))
                    {
                        co2NormalList[coId] = new List<double[]>();
                    }
                    co2NormalList[coId].Add(normal);
                }
            }
            foreach (var pair in co2NormalList)
            {
                int coId = pair.Key;
                IList<double[]> normalList = pair.Value;
                OpenTK.Vector2d av = new OpenTK.Vector2d();
                foreach (double[] normal in normalList)
                {
                    av.X += normal[0];
                    av.Y += normal[1];
                }
                av = CadUtils.Normalize(av);
                co2Normal[coId] = new double[] { av.X, av.Y };
            }
            return co2Normal;
        }

        private Dictionary<uint, IList<double>> GetSlavePointFromMasterNodes(
            Dictionary<int, double[]> co2Normal,
            uint uQuantityId, int uDof, uint cQuantityId)
        {
            Dictionary<uint, IList<double>> slaveFEL2s = new Dictionary<uint, IList<double>>();
            IList<uint> masterFEIds = World.GetContactMasterLineFEIds(cQuantityId);
            IList<int> masterCoIds = new List<int>();
            foreach (uint masterFEId in masterFEIds)
            {
                LineFE masterLineFE = World.GetLineFE(uQuantityId, masterFEId);
                uint elemNodeCnt = masterLineFE.NodeCount;
                int[] masterNodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = masterLineFE.NodeCoordIds[iNode];
                    int nodeId = World.Coord2Node(uQuantityId, coId);
                    masterNodes[iNode] = nodeId;
                }
                // 現在の頂点の位置
                double[][] masterCurNodeCoords = new double[elemNodeCnt][];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = masterLineFE.NodeCoordIds[iNode];
                    double[] coord = World.GetCoord(uQuantityId, coId);
                    int iNodeId = masterNodes[iNode];
                    masterCurNodeCoords[iNode] = new double[uDof];
                    if (iNodeId == -1)
                    {
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            masterCurNodeCoords[iNode][iDof] = coord[iDof];
                        }
                    }
                    else
                    {
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            masterCurNodeCoords[iNode][iDof] = coord[iDof] + U[iNodeId * uDof + iDof];
                        }
                    }
                }

                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = masterLineFE.NodeCoordIds[iNode];
                    if (masterCoIds.IndexOf(coId) == -1)
                    {
                        masterCoIds.Add(coId);
                    }
                    else
                    {
                        continue;
                    }
                    uint feId;
                    double[] L;
                    GetSlaveLineFEPoint(
                        masterCurNodeCoords[iNode], co2Normal,
                        uQuantityId, uDof, cQuantityId, 
                        out feId, out L);
                    if (feId == 0)
                    {
                        continue;
                    }
                    if (!slaveFEL2s.ContainsKey(feId))
                    {
                        slaveFEL2s[feId] = new List<double>();
                    }
                    double L2 = L[1];
                    slaveFEL2s[feId].Add(L2);
                }
            }

            uint[] feIds = slaveFEL2s.Keys.ToArray();
            foreach (uint feId in feIds)
            {
                List<double> L2s = slaveFEL2s[feId].ToList();
                L2s.Sort();
                slaveFEL2s[feId] = L2s;
            }
            return slaveFEL2s;
        }

        private void GetSlaveLineFEPoint(
            double[] masterCurCoord, Dictionary<int, double[]> co2Normal,
            uint uQuantityId, int uDof, uint cQuantityId,
            out uint slaveFEId, out double[] slaveL)
        {
            slaveFEId = 0;
            slaveL = null;

            OpenTK.Vector2d masterX = new OpenTK.Vector2d(masterCurCoord[0], masterCurCoord[1]);
            IList<uint> feIds = World.GetContactSlaveLineFEIds(cQuantityId);
            double minGap = double.MaxValue;
            foreach (uint feId in feIds)
            {
                LineFE lineFE = World.GetLineFE(uQuantityId, feId);
                uint elemNodeCnt = lineFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = lineFE.NodeCoordIds[iNode];
                    int nodeId = World.Coord2Node(uQuantityId, coId);
                    nodes[iNode] = nodeId;
                }

                // 現在の頂点の位置
                double[][] curNodeCoords = new double[elemNodeCnt][];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = lineFE.NodeCoordIds[iNode];
                    double[] coord = World.GetCoord(uQuantityId, coId);
                    int iNodeId = nodes[iNode];
                    curNodeCoords[iNode] = new double[uDof];
                    if (iNodeId == -1)
                    {
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            curNodeCoords[iNode][iDof] = coord[iDof];
                        }
                    }
                    else
                    {
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            curNodeCoords[iNode][iDof] = coord[iDof] + U[iNodeId * uDof + iDof];
                        }
                    }
                }
                OpenTK.Vector2d slaveX1 = new OpenTK.Vector2d(
                    curNodeCoords[0][0], curNodeCoords[0][1]);
                OpenTK.Vector2d slaveX2 = new OpenTK.Vector2d(
                    curNodeCoords[1][0], curNodeCoords[1][1]);
                var v1 = slaveX2 - slaveX1;
                var v2 = masterX - slaveX1;
                double[][] nodeNormals = new double[elemNodeCnt][];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = lineFE.NodeCoordIds[iNode];
                    nodeNormals[iNode] = co2Normal[coId];
                }

                double sqNorm = 0;
                double sqInvNorm0 = 0;
                double convRatio = ConvRatioToleranceForNewtonRaphson;
                double tolerance = convRatio;
                const int maxIter = IvyFEM.Linear.Constants.MaxIter;
                int iter = 0;
                double L2 = 0.5;
                OpenTK.Vector2d n = new OpenTK.Vector2d();
                for (iter = 0; iter < maxIter; iter++)
                {
                    double[] L = { 1.0 - L2, L2 };
                    double[] N = lineFE.CalcN(L);
                    double[] normal = new double[uDof];
                    if (L2 >= 0 && L2 <= 1.0)
                    {
                        // 連続な近似法線ベクトルを計算する
                        for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                        {
                            for (int iDof = 0; iDof < uDof; iDof++)
                            {
                                normal[iDof] += nodeNormals[iNode][iDof] * N[iNode];
                            }
                        }
                        // 規格化しない
                        //normal = IvyFEM.Lapack.Utils.NormalizeDoubleVector(normal);
                    }
                    else if (L2 < 0)
                    {
                        // 節点1の法線ベクトルを採用
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            normal[iDof] = nodeNormals[0][iDof];
                        }
                    }
                    else if (L2 > 1)
                    {
                        // 節点2の法線ベクトルを採用
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            normal[iDof] = nodeNormals[1][iDof];
                        }
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    n = new OpenTK.Vector2d(normal[0], normal[1]);
                    n = CadUtils.Normalize(n); // こちらは規格化する

                    var g = v2 - L2 * v1;
                    double[] gap = { g.X, g.Y };
                    // dn/dL2
                    double[] normalL2 = {
                        -nodeNormals[0][0] + nodeNormals[1][0],
                        -nodeNormals[0][1] + nodeNormals[1][1]
                    };
                    // dg/dL2
                    double[] gapL2 = { slaveX1.X - slaveX2.X, slaveX1.Y - slaveX2.Y };
                    double R = normal[0] * gap[1] - normal[1] * gap[0];
                    sqNorm = R * R;
                    if (iter == 0)
                    {
                        if (sqNorm < IvyFEM.Constants.PrecisionLowerLimit)
                        {
                            convRatio = 0;
                            break;
                        }
                        sqInvNorm0 = 1.0 / sqNorm;
                    }
                    else
                    {
                        convRatio = Math.Sqrt(sqNorm * sqInvNorm0);
                        if (sqNorm * sqInvNorm0 < tolerance * tolerance)
                        {
                            break;
                        }
                    }

                    double K = normalL2[0] * gap[1] + normal[0] * gapL2[1]
                        - normalL2[1] * gap[0] - normal[1] * gapL2[0];
                    L2 += -R / K;
                }
                //System.Diagnostics.Debug.WriteLine("GetSlaveLineFEPoint Newton Raphson iter = " + iter + " norm = " + convRatio);
                //System.Diagnostics.Debug.Assert(iter < maxIter);
                if (iter >= maxIter)
                {
                    continue;
                }
                if (L2 >= 0 && L2 <= 1.0)
                {
                    double[] L = { 1.0 - L2, L2 };
                    double[] N = lineFE.CalcN(L);
                    OpenTK.Vector2d slaveX = new OpenTK.Vector2d();
                    for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                    {
                        slaveX.X += N[iNode] * curNodeCoords[iNode][0];
                        slaveX.Y += N[iNode] * curNodeCoords[iNode][1];
                    }
                    double gap = OpenTK.Vector2d.Dot(n, masterX - slaveX);
                    if (gap < minGap)
                    {
                        minGap = gap;
                        slaveFEId = feId;
                        slaveL = L;
                    }
                }
            }
        }

        private void GetMasterLineFEPoint(
            double[] curCoord, double[] normal,
            uint uQuantityId, int uDof, uint cQuantityId,
            out uint masterFEId, out double[] masterL)
        {
            masterFEId = 0;
            masterL = null;
            IList<uint> masterFEIds = World.GetContactMasterLineFEIds(cQuantityId);
            OpenTK.Vector2d slaveX = new OpenTK.Vector2d(curCoord[0], curCoord[1]);
            OpenTK.Vector2d n = new OpenTK.Vector2d(normal[0], normal[1]);
            // t = e3 x n
            OpenTK.Vector2d t = new OpenTK.Vector2d(-normal[1], normal[0]);

            double minGap = double.MaxValue;
            foreach (uint feId in masterFEIds)
            {
                LineFE lineFE = World.GetLineFE(uQuantityId, feId);
                uint elemNodeCnt = lineFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = lineFE.NodeCoordIds[iNode];
                    int nodeId = World.Coord2Node(uQuantityId, coId);
                    nodes[iNode] = nodeId;
                }

                // 現在の頂点の位置
                double[][] masterCurNodeCoords = new double[elemNodeCnt][];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    masterCurNodeCoords[iNode] = new double[uDof];
                    int coId = lineFE.NodeCoordIds[iNode];
                    double[] coord = World.GetCoord(uQuantityId, coId);
                    int iNodeId = nodes[iNode];
                    if (iNodeId == -1)
                    {
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            masterCurNodeCoords[iNode][iDof] = coord[iDof];
                        }
                    }
                    else
                    {
                        for (int iDof = 0; iDof < uDof; iDof++)
                        {
                            masterCurNodeCoords[iNode][iDof] = coord[iDof] + U[iNodeId * uDof + iDof];
                        }
                    }
                }
                OpenTK.Vector2d masterX1 = new OpenTK.Vector2d(
                    masterCurNodeCoords[0][0], masterCurNodeCoords[0][1]);
                OpenTK.Vector2d masterX2 = new OpenTK.Vector2d(
                    masterCurNodeCoords[1][0], masterCurNodeCoords[1][1]);
                var v1 = masterX2 - masterX1;
                var v2 = slaveX - masterX1;
                if (Math.Abs(OpenTK.Vector2d.Dot(v1, t)) < IvyFEM.Constants.PrecisionLowerLimit)
                {
                    continue;
                }
                double L2 = OpenTK.Vector2d.Dot(v2, t) / OpenTK.Vector2d.Dot(v1, t);
                if (L2 >= 0.0 && L2 <= 1.0)
                {
                    // 対象となる要素
                    double[] L = { 1.0 - L2, L2 };
                    OpenTK.Vector2d masterX = new OpenTK.Vector2d(
                        L[0] * masterCurNodeCoords[0][0] + L[1] * masterCurNodeCoords[1][0],
                        L[0] * masterCurNodeCoords[0][1] + L[1] * masterCurNodeCoords[1][1]);
                    double gap = OpenTK.Vector2d.Dot(n, masterX - slaveX);
                    if (gap < minGap) // ギャップの小さい方を採用する
                    {
                        minGap = gap;
                        masterFEId = feId;
                        masterL = L;
                    }
                }
            }
        }
    }
}
