using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class Elastic2DBaseFEM
    {
        protected void SetTwoBodyContactSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            for (uint quantityId = 0; quantityId < World.GetQuantityCount(); quantityId++)
            {
                int slaveCnt = World.GetContactSlaveEIds(quantityId).Count;
                int masterCnt = World.GetContactMasterEIds(quantityId).Count;
                if (slaveCnt > 0 && masterCnt > 0)
                {
                    SetTwoBodyContactMortarSegmentationQuantitySpecialBC(quantityId, A, B);
                }
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
                av = OpenTK.Vector2d.Normalize(av);
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
            System.Diagnostics.Debug.Assert(masterFEIds.Count > 0);

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
                double convRatio = ConvRatioToleranceForNonlinearIter;
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
                    n = OpenTK.Vector2d.Normalize(n); // こちらは規格化する

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
            double[] slaveCurCoord, double[] normal,
            uint uQuantityId, int uDof, uint cQuantityId,
            out uint masterFEId, out double[] masterL)
        {
            masterFEId = 0;
            masterL = null;
            IList<uint> masterFEIds = World.GetContactMasterLineFEIds(cQuantityId);
            OpenTK.Vector2d slaveX = new OpenTK.Vector2d(slaveCurCoord[0], slaveCurCoord[1]);
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
