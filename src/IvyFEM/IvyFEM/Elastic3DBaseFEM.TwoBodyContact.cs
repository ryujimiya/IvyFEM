using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic3DBaseFEM
    {
        protected void SetTwoBodyContactSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            for (uint quantityId = 0; quantityId < World.GetQuantityCount(); quantityId++)
            {
                int slaveCnt = World.GetContactSlaveCadIds(quantityId).Count;
                int masterCnt = World.GetContactMasterCadIds(quantityId).Count;
                if (slaveCnt > 0 && masterCnt > 0)
                {
                    //SetTwoBodyContactNTSSegmentationQuantitySpecialBC(quantityId, A, B);
                    SetTwoBodyContactMortarSegmentationQuantitySpecialBC(quantityId, A, B);
                }
            }
        }

        private Dictionary<int, double[]> GetSlaveTriangleFECo2Normal(uint uQuantityId, int uDof, uint cQuantityId)
        {
            Dictionary<int, double[]> co2Normal = new Dictionary<int, double[]>();
            IList<uint> slaveFEIds = World.GetContactSlaveFEIds(cQuantityId);
            Dictionary<int, IList<double[]>> co2NormalList = new Dictionary<int, IList<double[]>>();
            foreach (uint slaveFEId in slaveFEIds)
            {
                TriangleFE triFE = World.GetTriangleFE(uQuantityId, slaveFEId);
                uint elemNodeCnt = triFE.NodeCount;
                double[] curInsideCo = GetTetrahedronPointNotSharedByTriangle(triFE, uQuantityId, uDof);
                double[] normal = triFE.GetNormal(curInsideCo);

                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = triFE.NodeCoordIds[iNode];
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
                OpenTK.Vector3d av = new OpenTK.Vector3d();
                foreach (double[] normal in normalList)
                {
                    av.X += normal[0];
                    av.Y += normal[1];
                    av.Z += normal[2];
                }
                av = OpenTK.Vector3d.Normalize(av);
                co2Normal[coId] = new double[] { av.X, av.Y, av.Z };
            }
            return co2Normal;
        }

        private double[] GetTetrahedronPointNotSharedByTriangle(TriangleFE triFE, uint uQuantityId, int uDof)
        {
            double[] curInsideCo;
            {
                var triCoIds = triFE.VertexCoordIds;
                int triCoId = triCoIds[0];
                IList<uint> tetFEIds = World.GetTetrahedronFEIdsFromCoord(uQuantityId, triCoId);
                int insideCoId = -1;
                foreach (uint tetFEId in tetFEIds)
                {
                    TetrahedronFE tetFE = World.GetTetrahedronFE(uQuantityId, tetFEId);
                    var tetCoIds = tetFE.VertexCoordIds;
                    IList<int> noSharedCoIds = new List<int>();
                    for (int i = 0; i < tetCoIds.Length; i++)
                    {
                        int tetCoId = tetCoIds[i];
                        if (!triCoIds.Contains(tetCoId))
                        {
                            noSharedCoIds.Add(tetCoId);
                        }
                    }
                    if (noSharedCoIds.Count == 1)
                    {
                        // 三角形2点を共有している四面体
                        insideCoId = noSharedCoIds[0];
                        break;
                    }
                }
                System.Diagnostics.Debug.Assert(insideCoId != -1);
                int insideNodeId = World.Coord2Node(uQuantityId, insideCoId);
                curInsideCo = World.GetCoord(uQuantityId, insideCoId);
                if (insideNodeId != -1)
                {
                    for (int iDof = 0; iDof < uDof; iDof++)
                    {
                        curInsideCo[iDof] += U[insideNodeId * uDof + iDof];
                    }
                }
            }
            return curInsideCo;
        }

        private Dictionary<uint, IList<OpenTK.Vector3d[]>> GetSlaveSegments(
            OpenTK.Vector3d[] slaveTriPts, OpenTK.Vector3d slaveNormal0, 
            uint uQuantityId, int uDof, uint cQuantityId)
        {
            Dictionary<uint, IList<OpenTK.Vector3d[]>> masterFEIdSlaveSegments =
                new Dictionary<uint, IList<OpenTK.Vector3d[]>>();

            OpenTK.Vector3d slaveOrigin = (slaveTriPts[0] + slaveTriPts[1] + slaveTriPts[2]) / 3.0;
            OpenTK.Vector3d slaveNormal = CadUtils3D.TriNormal(slaveTriPts[0], slaveTriPts[1], slaveTriPts[2]);
            OpenTK.Vector3d slaveXDir = CadUtils3D.GetVerticalUnitVector(slaveNormal);

            OpenTK.Vector2d[] slaveTriPt2Ds = new OpenTK.Vector2d[3];
            for (int i = 0; i < 3; i++)
            {
                slaveTriPt2Ds[i] = CadUtils3D.ProjectToPlane(slaveTriPts[i], slaveOrigin, slaveNormal, slaveXDir);
            }

            IList<uint> masterFEIds = World.GetContactMasterFEIds(cQuantityId);
            foreach (uint feId in masterFEIds)
            {
                TriangleFE triFE = World.GetTriangleFE(uQuantityId, feId);
                uint elemNodeCnt = triFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = triFE.NodeCoordIds[iNode];
                    int nodeId = World.Coord2Node(uQuantityId, coId);
                    nodes[iNode] = nodeId;
                }

                double[][] curNodeCoord = new double[elemNodeCnt][];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = triFE.NodeCoordIds[iNode];
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

                OpenTK.Vector3d[] masterTriPts = new OpenTK.Vector3d[] {
                    new OpenTK.Vector3d(curNodeCoord[0][0], curNodeCoord[0][1], curNodeCoord[0][2]),
                    new OpenTK.Vector3d(curNodeCoord[1][0], curNodeCoord[1][1], curNodeCoord[1][2]),
                    new OpenTK.Vector3d(curNodeCoord[2][0], curNodeCoord[2][1], curNodeCoord[2][2]),
                };

                OpenTK.Vector2d[] masterTriPt2Ds = new OpenTK.Vector2d[3];
                for (int i = 0; i < 3; i++)
                {
                    OpenTK.Vector3d lineOrigin = masterTriPts[i];
                    OpenTK.Vector3d lineDir = slaveNormal0;
                    OpenTK.Vector3d planeOrigin = slaveOrigin;
                    OpenTK.Vector3d planeNormal = slaveNormal;
                    OpenTK.Vector3d intersectP;
                    bool isIntersect = CadUtils3D.GetLinePlaneIntersectPoint(
                        lineOrigin, lineDir, planeOrigin, planeNormal, out intersectP);
                    System.Diagnostics.Debug.Assert(isIntersect);
                    masterTriPt2Ds[i] = CadUtils3D.ProjectToPlane(intersectP, slaveOrigin, slaveNormal, slaveXDir);
                }

                OpenTK.Vector2d[] clippedPt2Ds = CadUtils2D.GetTwoPolygonsIntersectPoints(slaveTriPt2Ds, masterTriPt2Ds);
                if (clippedPt2Ds.Length == 0)
                {
                    // 交わり無し
                    continue;
                }
                if (clippedPt2Ds.Length < 3)
                {
                    // ?
                    continue;
                }

                OpenTK.Vector3d[] clippedPts = new OpenTK.Vector3d[clippedPt2Ds.Length];
                for (int i = 0; i < clippedPt2Ds.Length; i++)
                {
                    clippedPts[i] = CadUtils3D.UnProjectFromPlane(clippedPt2Ds[i], slaveOrigin, slaveNormal, slaveXDir);
                }

                IList<OpenTK.Vector3d[]> triSegPtss = GetTriangleSegmentsOfPolygon(clippedPts);

                masterFEIdSlaveSegments.Add(feId, triSegPtss);
            }
            return masterFEIdSlaveSegments;
        }

        private IList<OpenTK.Vector3d[]> GetTriangleSegmentsOfPolygon(OpenTK.Vector3d[] polyPoints)
        {
            IList<OpenTK.Vector3d[]> triSegPtss = new List<OpenTK.Vector3d[]>();
            int ptCnt = polyPoints.Length;
            System.Diagnostics.Debug.Assert(ptCnt >= 3);
            if (ptCnt == 3)
            {
                OpenTK.Vector3d[] triPts = new OpenTK.Vector3d[] {
                    polyPoints[0], polyPoints[1], polyPoints[2]
                };
                triSegPtss.Add(triPts);
                return triSegPtss;
            }
            
            OpenTK.Vector3d mPt = new OpenTK.Vector3d();
            for (int i = 0; i < ptCnt; i++)
            {
                mPt.X += polyPoints[i].X;
                mPt.Y += polyPoints[i].Y;
                mPt.Z += polyPoints[i].Z;
            }
            mPt.X /= ptCnt;
            mPt.Y /= ptCnt;
            mPt.Z /= ptCnt;

            for (int i = 0; i < ptCnt; i++)
            {
                OpenTK.Vector3d[] triPts = new OpenTK.Vector3d[] {
                    polyPoints[i], polyPoints[(i +  1) % ptCnt], mPt
                };
                triSegPtss.Add(triPts);
            }
            return triSegPtss;
        }

        private double[] ModifyL(double[] L)
        {
            System.Diagnostics.Debug.Assert(L.Length == 3);
            double[] modL = new double[3];
            L.CopyTo(modL, 0);

            // 補正
            for (int iL = 0; iL < 3; iL++)
            {
                if (modL[iL] < 0.0)
                {
                    modL[iL] = 0.0;
                }
                if (modL[iL] > 1.0)
                {
                    modL[iL] = 1.0;
                }
            }
            for (int iL = 0; iL < 3; iL++)
            {
                if (modL[iL] + modL[(iL + 1) % 3] > 1.0)
                {
                    modL[(iL + 1) % 3] = 1.0 - modL[iL];
                    modL[(iL + 2) % 3] = 0.0;
                }
            }
            if (modL[0] + modL[1] + modL[2] > 1.0)
            {
                modL[2] = 1.0 - modL[0] - modL[1];
            }
            return modL;
        }

        /*
        private bool GetSlaveTriangleFEPlanePoint(
            uint slaveFEId, Dictionary<int, double[]> co2Normal,
            OpenTK.Vector3d masterPt,
            uint uQuantityId, int uDof, uint cQuantityId,
            out OpenTK.Vector3d slavePt)
        {
            slavePt = new OpenTK.Vector3d();

            TriangleFE triFE = World.GetTriangleFE(uQuantityId, slaveFEId);
            uint elemNodeCnt = triFE.NodeCount;
            int[] nodes = new int[elemNodeCnt];
            for (int iNode = 0; iNode < elemNodeCnt; iNode++)
            {
                int coId = triFE.NodeCoordIds[iNode];
                int nodeId = World.Coord2Node(uQuantityId, coId);
                nodes[iNode] = nodeId;
            }
            // 現在の頂点の位置
            double[][] curNodeCoords = new double[elemNodeCnt][];
            for (int iNode = 0; iNode < elemNodeCnt; iNode++)
            {
                int coId = triFE.NodeCoordIds[iNode];
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
            OpenTK.Vector3d[] triPts = new OpenTK.Vector3d[3] {
                new OpenTK.Vector3d(curNodeCoords[0][0], curNodeCoords[0][1], curNodeCoords[0][2]),
                new OpenTK.Vector3d(curNodeCoords[1][0], curNodeCoords[1][1], curNodeCoords[1][2]),
                new OpenTK.Vector3d(curNodeCoords[2][0], curNodeCoords[2][1], curNodeCoords[2][2])
            };
            double[][] nodeNormals = new double[elemNodeCnt][];
            for (int iNode = 0; iNode < elemNodeCnt; iNode++)
            {
                int coId = triFE.NodeCoordIds[iNode];
                nodeNormals[iNode] = co2Normal[coId];
            }

            double sqNorm = 0;
            double sqInvNorm0 = 0;
            double convRatio = ConvRatioToleranceForNonlinearIter;
            double tolerance = convRatio;
            const int maxIter = IvyFEM.Linear.Constants.MaxIter;
            int iter = 0;
            double[] X = new double[2];
            for (iter = 0; iter < maxIter; iter++)
            {
                double[] L = { X[0], X[1], 1.0 - X[0] - X[1] };
                // 三角形の外の法線ベクトルは一定値とする
                double[] modL = ModifyL(L);

                double[] N = triFE.CalcN(L);
                double[] modN = triFE.CalcN(modL);

                double[] normal = new double[uDof];
                // 連続な近似法線ベクトルを計算する
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    for (int iDof = 0; iDof < uDof; iDof++)
                    {
                        normal[iDof] += nodeNormals[iNode][iDof] * modN[iNode];
                    }
                }
                // 規格化しない
                //normal = IvyFEM.Lapack.Utils.NormalizeDoubleVector(normal);

                OpenTK.Vector3d n = new OpenTK.Vector3d(normal[0], normal[1], normal[2]);
                n = OpenTK.Vector3d.Normalize(n); // こちらは規格化する

                double[,] A0 = new double[3, 2];
                double[] B0 = new double[3];
                OpenTK.Vector3d vec1 = OpenTK.Vector3d.Cross(n, (triPts[0] - triPts[2]));
                OpenTK.Vector3d vec2 = OpenTK.Vector3d.Cross(n, (triPts[1] - triPts[2]));
                OpenTK.Vector3d vec3 = OpenTK.Vector3d.Cross(n, (masterPt - triPts[2]));
                {
                    A0[0, 0] = -vec1.X;
                    A0[0, 1] = -vec2.X;
                    B0[0] = -vec3.X;

                    A0[1, 0] = -vec1.Y;
                    A0[1, 1] = -vec2.Y;
                    B0[1] = -vec3.Y;

                    A0[2, 0] = -vec1.Z;
                    A0[2, 1] = -vec2.Z;
                    B0[2] = -vec3.Z;
                }
                IvyFEM.Linear.DoubleSparseMatrix A = new IvyFEM.Linear.DoubleSparseMatrix(2, 2);
                double[] B = new double[2];
                {
                    int counter = 0;
                    for (int i = 0; i < 3; i++)
                    {
                        bool isZero = true;
                        for (int j = 0; j < 2; j++)
                        {
                            if (Math.Abs(A0[i, j]) >= Constants.PrecisionLowerLimit)
                            {
                                isZero = false;
                                break;
                            }
                        }
                        if (isZero)
                        {
                            continue;
                        }
                        for (int j = 0; j < 2; j++)
                        {
                            A[counter, j] = A0[i, j];
                        }
                        B[counter] = B0[i];
                        counter++;
                        if (counter == 2)
                        {
                            break;
                        }
                    }
                    System.Diagnostics.Debug.Assert(counter == 2);
                }

                double[] AX = A * X;
                double[] R = IvyFEM.Lapack.Functions.daxpy(-1.0, B, AX);
                sqNorm = IvyFEM.Lapack.Functions.ddot(R, R);
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
                var solver = new IvyFEM.Linear.LapackEquationSolver();
                solver.Method = IvyFEM.Linear.LapackEquationSolverMethod.Dense;
                solver.DoubleSolve(out X, A, B);

                L = new double[] { X[0], X[1], 1.0 - X[0] - X[1] };
                N = triFE.CalcN(L);
                slavePt = new OpenTK.Vector3d();
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    slavePt.X += N[iNode] * curNodeCoords[iNode][0];
                    slavePt.Y += N[iNode] * curNodeCoords[iNode][1];
                    slavePt.Z += N[iNode] * curNodeCoords[iNode][2];
                }
            }
            if (iter >= maxIter)
            {
                System.Diagnostics.Debug.Assert(false);
                return false;
            }
            return true;
        }
        */

        private void GetMasterTriangleFEPoint1(
            OpenTK.Vector3d slavePt, OpenTK.Vector3d slaveNormal,
            uint uQuantityId, int uDof, uint cQuantityId,
            out uint masterFEId, out OpenTK.Vector3d masterPt)
        {
            masterFEId = 0;
            masterPt = new OpenTK.Vector3d();
            double minGap = double.MaxValue;
            IList<uint> masterFEIds = World.GetContactMasterFEIds(cQuantityId);
            foreach (uint feId in masterFEIds)
            {
                TriangleFE triFE = World.GetTriangleFE(uQuantityId, feId);
                uint elemNodeCnt = triFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = triFE.NodeCoordIds[iNode];
                    int nodeId = World.Coord2Node(uQuantityId, coId);
                    nodes[iNode] = nodeId;
                }

                double[][] curNodeCoord = new double[elemNodeCnt][];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = triFE.NodeCoordIds[iNode];
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

                OpenTK.Vector3d[] masterTriPts = new OpenTK.Vector3d[] {
                    new OpenTK.Vector3d(curNodeCoord[0][0], curNodeCoord[0][1], curNodeCoord[0][2]),
                    new OpenTK.Vector3d(curNodeCoord[1][0], curNodeCoord[1][1], curNodeCoord[1][2]),
                    new OpenTK.Vector3d(curNodeCoord[2][0], curNodeCoord[2][1], curNodeCoord[2][2]),
                };
                OpenTK.Vector3d intersectPt;
                bool isIntersect = CadUtils3D.GetLineTriangleIntersectPoint(
                    slavePt, slaveNormal, masterTriPts[0], masterTriPts[1], masterTriPts[2], out intersectPt);
                if (isIntersect)
                {
                    double gap = OpenTK.Vector3d.Dot(slaveNormal, masterPt - slavePt);
                    if (gap < minGap) // ギャップの小さい方を採用する
                    {
                        minGap = gap;
                        masterFEId = feId;
                        masterPt = intersectPt;
                    }
                }
            }
        }

        private bool GetMasterTriangleFEPoint2(
            uint masterFEId,
            OpenTK.Vector3d slavePt, OpenTK.Vector3d slaveNormal,
            uint uQuantityId, int uDof, uint cQuantityId,
            out OpenTK.Vector3d masterPt)
        {
            masterPt = new OpenTK.Vector3d();

            TriangleFE triFE = World.GetTriangleFE(uQuantityId, masterFEId);
            uint elemNodeCnt = triFE.NodeCount;
            int[] nodes = new int[elemNodeCnt];
            for (int iNode = 0; iNode < elemNodeCnt; iNode++)
            {
                int coId = triFE.NodeCoordIds[iNode];
                int nodeId = World.Coord2Node(uQuantityId, coId);
                nodes[iNode] = nodeId;
            }

            double[][] curNodeCoord = new double[elemNodeCnt][];
            for (int iNode = 0; iNode < elemNodeCnt; iNode++)
            {
                int coId = triFE.NodeCoordIds[iNode];
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

            OpenTK.Vector3d[] masterTriPts = new OpenTK.Vector3d[] {
                    new OpenTK.Vector3d(curNodeCoord[0][0], curNodeCoord[0][1], curNodeCoord[0][2]),
                    new OpenTK.Vector3d(curNodeCoord[1][0], curNodeCoord[1][1], curNodeCoord[1][2]),
                    new OpenTK.Vector3d(curNodeCoord[2][0], curNodeCoord[2][1], curNodeCoord[2][2]),
                };
            OpenTK.Vector3d intersectPt;
            bool isIntersect = CadUtils3D.GetLineTriangleIntersectPoint(
                slavePt, slaveNormal, masterTriPts[0], masterTriPts[1], masterTriPts[2], out intersectPt);
            if (isIntersect)
            {
                masterPt = intersectPt;
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return isIntersect;
        }
    }
}
