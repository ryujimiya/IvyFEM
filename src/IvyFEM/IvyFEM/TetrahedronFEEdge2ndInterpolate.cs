using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TetrahedronFEEdge2ndInterpolate : IEdgeInterpolate3D
    {
        public TetrahedronFE Owner { get; set; }

        public TetrahedronFEEdge2ndInterpolate()
        {

        }

        public TetrahedronFEEdge2ndInterpolate(TetrahedronFE owner)
        {
            Owner = owner;
        }


        public uint GetNodeCount()
        {
            // 10の頂点
            return 10;
        }

        public uint GetEdgeCount()
        {
            // 20の辺
            return 20;
        }

        public double[] GetNodeL(int nodeId)
        {
            throw new NotImplementedException();
        }

        public double[] CalcN(double[] L)
        {
            throw new NotImplementedException();
        }

        public double[][] CalcNu(double[] L)
        {
            throw new NotImplementedException();
        }

        public double[,][] CalcNuv(double[] L)
        {
            throw new NotImplementedException();
        }

        public double[] CalcSN()
        {
            throw new NotImplementedException();
        }

        public double[,] CalcSNN()
        {
            throw new NotImplementedException();
        }

        public double[,][,] CalcSNuNv()
        {
            throw new NotImplementedException();
        }

        public double[][,] CalcSNuN()
        {
            throw new NotImplementedException();
        }

        public int[][] GetEdgePointIdss()
        {
            const int faceCnt = 4;
            int[][] facePtNo = new int[faceCnt][]
            {
                new int[] { 1, 2, 3 },
                new int[] { 0, 2, 3 },
                new int[] { 0, 1, 3 },
                new int[] { 0, 1, 2 }
            };

            IList<int>[] vNoss = new List<int>[faceCnt];
            for (int iface = 0; iface < faceCnt; iface++)
            {
                var list = new List<KeyValuePair<int, int>>();
                for (int i = 0; i < 3; i++)
                {
                    int ptNo = facePtNo[iface][i];
                    list.Add(new KeyValuePair<int, int>(Owner.VertexCoordIds[ptNo], ptNo));
                }
                list.Sort((a, b) =>
                {
                    double diff = a.Key - b.Key;
                    if (diff > 0)
                    {
                        return 1;
                    }
                    else if (diff < 0)
                    {
                        return -1;
                    }
                    return 0;
                });

                IList<int> vNos = new List<int>();
                vNoss[iface] = vNos;
                foreach (var pair in list)
                {
                    int coId = pair.Key;
                    int vNo = pair.Value;
                    vNos.Add(vNo);
                }
            }

            int[][] edgePointId = new int[20][]
            {
                new int[]{ 0, 4 },
                new int[]{ 1, 4 },
                new int[]{ 0, 6 },
                new int[]{ 2, 6 },
                new int[]{ 0, 7 },
                new int[]{ 3, 7 },
                new int[]{ 1, 5 },
                new int[]{ 2, 5 },
                new int[]{ 3, 8 },
                new int[]{ 1, 8 },
                new int[]{ 2, 9 },
                new int[]{ 3, 9 },
                null,
                null,
                null,
                null,
                null,
                null,
                null,
                null
            };

            {
                int iface = 0;
                var vNos = vNoss[iface];
                for (int eIndex = 12; eIndex < 14; eIndex++)
                {
                    int vNo = vNos[eIndex - 12];
                    int[] pointIds = null;
                    if (vNo == 1)
                    {
                        pointIds = new int[] { 9, 1 };
                    }
                    else if (vNo == 2)
                    {
                        pointIds = new int[] { 8, 2 };
                    }
                    else if (vNo == 3)
                    {
                        pointIds = new int[] { 5, 3 };
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    edgePointId[eIndex] = pointIds;
                }
            }
            {
                int iface = 1;
                var vNos = vNoss[iface];
                for (int eIndex = 14; eIndex < 16; eIndex++)
                {
                    int vNo = vNos[eIndex - 14];
                    int[] pointIds = null;
                    if (vNo == 0)
                    {
                        pointIds = new int[] { 9, 0 };
                    }
                    else if (vNo == 2)
                    {
                        pointIds = new int[] { 7, 2 };
                    }
                    else if (vNo == 3)
                    {
                        pointIds = new int[] { 6, 3 };
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    edgePointId[eIndex] = pointIds;
                }
            }
            {
                int iface = 2;
                var vNos = vNoss[iface];
                for (int eIndex = 16; eIndex < 18; eIndex++)
                {
                    int vNo = vNos[eIndex - 16];
                    int[] pointIds = null;
                    if (vNo == 0)
                    {
                        pointIds = new int[] { 8, 0 };
                    }
                    else if (vNo == 1)
                    {
                        pointIds = new int[] { 7, 1 };
                    }
                    else if (vNo == 3)
                    {
                        pointIds = new int[] { 4, 3 };
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    edgePointId[eIndex] = pointIds;
                }
            }
            {
                int iface = 3;
                var vNos = vNoss[iface];
                for (int eIndex = 18; eIndex < 20; eIndex++)
                {
                    int vNo = vNos[eIndex - 18];
                    int[] pointIds = null;
                    if (vNo == 0)
                    {
                        pointIds = new int[] { 5, 0 };
                    }
                    else if (vNo == 1)
                    {
                        pointIds = new int[] { 6, 1 };
                    }
                    else if (vNo == 2)
                    {
                        pointIds = new int[] { 4, 2 };
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    edgePointId[eIndex] = pointIds;
                }
            }

            return edgePointId;
        }

        private int[][] GetLongEdgePointIdss()
        {
            int[][] edgePointId = new int[6][]
            {
                new int[]{ 0, 1 },
                new int[]{ 0, 2 },
                new int[]{ 0, 3 },
                new int[]{ 1, 2 },
                new int[]{ 3, 1 },
                new int[]{ 2, 3 }
            };
            return edgePointId;
        }

        private void CalcLens(out double[] lens)
        {
            int[][] longEdgePointId = GetLongEdgePointIdss();

            lens = new double[6];
            for (int eIndex = 0; eIndex < 6; eIndex++)
            {
                int[] edgePointId1 = longEdgePointId[eIndex];
                int iPtId = edgePointId1[0];
                int jPtId = edgePointId1[1];
                int coId1 = Owner.VertexCoordIds[iPtId];
                int coId2 = Owner.VertexCoordIds[jPtId];
                double[] co1 = Owner.World.GetCoord((uint)Owner.QuantityId, coId1);
                double[] co2 = Owner.World.GetCoord((uint)Owner.QuantityId, coId2);
                lens[eIndex] = Math.Sqrt(
                    (co2[0] - co1[0]) * (co2[0] - co1[0]) +
                    (co2[1] - co1[1]) * (co2[1] - co1[1]) +
                    (co2[2] - co1[2]) * (co2[2] - co1[2]));
            }
        }

        private void CalcGradLs(
            double[] a, double[] b, double[]c, double[] d,
            out double[][] gradLs)
        {
            gradLs = new double[4][];
            for (int i = 0; i < 3; i++)
            {
                double[] gradL = new double[3];
                gradLs[i] = gradL;

                gradL[0] = b[i];
                gradL[1] = c[i];
                gradL[2] = d[i];
            }
            {
                int i = 3;
                double[] gradL = new double[3];
                gradLs[i] = gradL;

                for (int idim = 0; idim < 3; idim++)
                {
                    double[] gradL0 = gradLs[0];
                    double[] gradL1 = gradLs[1];
                    double[] gradL2 = gradLs[2];
                    gradL[idim] = -gradL0[idim] - gradL1[idim] - gradL2[idim];
                }
            }
        }

        // gradLを使う
        // FIXME: 三角形要素と整合した規格化の実施
        private void CalcNVecss(
            double[][] gradLs,
            out OpenTK.Vector3d[][] nVecss)
        {
            var gradL1Vec = new OpenTK.Vector3d(gradLs[0][0], gradLs[0][1], gradLs[0][2]);
            var gradL2Vec = new OpenTK.Vector3d(gradLs[1][0], gradLs[1][1], gradLs[1][2]);
            var gradL3Vec = new OpenTK.Vector3d(gradLs[2][0], gradLs[2][1], gradLs[2][2]);
            var gradL4Vec = new OpenTK.Vector3d(gradLs[3][0], gradLs[3][1], gradLs[3][2]);

            nVecss = new OpenTK.Vector3d[4][];
            {
                int iface = 0;
                var n2 = new OpenTK.Vector3d(gradL2Vec);
                var n3 = new OpenTK.Vector3d(gradL3Vec);
                var n4 = new OpenTK.Vector3d(gradL4Vec);
                nVecss[iface] = new OpenTK.Vector3d[] { n2, n3, n4 };
            }
            {
                int iface = 1;
                var n1 = new OpenTK.Vector3d(gradL1Vec);
                var n3 = new OpenTK.Vector3d(gradL3Vec);
                var n4 = new OpenTK.Vector3d(gradL4Vec);
                nVecss[iface] = new OpenTK.Vector3d[] { n1, n3, n4 };
            }
            {
                int iface = 2;
                var n1 = new OpenTK.Vector3d(gradL1Vec);
                var n2 = new OpenTK.Vector3d(gradL2Vec);
                var n4 = new OpenTK.Vector3d(gradL4Vec);
                nVecss[iface] = new OpenTK.Vector3d[] { n1, n2, n4 };
            }
            {
                int iface = 3;
                var n1 = new OpenTK.Vector3d(gradL1Vec);
                var n2 = new OpenTK.Vector3d(gradL2Vec);
                var n3 = new OpenTK.Vector3d(gradL3Vec);
                nVecss[iface] = new OpenTK.Vector3d[] { n1, n2, n3 };
            }
        }

        public double[][] CalcEdgeN(double[] L)
        {
            double[][] edgeNs = new double[20][];
            double[] a;
            double[] b;
            double[] c;
            double[] d;
            Owner.CalcTransMatrix(out a, out b, out c, out d);

            int[][] edgePointId = GetEdgePointIdss();

            double[] lens;
            CalcLens(out lens);

            double[][] gradLs;
            CalcGradLs(a, b, c, d, out gradLs);

            OpenTK.Vector3d[][] nVecss;
            CalcNVecss(gradLs, out nVecss);

            {
                int eIndex = 0;
                double[] edgeN = new double[3];
                edgeNs[eIndex] = edgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    edgeN[idim] = lens[0] * L[0] * gradLs[1][idim];
                }
            }
            {
                int eIndex = 1;
                double[] edgeN = new double[3];
                edgeNs[eIndex] = edgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    edgeN[idim] = lens[0] * L[1] * gradLs[0][idim];
                }
            }
            {
                int eIndex = 2;
                double[] edgeN = new double[3];
                edgeNs[eIndex] = edgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    edgeN[idim] = lens[1] * L[0] * gradLs[2][idim];
                }
            }
            {
                int eIndex = 3;
                double[] edgeN = new double[3];
                edgeNs[eIndex] = edgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    edgeN[idim] = lens[1] * L[2] * gradLs[0][idim];
                }
            }
            {
                int eIndex = 4;
                double[] edgeN = new double[3];
                edgeNs[eIndex] = edgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    edgeN[idim] = lens[2] * L[0] * gradLs[3][idim];
                }
            }
            {
                int eIndex = 5;
                double[] edgeN = new double[3];
                edgeNs[eIndex] = edgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    edgeN[idim] = lens[2] * L[3] * gradLs[0][idim];
                }
            }
            {
                int eIndex = 6;
                double[] edgeN = new double[3];
                edgeNs[eIndex] = edgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    edgeN[idim] = lens[3] * L[1] * gradLs[2][idim];
                }
            }
            {
                int eIndex = 7;
                double[] edgeN = new double[3];
                edgeNs[eIndex] = edgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    edgeN[idim] = lens[3] * L[2] * gradLs[1][idim];
                }
            }
            {
                int eIndex = 8;
                double[] edgeN = new double[3];
                edgeNs[eIndex] = edgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    edgeN[idim] = lens[4] * L[3] * gradLs[1][idim];
                }
            }
            {
                int eIndex = 9;
                double[] edgeN = new double[3];
                edgeNs[eIndex] = edgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    edgeN[idim] = lens[4] * L[1] * gradLs[3][idim];
                }
            }
            {
                int eIndex = 10;
                double[] edgeN = new double[3];
                edgeNs[eIndex] = edgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    edgeN[idim] = lens[5] * L[2] * gradLs[3][idim];
                }
            }
            {
                int eIndex = 11;
                double[] edgeN = new double[3];
                edgeNs[eIndex] = edgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    edgeN[idim] = lens[5] * L[3] * gradLs[2][idim];
                }
            }

            {
                for (int iVNo = 0; iVNo < 2; iVNo++)
                {
                    int eIndex = 12 + iVNo;
                    int vIndex = edgePointId[eIndex][1];
                    System.Diagnostics.Debug.Assert(vIndex >= 0 && vIndex < 4);
                    double[] edgeN = new double[3];
                    edgeNs[eIndex] = edgeN;

                    int iface = 0;
                    double[] n2 = new double[] { nVecss[iface][0].X, nVecss[iface][0].Y, nVecss[iface][0].Z };
                    double[] n3 = new double[] { nVecss[iface][1].X, nVecss[iface][1].Y, nVecss[iface][1].Z };
                    double[] n4 = new double[] { nVecss[iface][2].X, nVecss[iface][2].Y, nVecss[iface][2].Z };

                    for (int idim = 0; idim < 3; idim++)
                    {
                        if (vIndex == 1)
                        {
                            edgeN[idim] = 4.0 * L[2] * L[3] * n2[idim];
                        }
                        else if (vIndex == 2)
                        {
                            edgeN[idim] = 4.0 * L[1] * L[3] * n3[idim];
                        }
                        else if (vIndex == 3)
                        {
                            edgeN[idim] = 4.0 * L[1] * L[2] * n4[idim];
                        }
                        else
                        {
                            System.Diagnostics.Debug.Assert(false);
                        }
                    }
                }
            }
            {
                for (int iVNo = 0; iVNo < 2; iVNo++)
                {
                    int eIndex = 14 + iVNo;
                    int vIndex = edgePointId[eIndex][1];
                    System.Diagnostics.Debug.Assert(vIndex >= 0 && vIndex < 4);
                    double[] edgeN = new double[3];
                    edgeNs[eIndex] = edgeN;

                    int iface = 1;
                    double[] n1 = new double[] { nVecss[iface][0].X, nVecss[iface][0].Y, nVecss[iface][0].Z };
                    double[] n3 = new double[] { nVecss[iface][1].X, nVecss[iface][1].Y, nVecss[iface][1].Z };
                    double[] n4 = new double[] { nVecss[iface][2].X, nVecss[iface][2].Y, nVecss[iface][2].Z };

                    for (int idim = 0; idim < 3; idim++)
                    {
                        if (vIndex == 0)
                        {
                            edgeN[idim] = 4.0 * L[2] * L[3] * n1[idim];
                        }
                        else if (vIndex == 2)
                        {
                            edgeN[idim] = 4.0 * L[0] * L[3] * n3[idim];
                        }
                        else if (vIndex == 3)
                        {
                            edgeN[idim] = 4.0 * L[0] * L[2] * n4[idim];
                        }
                        else
                        {
                            System.Diagnostics.Debug.Assert(false);
                        }
                    }
                }
            }
            {
                for (int iVNo = 0; iVNo < 2; iVNo++)
                {
                    int eIndex = 16 + iVNo;
                    int vIndex = edgePointId[eIndex][1];
                    System.Diagnostics.Debug.Assert(vIndex >= 0 && vIndex < 4);
                    double[] edgeN = new double[3];
                    edgeNs[eIndex] = edgeN;

                    int iface = 2;
                    double[] n1 = new double[] { nVecss[iface][0].X, nVecss[iface][0].Y, nVecss[iface][0].Z };
                    double[] n2 = new double[] { nVecss[iface][1].X, nVecss[iface][1].Y, nVecss[iface][1].Z };
                    double[] n4 = new double[] { nVecss[iface][2].X, nVecss[iface][2].Y, nVecss[iface][2].Z };

                    for (int idim = 0; idim < 3; idim++)
                    {
                        if (vIndex == 0)
                        {
                            edgeN[idim] = 4.0 * L[1] * L[3] * n1[idim];
                        }
                        else if (vIndex == 1)
                        {
                            edgeN[idim] = 4.0 * L[0] * L[3] * n2[idim];
                        }
                        else if (vIndex == 3)
                        {
                            edgeN[idim] = 4.0 * L[0] * L[1] * n4[idim];
                        }
                        else
                        {
                            System.Diagnostics.Debug.Assert(false);
                        }
                    }
                }
            }
            {
                for (int iVNo = 0; iVNo < 2; iVNo++)
                {
                    int eIndex = 18 + iVNo;
                    int vIndex = edgePointId[eIndex][1];
                    System.Diagnostics.Debug.Assert(vIndex >= 0 && vIndex < 4);
                    double[] edgeN = new double[3];
                    edgeNs[eIndex] = edgeN;

                    int iface = 3;
                    double[] n1 = new double[] { nVecss[iface][0].X, nVecss[iface][0].Y, nVecss[iface][0].Z };
                    double[] n2 = new double[] { nVecss[iface][1].X, nVecss[iface][1].Y, nVecss[iface][1].Z };
                    double[] n3 = new double[] { nVecss[iface][2].X, nVecss[iface][2].Y, nVecss[iface][2].Z };

                    for (int idim = 0; idim < 3; idim++)
                    {
                        if (vIndex == 0)
                        {
                            edgeN[idim] = 4.0 * L[1] * L[2] * n1[idim];
                        }
                        else if (vIndex == 1)
                        {
                            edgeN[idim] = 4.0 * L[0] * L[2] * n2[idim];
                        }
                        else if (vIndex == 2)
                        {
                            edgeN[idim] = 4.0 * L[0] * L[1] * n3[idim];
                        }
                        else
                        {
                            System.Diagnostics.Debug.Assert(false);
                        }
                    }
                }
            }

            return edgeNs;
        }

        public double[][] CalcRotEdgeN(double[] L)
        {
            double[][] rotEdgeNs = new double[20][];
            double[] a;
            double[] b;
            double[] c;
            double[] d;
            Owner.CalcTransMatrix(out a, out b, out c, out d);

            int[][] edgePointId = GetEdgePointIdss();

            double[] lens;
            CalcLens(out lens);

            double[][] gradLs;
            CalcGradLs(a, b, c, d, out gradLs);

            OpenTK.Vector3d[][] nVecss;
            CalcNVecss(gradLs, out nVecss);

            {
                int eIndex = 0;
                int iPt1 = 0;
                int iPt2 = 1;
                OpenTK.Vector3d gradLVec1 = new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                OpenTK.Vector3d gradLVec2 = new OpenTK.Vector3d(gradLs[iPt2][0], gradLs[iPt2][1], gradLs[iPt2][2]);
                var gradL1xgrapdL2Vec = OpenTK.Vector3d.Cross(gradLVec1, gradLVec2);
                double[] gradL1xgradL2 = { gradL1xgrapdL2Vec.X, gradL1xgrapdL2Vec.Y, gradL1xgrapdL2Vec.Z };

                double[] rotEdgeN = new double[3];
                rotEdgeNs[eIndex] = rotEdgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    rotEdgeN[idim] = lens[0] * gradL1xgradL2[idim];
                }
            }
            {
                int eIndex = 1;
                int iPt1 = 1;
                int iPt2 = 0;
                OpenTK.Vector3d gradLVec1 = new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                OpenTK.Vector3d gradLVec2 = new OpenTK.Vector3d(gradLs[iPt2][0], gradLs[iPt2][1], gradLs[iPt2][2]);
                var gradL1xgrapdL2Vec = OpenTK.Vector3d.Cross(gradLVec1, gradLVec2);
                double[] gradL1xgradL2 = { gradL1xgrapdL2Vec.X, gradL1xgrapdL2Vec.Y, gradL1xgrapdL2Vec.Z };

                double[] rotEdgeN = new double[3];
                rotEdgeNs[eIndex] = rotEdgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    rotEdgeN[idim] = lens[0] * gradL1xgradL2[idim];
                }
            }
            {
                int eIndex = 2;
                int iPt1 = 0;
                int iPt2 = 2;
                OpenTK.Vector3d gradLVec1 = new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                OpenTK.Vector3d gradLVec2 = new OpenTK.Vector3d(gradLs[iPt2][0], gradLs[iPt2][1], gradLs[iPt2][2]);
                var gradL1xgrapdL2Vec = OpenTK.Vector3d.Cross(gradLVec1, gradLVec2);
                double[] gradL1xgradL2 = { gradL1xgrapdL2Vec.X, gradL1xgrapdL2Vec.Y, gradL1xgrapdL2Vec.Z };

                double[] rotEdgeN = new double[3];
                rotEdgeNs[eIndex] = rotEdgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    rotEdgeN[idim] = lens[1] * gradL1xgradL2[idim];
                }
            }
            {
                int eIndex = 3;
                int iPt1 = 2;
                int iPt2 = 0;
                OpenTK.Vector3d gradLVec1 = new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                OpenTK.Vector3d gradLVec2 = new OpenTK.Vector3d(gradLs[iPt2][0], gradLs[iPt2][1], gradLs[iPt2][2]);
                var gradL1xgrapdL2Vec = OpenTK.Vector3d.Cross(gradLVec1, gradLVec2);
                double[] gradL1xgradL2 = { gradL1xgrapdL2Vec.X, gradL1xgrapdL2Vec.Y, gradL1xgrapdL2Vec.Z };

                double[] rotEdgeN = new double[3];
                rotEdgeNs[eIndex] = rotEdgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    rotEdgeN[idim] = lens[1] * gradL1xgradL2[idim];
                }
            }
            {
                int eIndex = 4;
                int iPt1 = 0;
                int iPt2 = 3;
                OpenTK.Vector3d gradLVec1 = new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                OpenTK.Vector3d gradLVec2 = new OpenTK.Vector3d(gradLs[iPt2][0], gradLs[iPt2][1], gradLs[iPt2][2]);
                var gradL1xgrapdL2Vec = OpenTK.Vector3d.Cross(gradLVec1, gradLVec2);
                double[] gradL1xgradL2 = { gradL1xgrapdL2Vec.X, gradL1xgrapdL2Vec.Y, gradL1xgrapdL2Vec.Z };

                double[] rotEdgeN = new double[3];
                rotEdgeNs[eIndex] = rotEdgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    rotEdgeN[idim] = lens[2] * gradL1xgradL2[idim];
                }
            }
            {
                int eIndex = 5;
                int iPt1 = 3;
                int iPt2 = 0;
                OpenTK.Vector3d gradLVec1 = new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                OpenTK.Vector3d gradLVec2 = new OpenTK.Vector3d(gradLs[iPt2][0], gradLs[iPt2][1], gradLs[iPt2][2]);
                var gradL1xgrapdL2Vec = OpenTK.Vector3d.Cross(gradLVec1, gradLVec2);
                double[] gradL1xgradL2 = { gradL1xgrapdL2Vec.X, gradL1xgrapdL2Vec.Y, gradL1xgrapdL2Vec.Z };

                double[] rotEdgeN = new double[3];
                rotEdgeNs[eIndex] = rotEdgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    rotEdgeN[idim] = lens[2] * gradL1xgradL2[idim];
                }
            }
            {
                int eIndex = 6;
                int iPt1 = 1;
                int iPt2 = 2;
                OpenTK.Vector3d gradLVec1 = new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                OpenTK.Vector3d gradLVec2 = new OpenTK.Vector3d(gradLs[iPt2][0], gradLs[iPt2][1], gradLs[iPt2][2]);
                var gradL1xgrapdL2Vec = OpenTK.Vector3d.Cross(gradLVec1, gradLVec2);
                double[] gradL1xgradL2 = { gradL1xgrapdL2Vec.X, gradL1xgrapdL2Vec.Y, gradL1xgrapdL2Vec.Z };

                double[] rotEdgeN = new double[3];
                rotEdgeNs[eIndex] = rotEdgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    rotEdgeN[idim] = lens[3] * gradL1xgradL2[idim];
                }
            }
            {
                int eIndex = 7;
                int iPt1 = 2;
                int iPt2 = 1;
                OpenTK.Vector3d gradLVec1 = new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                OpenTK.Vector3d gradLVec2 = new OpenTK.Vector3d(gradLs[iPt2][0], gradLs[iPt2][1], gradLs[iPt2][2]);
                var gradL1xgrapdL2Vec = OpenTK.Vector3d.Cross(gradLVec1, gradLVec2);
                double[] gradL1xgradL2 = { gradL1xgrapdL2Vec.X, gradL1xgrapdL2Vec.Y, gradL1xgrapdL2Vec.Z };

                double[] rotEdgeN = new double[3];
                rotEdgeNs[eIndex] = rotEdgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    rotEdgeN[idim] = lens[3] * gradL1xgradL2[idim];
                }
            }
            {
                int eIndex = 8;
                int iPt1 = 3;
                int iPt2 = 1;
                OpenTK.Vector3d gradLVec1 = new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                OpenTK.Vector3d gradLVec2 = new OpenTK.Vector3d(gradLs[iPt2][0], gradLs[iPt2][1], gradLs[iPt2][2]);
                var gradL1xgrapdL2Vec = OpenTK.Vector3d.Cross(gradLVec1, gradLVec2);
                double[] gradL1xgradL2 = { gradL1xgrapdL2Vec.X, gradL1xgrapdL2Vec.Y, gradL1xgrapdL2Vec.Z };

                double[] rotEdgeN = new double[3];
                rotEdgeNs[eIndex] = rotEdgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    rotEdgeN[idim] = lens[4] * gradL1xgradL2[idim];
                }
            }
            {
                int eIndex = 9;
                int iPt1 = 1;
                int iPt2 = 3;
                OpenTK.Vector3d gradLVec1 = new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                OpenTK.Vector3d gradLVec2 = new OpenTK.Vector3d(gradLs[iPt2][0], gradLs[iPt2][1], gradLs[iPt2][2]);
                var gradL1xgrapdL2Vec = OpenTK.Vector3d.Cross(gradLVec1, gradLVec2);
                double[] gradL1xgradL2 = { gradL1xgrapdL2Vec.X, gradL1xgrapdL2Vec.Y, gradL1xgrapdL2Vec.Z };

                double[] rotEdgeN = new double[3];
                rotEdgeNs[eIndex] = rotEdgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    rotEdgeN[idim] = lens[4] * gradL1xgradL2[idim];
                }
            }
            {
                int eIndex = 10;
                int iPt1 = 2;
                int iPt2 = 3;
                OpenTK.Vector3d gradLVec1 = new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                OpenTK.Vector3d gradLVec2 = new OpenTK.Vector3d(gradLs[iPt2][0], gradLs[iPt2][1], gradLs[iPt2][2]);
                var gradL1xgrapdL2Vec = OpenTK.Vector3d.Cross(gradLVec1, gradLVec2);
                double[] gradL1xgradL2 = { gradL1xgrapdL2Vec.X, gradL1xgrapdL2Vec.Y, gradL1xgrapdL2Vec.Z };

                double[] rotEdgeN = new double[3];
                rotEdgeNs[eIndex] = rotEdgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    rotEdgeN[idim] = lens[5] * gradL1xgradL2[idim];
                }
            }
            {
                int eIndex = 11;
                int iPt1 = 3;
                int iPt2 = 2;
                OpenTK.Vector3d gradLVec1 = new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                OpenTK.Vector3d gradLVec2 = new OpenTK.Vector3d(gradLs[iPt2][0], gradLs[iPt2][1], gradLs[iPt2][2]);
                var gradL1xgrapdL2Vec = OpenTK.Vector3d.Cross(gradLVec1, gradLVec2);
                double[] gradL1xgradL2 = { gradL1xgrapdL2Vec.X, gradL1xgrapdL2Vec.Y, gradL1xgrapdL2Vec.Z };

                double[] rotEdgeN = new double[3];
                rotEdgeNs[eIndex] = rotEdgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    rotEdgeN[idim] = lens[5] * gradL1xgradL2[idim];
                }
            }

            {
                for (int iVNo = 0; iVNo < 2; iVNo++)
                {
                    int eIndex = 12 + iVNo;
                    int vIndex = edgePointId[eIndex][1];
                    System.Diagnostics.Debug.Assert(vIndex >= 0 && vIndex < 4);
                    double[] rotEdgeN = new double[3];
                    rotEdgeNs[eIndex] = rotEdgeN;
                    double L1 = 0.0;
                    double L2 = 0.0;
                    double[] gradL3xn4 = null;
                    double[] gradL5xn6 = null;

                    int iface = 0;
                    var n2 = nVecss[iface][0];
                    var n3 = nVecss[iface][1];
                    var n4 = nVecss[iface][2];

                    if (vIndex == 1)
                    {
                        L1 = L[2];
                        L2 = L[3];
                        {
                            int iPt1 = 3;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n2);
                            gradL3xn4 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                        {
                            int iPt1 = 2;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n2);
                            gradL5xn6 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                    }
                    else if (vIndex == 2)
                    {
                        L1 = L[1];
                        L2 = L[3];
                        {
                            int iPt1 = 3;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n3);
                            gradL3xn4 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                        {
                            int iPt1 = 1;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n3);
                            gradL5xn6 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                    }
                    else if (vIndex == 3)
                    {
                        L1 = L[1];
                        L2 = L[2];
                        {
                            int iPt1 = 2;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n4);
                            gradL3xn4 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                        {
                            int iPt1 = 1;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n4);
                            gradL5xn6 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }

                    for (int idim = 0; idim < 3; idim++)
                    {
                        rotEdgeN[idim] = 4.0 * L1 * gradL3xn4[idim] + 4.0 * L2 * gradL5xn6[idim];
                    }
                }
            }
            {
                for (int iVNo = 0; iVNo < 2; iVNo++)
                {
                    int eIndex = 14 + iVNo;
                    int vIndex = edgePointId[eIndex][1];
                    System.Diagnostics.Debug.Assert(vIndex >= 0 && vIndex < 4);
                    double[] rotEdgeN = new double[3];
                    rotEdgeNs[eIndex] = rotEdgeN;
                    double L1 = 0.0;
                    double L2 = 0.0;
                    double[] gradL3xn4 = null;
                    double[] gradL5xn6 = null;

                    int iface = 1;
                    var n1 = nVecss[iface][0];
                    var n3 = nVecss[iface][1];
                    var n4 = nVecss[iface][2];

                    if (vIndex == 0)
                    {
                        L1 = L[2];
                        L2 = L[3];
                        {
                            int iPt1 = 3;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n1);
                            gradL3xn4 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                        {
                            int iPt1 = 2;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n1);
                            gradL5xn6 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                    }
                    else if (vIndex == 2)
                    {
                        L1 = L[0];
                        L2 = L[3];
                        {
                            int iPt1 = 3;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n3);
                            gradL3xn4 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                        {
                            int iPt1 = 0;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n3);
                            gradL5xn6 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                    }
                    else if (vIndex == 3)
                    {
                        L1 = L[0];
                        L2 = L[2];
                        {
                            int iPt1 = 2;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n4);
                            gradL3xn4 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                        {
                            int iPt1 = 0;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n4);
                            gradL5xn6 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }

                    for (int idim = 0; idim < 3; idim++)
                    {
                        rotEdgeN[idim] = 4.0 * L1 * gradL3xn4[idim] + 4.0 * L2 * gradL5xn6[idim];
                    }
                }
            }
            {
                for (int iVNo = 0; iVNo < 2; iVNo++)
                {
                    int eIndex = 16 + iVNo;
                    int vIndex = edgePointId[eIndex][1];
                    System.Diagnostics.Debug.Assert(vIndex >= 0 && vIndex < 4);
                    double[] rotEdgeN = new double[3];
                    rotEdgeNs[eIndex] = rotEdgeN;
                    double L1 = 0.0;
                    double L2 = 0.0;
                    double[] gradL3xn4 = null;
                    double[] gradL5xn6 = null;

                    int iface = 2;
                    var n1 = nVecss[iface][0];
                    var n2 = nVecss[iface][1];
                    var n4 = nVecss[iface][2];

                    if (vIndex == 0)
                    {
                        L1 = L[1];
                        L2 = L[3];
                        {
                            int iPt1 = 3;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n1);
                            gradL3xn4 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                        {
                            int iPt1 = 1;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n1);
                            gradL5xn6 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                    }
                    else if (vIndex == 1)
                    {
                        L1 = L[0];
                        L2 = L[3];
                        {
                            int iPt1 = 3;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n2);
                            gradL3xn4 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                        {
                            int iPt1 = 0;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n2);
                            gradL5xn6 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                    }
                    else if (vIndex == 3)
                    {
                        L1 = L[0];
                        L2 = L[1];
                        {
                            int iPt1 = 1;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n4);
                            gradL3xn4 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                        {
                            int iPt1 = 0;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n4);
                            gradL5xn6 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }

                    for (int idim = 0; idim < 3; idim++)
                    {
                        rotEdgeN[idim] = 4.0 * L1 * gradL3xn4[idim] + 4.0 * L2 * gradL5xn6[idim];
                    }
                }
            }
            {
                for (int iVNo = 0; iVNo < 2; iVNo++)
                {
                    int eIndex = 18 + iVNo;
                    int vIndex = edgePointId[eIndex][1];
                    System.Diagnostics.Debug.Assert(vIndex >= 0 && vIndex < 4);
                    double[] rotEdgeN = new double[3];
                    rotEdgeNs[eIndex] = rotEdgeN;
                    double L1 = 0.0;
                    double L2 = 0.0;
                    double[] gradL3xn4 = null;
                    double[] gradL5xn6 = null;

                    int iface = 3;
                    var n1 = nVecss[iface][0];
                    var n2 = nVecss[iface][1];
                    var n3 = nVecss[iface][2];

                    if (vIndex == 0)
                    {
                        L1 = L[1];
                        L2 = L[2];
                        {
                            int iPt1 = 2;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n1);
                            gradL3xn4 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                        {
                            int iPt1 = 1;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n1);
                            gradL5xn6 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                    }
                    else if (vIndex == 1)
                    {
                        L1 = L[0];
                        L2 = L[2];
                        {
                            int iPt1 = 2;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n2);
                            gradL3xn4 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                        {
                            int iPt1 = 0;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n2);
                            gradL5xn6 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                    }
                    else if (vIndex == 2)
                    {
                        L1 = L[0];
                        L2 = L[1];
                        {
                            int iPt1 = 1;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n3);
                            gradL3xn4 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                        {
                            int iPt1 = 0;
                            OpenTK.Vector3d gradLVec1 =
                                new OpenTK.Vector3d(gradLs[iPt1][0], gradLs[iPt1][1], gradLs[iPt1][2]);
                            var gradL1xn2Vec = OpenTK.Vector3d.Cross(gradLVec1, n3);
                            gradL5xn6 = new double[] { gradL1xn2Vec.X, gradL1xn2Vec.Y, gradL1xn2Vec.Z };
                        }
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }

                    for (int idim = 0; idim < 3; idim++)
                    {
                        rotEdgeN[idim] = 4.0 * L1 * gradL3xn4[idim] + 4.0 * L2 * gradL5xn6[idim];
                    }
                }
            }

            return rotEdgeNs;
        }
    }
}
