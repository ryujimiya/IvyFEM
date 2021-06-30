using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TriangleFEEdge2ndInterpolate : IEdgeInterpolate
    {
        public TriangleFE Owner { get; set; }

        public TriangleFEEdge2ndInterpolate()
        {

        }

        public TriangleFEEdge2ndInterpolate(TriangleFE owner)
        {
            Owner = owner;
        }


        public uint GetNodeCount()
        {
            // 6つの頂点
            return 6;
        }

        public uint GetEdgeCount()
        {
            // 8つの辺
            return 8;
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
            IList<int> vNos;
            {
                var list = new List<KeyValuePair<int, int>>();
                for (int i = 0; i < 3; i++)
                {
                    list.Add(new KeyValuePair<int, int>(Owner.VertexCoordIds[i], i));
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

                vNos = new List<int>();
                foreach (var pair in list)
                {
                    int coId = pair.Key;
                    int vNo = pair.Value;
                    vNos.Add(vNo);
                }
            }

            int[][] edgePointId = new int[8][]
            {
                new int[]{ 1, 4 },
                new int[]{ 4, 2 },
                new int[]{ 2, 5 },
                new int[]{ 5, 0 },
                new int[]{ 0, 3 },
                new int[]{ 3, 1 },
                null,
                null
            };

            for (int eIndex = 6; eIndex < 8; eIndex++)
            {
                int vNo = vNos[eIndex - 6];
                int[] pointIds = null;
                if (vNo == 0)
                {
                    pointIds = new int[] { 4, 0 };
                }
                else if (vNo == 1)
                {
                    pointIds = new int[] { 5, 1 };
                }
                else if (vNo == 2)
                {
                    pointIds = new int[] { 3, 2 };
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
                edgePointId[eIndex] = pointIds;
            }
            return edgePointId;
        }

        public double[][] CalcEdgeN(double[] L)
        {
            double[][] edgeNs = new double[8][];
            double[] a;
            double[] b;
            double[] c;
            Owner.CalcTransMatrix(out a, out b, out c);

            double A = Owner.GetArea();
            double[] barA = new double[3];
            double[] barB = new double[3];
            double[] barC = new double[3];
            for (int i = 0; i < 3; i++)
            {
                barA[i] = a[i] * (2.0 * A);
                barB[i] = b[i] * (2.0 * A);
                barC[i] = c[i] * (2.0 * A);
            }
            double[] lens = new double[3];
            for (int i = 0; i < 3; i++)
            {
                lens[i] = Math.Sqrt(barB[i] * barB[i] + barC[i] * barC[i]);
            }
            
            double[][] gradLs = new double[3][];
            for (int vIndex = 0; vIndex < 2; vIndex++)
            {
                double[] gradL = new double[2];
                gradLs[vIndex] = gradL;

                gradL[0] = b[vIndex];
                gradL[1] = c[vIndex];
            }
            {
                int vIndex = 2;
                double[] gradL = new double[2];
                gradLs[vIndex] = gradL;

                for (int idim = 0; idim < 2; idim++)
                {
                    double[] gradL0 = gradLs[0];
                    double[] gradL1 = gradLs[1];
                    gradL[idim] = -gradL0[idim] - gradL1[idim];
                }
            }

            int[][] edgePointId = GetEdgePointIdss();

            // t1...t6
            for (int vIndex = 0; vIndex < 3; vIndex++)
            {
                for (int i = 0; i < 2; i++)
                {
                    int eIndex = vIndex * 2 + i;
                    double[] edgeN = new double[2];
                    edgeNs[eIndex] = edgeN;
                    // x,y成分
                    if (i % 2 == 0)
                    {
                        for (int idim = 0; idim < 2; idim++)
                        {
                            edgeN[idim] = lens[vIndex] * L[(vIndex + 1) % 3] * gradLs[(vIndex + 2) % 3][idim]; 
                        }
                    }
                    else
                    {
                        for (int idim = 0; idim < 2; idim++)
                        {
                            edgeN[idim] = -lens[vIndex] * L[(vIndex + 2) % 3] * gradLs[(vIndex + 1) % 3][idim];
                        }
                    }
                }
            }
            
            // t7,t8
            for (int iVNo = 0; iVNo < 2; iVNo++)
            {
                int eIndex = 6 + iVNo;
                int vIndex = edgePointId[eIndex][1];
                System.Diagnostics.Debug.Assert(vIndex >= 0 && vIndex < 3);
                double[] edgeN = new double[2];
                edgeNs[eIndex] = edgeN;

                double[] n = CadUtils2D.Normalize2D(gradLs[vIndex]);

                // x,y成分
                for (int idim = 0; idim < 2; idim++)
                {
                    edgeN[idim] = 4.0 * L[(vIndex + 1) % 3] * L[(vIndex + 2) % 3] * n[idim];
                }
            }

            return edgeNs;
        }

        public double[] CalcRotEdgeN(double[] L)
        {
            double[] rotEdgeNs = new double[8];
            double[] a;
            double[] b;
            double[] c;
            Owner.CalcTransMatrix(out a, out b, out c);

            double A = Owner.GetArea();
            double[] barA = new double[3];
            double[] barB = new double[3];
            double[] barC = new double[3];
            for (int i = 0; i < 3; i++)
            {
                barA[i] = a[i] * (2.0 * A);
                barB[i] = b[i] * (2.0 * A);
                barC[i] = c[i] * (2.0 * A);
            }
            double[] lens = new double[3];
            for (int i = 0; i < 3; i++)
            {
                lens[i] = Math.Sqrt(barB[i] * barB[i] + barC[i] * barC[i]);
            }

            double[][] gradLs = new double[3][];
            for (int vIndex = 0; vIndex < 2; vIndex++)
            {
                double[] gradL = new double[2];
                gradLs[vIndex] = gradL;

                gradL[0] = b[vIndex];
                gradL[1] = c[vIndex];
            }
            {
                int vIndex = 2;
                double[] gradL = new double[2];
                gradLs[vIndex] = gradL;

                for (int idim = 0; idim < 2; idim++)
                {
                    double[] gradL0 = gradLs[0];
                    double[] gradL1 = gradLs[1];
                    gradL[idim] = -gradL0[idim] - gradL1[idim];
                }
            }

            int[][] edgePointId = GetEdgePointIdss();

            double gradLigradLiplus1 = 1.0 / (2.0 * A);
            {
                int eIndex = 0;
                // z成分
                rotEdgeNs[eIndex] = lens[0] * gradLigradLiplus1;
            }
            {
                int eIndex = 1;
                // z成分
                rotEdgeNs[eIndex] = lens[0] * gradLigradLiplus1;
            }
            {
                int eIndex = 2;
                // z成分
                rotEdgeNs[eIndex] = lens[1] * gradLigradLiplus1;
            }
            {
                int eIndex = 3;
                // z成分
                rotEdgeNs[eIndex] = lens[1] * gradLigradLiplus1;
            }
            {
                int eIndex = 4;
                // z成分
                rotEdgeNs[eIndex] = lens[2] * gradLigradLiplus1;
            }
            {
                int eIndex = 5;
                // z成分
                rotEdgeNs[eIndex] = lens[2] * gradLigradLiplus1;
            }
            for (int iVNo = 0; iVNo < 2; iVNo++)
            {
                int eIndex = 6 + iVNo;
                int vIndex = edgePointId[eIndex][1];
                System.Diagnostics.Debug.Assert(vIndex >= 0 && vIndex < 3);
                double value = 0.0;
                double gradLAbs = Math.Sqrt(
                    gradLs[vIndex][0] * gradLs[vIndex][0] +
                    gradLs[vIndex][1] * gradLs[vIndex][1]);
                if (vIndex == 0)
                {
                    value =
                        4.0 * L[2] * (-gradLigradLiplus1) +
                        4.0 * L[1] * gradLigradLiplus1;
                    value /= gradLAbs;
                }
                else if (vIndex == 1)
                {
                    value =
                        4.0 * L[0] * (-gradLigradLiplus1) +
                        4.0 * L[2] * gradLigradLiplus1;
                    value /= gradLAbs;
                }
                else if (vIndex == 2)
                {
                    value =
                        4.0 * L[1] * (-gradLigradLiplus1) +
                        4.0 * L[0] * gradLigradLiplus1;
                    value /= gradLAbs;
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
                rotEdgeNs[eIndex] = value;
            }
            return rotEdgeNs;
        }

        public double[][][] CalcEdgeNu(double[] L)
        {
            double[] a;
            double[] b;
            double[] c;
            Owner.CalcTransMatrix(out a, out b, out c);

            double A = Owner.GetArea();
            double[] barA = new double[3];
            double[] barB = new double[3];
            double[] barC = new double[3];
            for (int i = 0; i < 3; i++)
            {
                barA[i] = a[i] * (2.0 * A);
                barB[i] = b[i] * (2.0 * A);
                barC[i] = c[i] * (2.0 * A);
            }
            double[] lens = new double[3];
            for (int i = 0; i < 3; i++)
            {
                lens[i] = Math.Sqrt(barB[i] * barB[i] + barC[i] * barC[i]);
            }

            double[][] gradLs = new double[3][];
            for (int i = 0; i < 2; i++)
            {
                double[] gradL = new double[2];
                gradLs[i] = gradL;

                gradL[0] = b[i];
                gradL[1] = c[i];
            }
            {
                int i = 2;
                double[] gradL = new double[2];
                gradLs[i] = gradL;

                for (int idim = 0; idim < 2; idim++)
                {
                    double[] gradL0 = gradLs[0];
                    double[] gradL1 = gradLs[1];
                    gradL[idim] = -gradL0[idim] - gradL1[idim];
                }
            }

            int[][] edgePointId = GetEdgePointIdss();

            double[][][] edgeNus = new double[2][][];
            double[][] edgeNxs = new double[8][];
            edgeNus[0] = edgeNxs;
            // t1...t6
            for (int vIndex = 0; vIndex < 3; vIndex++)
            {
                for (int i = 0; i < 2; i++)
                {
                    int eIndex = vIndex * 2 + i;
                    double[] edgeNx = new double[2];
                    edgeNxs[eIndex] = edgeNx;
                    // x,y成分
                    if (i % 2 == 0)
                    {
                        for (int idim = 0; idim < 2; idim++)
                        {
                            edgeNx[idim] = lens[vIndex] * b[(vIndex + 1) % 3] * gradLs[(vIndex + 2) % 3][idim];
                        }
                    }
                    else
                    {
                        for (int idim = 0; idim < 2; idim++)
                        {
                            edgeNx[idim] = -lens[vIndex] * b[(vIndex + 2) % 3] * gradLs[(vIndex + 1) % 3][idim];
                        }
                    }
                }
            }
            // t7,t8
            for (int iVNo = 0; iVNo < 2; iVNo++)
            {
                int eIndex = 6 + iVNo;
                int vIndex = edgePointId[eIndex][1];
                System.Diagnostics.Debug.Assert(vIndex >= 0 && vIndex < 3);
                double[] edgeNx = new double[2];
                edgeNxs[eIndex] = edgeNx;

                double[] n = CadUtils2D.Normalize2D(gradLs[vIndex]);

                // x,y成分
                for (int idim = 0; idim < 2; idim++)
                {
                    edgeNx[idim] = 4.0 * (
                        b[(vIndex + 1) % 3] * L[(vIndex + 2) % 3] +
                        b[(vIndex + 2) % 3] * L[(vIndex + 1) % 3]
                        ) * n[idim];
                }
            }

            double[][] edgeNys = new double[8][];
            edgeNus[1] = edgeNys;
            // t1...t6
            for (int vIndex = 0; vIndex < 3; vIndex++)
            {
                for (int i = 0; i < 2; i++)
                {
                    int eIndex = vIndex * 2 + i;
                    double[] edgeNy = new double[2];
                    edgeNys[eIndex] = edgeNy;
                    // x,y成分
                    if (i % 2 == 0)
                    {
                        for (int idim = 0; idim < 2; idim++)
                        {
                            edgeNy[idim] = lens[vIndex] * c[(vIndex + 1) % 3] * gradLs[(vIndex + 2) % 3][idim];
                        }
                    }
                    else
                    {
                        for (int idim = 0; idim < 2; idim++)
                        {
                            edgeNy[idim] = -lens[vIndex] * c[(vIndex + 2) % 3] * gradLs[(vIndex + 1) % 3][idim];
                        }
                    }
                }
            }
            // t7,t8
            for (int iVNo = 0; iVNo < 2; iVNo++)
            {
                int eIndex = 6 + iVNo;
                int vIndex = edgePointId[eIndex][1];
                System.Diagnostics.Debug.Assert(vIndex >= 0 && vIndex < 3);
                double[] edgeNy = new double[2];
                edgeNys[eIndex] = edgeNy;

                double[] n = CadUtils2D.Normalize2D(gradLs[vIndex]);

                // x,y成分
                for (int idim = 0; idim < 2; idim++)
                {
                    edgeNy[idim] = 4.0 * (
                        c[(vIndex + 1) % 3] * L[(vIndex + 2) % 3] +
                        c[(vIndex + 2) % 3] * L[(vIndex + 1) % 3]
                        ) * n[idim];
                }
            }
            return edgeNus;
        }
    }
}
