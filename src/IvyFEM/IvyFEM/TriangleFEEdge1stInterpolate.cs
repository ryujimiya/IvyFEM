using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TriangleFEEdge1stInterpolate : IEdgeInterpolate
    {
        public TriangleFE Owner { get; set; }

        public TriangleFEEdge1stInterpolate()
        {

        }

        public TriangleFEEdge1stInterpolate(TriangleFE owner)
        {
            Owner = owner;
        }


        public uint GetNodeCount()
        {
            // 3つの頂点
            return 3;
        }

        public uint GetEdgeCount()
        {
            // 3つの辺
            return 3;
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

        public double[] GetEdgeL(int edgeId)
        {
            double[][] edgeL = new double[3][]
            {
                new double[] { 0.0, 0.5, 0.5 }, // 2-3の中点
                new double[] { 0.5, 0.0, 0.5 }, // 3-1の中点
                new double[] { 0.5, 0.5, 0.0 }  // 1-2の中点
            };
            return edgeL[edgeId];
        }

        public int[][] GetEdgePointIdss()
        {
            int[][] edgePointId = new int[3][]
            {
                new int[]{ 1, 2 },
                new int[]{ 2, 0 },
                new int[]{ 0, 1 }
            };
            return edgePointId;
        }

        public double[][] CalcEdgeN(double[] L)
        {
            double[][] edgeNs = new double[3][];
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
            for (int eIndex = 0; eIndex < 2; eIndex++)
            {
                double[] gradL = new double[2];
                gradLs[eIndex] = gradL;

                gradL[0] = b[eIndex];
                gradL[1] = c[eIndex];
            }
            {
                int eIndex = 2;
                double[] gradL = new double[2];
                gradLs[eIndex] = gradL;

                for (int idim = 0; idim < 2; idim++)
                {
                    double[] gradL0 = gradLs[0];
                    double[] gradL1 = gradLs[1];
                    gradL[idim] = -gradL0[idim] - gradL1[idim];
                }
            }

            for (int eIndex = 0; eIndex < 3; eIndex++)
            {
                double[] edgeN = new double[2];
                edgeNs[eIndex] = edgeN;

                for (int idim = 0; idim < 2; idim++)
                {
                    // x,y成分
                    edgeN[idim] = lens[eIndex] * (
                        L[(eIndex + 1) % 3] * gradLs[(eIndex + 2) % 3][idim] -
                        L[(eIndex + 2) % 3] * gradLs[(eIndex + 1) % 3][idim]);
                }
            }
            return edgeNs;
        }

        public double[] CalcRotEdgeN(double[] L)
        {
            double[] rotEdgeNs = new double[3];
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

            for (int eIndex = 0; eIndex < 3; eIndex++)
            {
                // z成分
                rotEdgeNs[eIndex] = lens[eIndex] / A;
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
            for (int eIndex = 0; eIndex < 2; eIndex++)
            {
                double[] gradL = new double[2];
                gradLs[eIndex] = gradL;

                gradL[0] = b[eIndex];
                gradL[1] = c[eIndex];
            }
            {
                int eIndex = 2;
                double[] gradL = new double[2];
                gradLs[eIndex] = gradL;

                for (int idim = 0; idim < 2; idim++)
                {
                    double[] gradL0 = gradLs[0];
                    double[] gradL1 = gradLs[1];
                    gradL[idim] = -gradL0[idim] - gradL1[idim];
                }
            }

            double[][][] edgeNus = new double[2][][];
            double[][] edgeNxs = new double[3][];
            edgeNus[0] = edgeNxs;
            for (int eIndex = 0; eIndex < 3; eIndex++)
            {
                double[] edgeNx = new double[2];
                edgeNxs[eIndex] = edgeNx;
                for (int idim = 0; idim < 2; idim++)
                {
                    edgeNx[idim] = lens[eIndex] * (
                        b[(eIndex + 1) % 3] * gradLs[(eIndex + 2) % 3][idim] -
                        b[(eIndex + 2) % 3] * gradLs[(eIndex + 1) % 3][idim]);
                }
            }
            double[][] edgeNys = new double[3][];
            edgeNus[1] = edgeNys;
            for (int eIndex = 0; eIndex < 3; eIndex++)
            {
                double[] edgeNy = new double[2];
                edgeNys[eIndex] = edgeNy;
                for (int idim = 0; idim < 2; idim++)
                {
                    edgeNy[idim] = lens[eIndex] * (
                        c[(eIndex + 1) % 3] * gradLs[(eIndex + 2) % 3][idim] -
                        c[(eIndex + 2) % 3] * gradLs[(eIndex + 1) % 3][idim]);
                }
            }
            return edgeNus;
        }
    }
}
