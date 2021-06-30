using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TetrahedronFEEdge1stInterpolate : IEdgeInterpolate3D
    {
        public TetrahedronFE Owner { get; set; }

        public TetrahedronFEEdge1stInterpolate()
        {

        }

        public TetrahedronFEEdge1stInterpolate(TetrahedronFE owner)
        {
            Owner = owner;
        }


        public uint GetNodeCount()
        {
            // 4つの頂点
            return 4;
        }

        public uint GetEdgeCount()
        {
            // 6つの辺
            return 6;
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

        public double[][] CalcEdgeN(double[] L)
        {
            double[][] edgeNs = new double[6][];
            double[] a;
            double[] b;
            double[] c;
            double[] d;
            Owner.CalcTransMatrix(out a, out b, out c, out d);

            int[][] edgePointId = GetEdgePointIdss();

            double[] lens = new double[6];
            for (int eIndex = 0; eIndex < 6; eIndex++)
            {
                int[] edgePointId1 = edgePointId[eIndex];
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

            double[][] gradLs = new double[4][];
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

            for (int eIndex = 0; eIndex < 6; eIndex++)
            {
                double[] edgeN = new double[3];
                edgeNs[eIndex] = edgeN;

                int[] edgePointId1 = edgePointId[eIndex];
                int iPtId = edgePointId1[0];
                int jPtId = edgePointId1[1];
                for (int idim = 0; idim < 3; idim++)
                {
                    // x,y,z成分
                    edgeN[idim] = lens[eIndex] * (
                        L[iPtId] * gradLs[jPtId][idim] -
                        L[jPtId] * gradLs[iPtId][idim]);
                }
            }
            return edgeNs;
        }

        public double[][] CalcRotEdgeN(double[] L)
        {
            double[][] rotEdgeNs = new double[6][];
            double[] a;
            double[] b;
            double[] c;
            double[] d;
            Owner.CalcTransMatrix(out a, out b, out c, out d);

            int[][] edgePointId = GetEdgePointIdss();

            double[] lens = new double[6];
            for (int eIndex = 0; eIndex < 6; eIndex++)
            {
                int[] edgePointId1 = edgePointId[eIndex];
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

            double[][] gradLs = new double[4][];
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

            for (int eIndex = 0; eIndex < 6; eIndex++)
            {
                int[] edgePointId1 = edgePointId[eIndex];
                int iPtId = edgePointId1[0];
                int jPtId = edgePointId1[1];
                OpenTK.Vector3d gradLVec1 = new OpenTK.Vector3d(gradLs[iPtId][0], gradLs[iPtId][1], gradLs[iPtId][2]);
                OpenTK.Vector3d gradLVec2 = new OpenTK.Vector3d(gradLs[jPtId][0], gradLs[jPtId][1], gradLs[jPtId][2]);
                var gradL1xgrapdL2Vec = OpenTK.Vector3d.Cross(gradLVec1, gradLVec2);
                double[] gradL1xgradL2 = { gradL1xgrapdL2Vec.X, gradL1xgrapdL2Vec.Y, gradL1xgrapdL2Vec.Z };

                double[] rotEdgeN = new double[3];
                rotEdgeNs[eIndex] = rotEdgeN;
                for (int idim = 0; idim < 3; idim++)
                {
                    rotEdgeN[idim] = 2.0 * lens[eIndex] * gradL1xgradL2[idim];
                }
            }
            return rotEdgeNs;
        }
    }
}
