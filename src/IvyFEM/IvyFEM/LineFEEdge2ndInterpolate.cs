using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class LineFEEdge2ndInterpolate : IEdgeInterpolate
    {
        public LineFE Owner { get; set; }

        public  LineFEEdge2ndInterpolate()
        {

        }

        public LineFEEdge2ndInterpolate(LineFE owner)
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
            // 2つの辺
            return 2;
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
            int[][] edgePointId = new int[2][]
            {
                new int[]{ 0, 2 },
                new int[]{ 2, 1 }
            };
            return edgePointId;
        }

        public double[][] CalcEdgeN(double[] L)
        {
            double[][] edgeNs = new double[2][];
            double l = Owner.GetLineLength();
            double[] normal = Owner.GetNormal();
            double[] tan = { normal[1], -normal[0] };
            double[] dir = { -tan[0], -tan[1] };

            for (int eIndex = 0; eIndex < 2; eIndex++)
            {
                double[] edgeN = new double[2];
                edgeNs[eIndex] = edgeN;

                for (int idim = 0; idim < 2; idim++)
                {
                    // x,y成分
                    edgeN[idim] = L[eIndex] * dir[idim];
                }
            }
            return edgeNs;
        }

        public double[] CalcRotEdgeN(double[] L)
        {
            throw new NotImplementedException();
        }

        public double[][][] CalcEdgeNu(double[] L)
        {
            throw new NotImplementedException();
        }
    }
}
