using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class LineFEEdge1stInterpolate : IEdgeInterpolate
    {
        public LineFE Owner { get; set; }

        public  LineFEEdge1stInterpolate()
        {

        }

        public LineFEEdge1stInterpolate(LineFE owner)
        {
            Owner = owner;
        }


        public uint GetNodeCount()
        {
            // 2つの頂点 
            return 2;
        }

        public uint GetEdgeCount()
        {
            // 1つの辺
            return 1;
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
            double[][] edgeL = new double[1][]
            {
                new double[] { 0.5, 0.5 } // 1-2の中点
            };
            return edgeL[edgeId];
        }

        public int[][] GetEdgePointIdss()
        {
            int[][] edgePointId = new int[1][]
            {
                new int[]{ 0, 1 }
            };
            return edgePointId;
        }

        public double[][] CalcEdgeN(double[] L)
        {
            double[][] edgeNs = new double[1][];
            double[] a;
            double[] b;
            Owner.CalcTransMatrix(out a, out b);

            double l = Owner.GetLineLength();
            double[] normal = Owner.GetNormal();
            double[] tan = { normal[1], -normal[0] };
            double[] dir = { -tan[0], -tan[1] };

            for (int eIndex = 0; eIndex < 1; eIndex++)
            {
                double[] edgeN = new double[2];
                edgeNs[eIndex] = edgeN;

                for (int idim = 0; idim < 2; idim++)
                {
                    // x,y成分
                    edgeN[idim] = tan[idim];
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
