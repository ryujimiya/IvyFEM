using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class LineFE : FE
    {
        public LineFE() : base()
        {
            Type = ElementType.Line;
            Order = 1;
            VertexCount = 2;
            CreateInterpolate();
            NodeCount = GetNodeCount();
            EdgeCount = GetEdgeCount();
        }

        public LineFE(int order, FiniteElementType feType) : base()
        {
            Type = ElementType.Line;
            Order = order;
            FEType = feType;
            VertexCount = 2;
            CreateInterpolate();
            NodeCount = GetNodeCount();
            EdgeCount = GetEdgeCount();
        }

        public LineFE(LineFE src)
        {
            Copy(src);
        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }

        protected uint GetNodeCount()
        {
            return Interpolate.GetNodeCount();
        }

        protected uint GetEdgeCount()
        {
            if (!(Interpolate is IEdgeInterpolate))
            {
                return 0;
            }
            IEdgeInterpolate edgeInterpolate = Interpolate as IEdgeInterpolate;
            return edgeInterpolate.GetEdgeCount();
        }

        private double[] AddDisplacement(int iNode, double[] co)
        {
            int dim = co.Length;
            double[] curCo = new double[dim];
            co.CopyTo(curCo, 0);
            if (Displacements != null)
            {
                double[] u = Displacements[iNode];
                for (int iDim = 0; iDim < dim; iDim++)
                {
                    curCo[iDim] += u[iDim];
                }
            }
            return curCo;
        }

        public double GetLineLength()
        {
            double[] co1 = World.GetVertexCoord(VertexCoordIds[0]);
            double[] co2 = World.GetVertexCoord(VertexCoordIds[1]);
            co1 = AddDisplacement(0, co1);
            co2 = AddDisplacement(1, co2);
            OpenTK.Vector2d v1 = new OpenTK.Vector2d(co1[0], co1[1]);
            OpenTK.Vector2d v2 = new OpenTK.Vector2d(co2[0], co2[1]);
            double l = (v2 - v1).Length;
            return l;
        }

        public double[] GetNormal()
        {
            double[] co1 = World.GetVertexCoord(VertexCoordIds[0]);
            double[] co2 = World.GetVertexCoord(VertexCoordIds[1]);
            co1 = AddDisplacement(0, co1);
            co2 = AddDisplacement(1, co2);
            double[] normal = IvyFEM.CadUtils.GetNormal2D(co1, co2);
            return normal;
        }

        public double[] GetNodeL(int nodeId)
        {
            double[] ret = Interpolate.GetNodeL(nodeId);
            return ret;
        }

        public void CalcTransMatrix(out double[] a, out double[] b)
        {
            a = new double[2];
            b = new double[2];
            double[] co1 = World.GetVertexCoord(VertexCoordIds[0]);
            double[] co2 = World.GetVertexCoord(VertexCoordIds[1]);
            co1 = AddDisplacement(0, co1);
            co2 = AddDisplacement(1, co2);
            OpenTK.Vector2d v1 = new OpenTK.Vector2d(co1[0], co1[1]);
            OpenTK.Vector2d v2 = new OpenTK.Vector2d(co2[0], co2[1]);
            var dir = v2 - v1;
            dir = CadUtils.Normalize(dir);
            double l = GetLineLength();
            {
                a[0] = (1.0 / l) * OpenTK.Vector2d.Dot(v2, dir);
                a[1] = (1.0 / l) * OpenTK.Vector2d.Dot(-v1, dir);
            }
            {
                b[0] = (1.0 / l) * (-1.0);
                b[1] = (1.0 / l) * 1.0;
            }
        }

        public double[] L2Coord(double[] L)
        {
            double[] pt = new double[2];
            double[][] ptValue = new double[NodeCount][];
            int[] coIds = NodeCoordIds;
            for (int iNode = 0; iNode < NodeCount; iNode++)
            {
                int coId = coIds[iNode];
                ptValue[iNode] = World.GetCoord((uint)QuantityId, coId);
            }

            double[] N = CalcN(L);
            for (int iNode = 0; iNode < NodeCount; iNode++)
            {
                for (int iDof = 0; iDof < 2; iDof++)
                {
                    pt[iDof] += N[iNode] * ptValue[iNode][iDof];
                }
            }
            return pt;
        }

        public double[] Coord2L(double[] pt)
        {
            System.Diagnostics.Debug.Assert(pt.Length == 2);
            double[] co1 = World.GetVertexCoord(VertexCoordIds[0]);
            double[] co2 = World.GetVertexCoord(VertexCoordIds[1]);
            OpenTK.Vector2d v1 = new OpenTK.Vector2d(co1[0], co1[1]);
            OpenTK.Vector2d v2 = new OpenTK.Vector2d(co2[0], co2[1]);
            OpenTK.Vector2d v = new OpenTK.Vector2d(pt[0], pt[1]);
            var dir12 = v2 - v1;
            dir12 = CadUtils.Normalize(dir12);
            var dir = v - v1;
            dir = CadUtils.Normalize(dir);
            double x = dir12[0] * dir[0] + dir12[1] * dir[1];

            double[] a;
            double[] b;
            CalcTransMatrix(out a, out b);
            double L1 = a[0] + b[0] * x;
            double L2 = a[1] + b[1] * x;
            double[] L = { L1, L2 };
            return L;
        }

        public double[] GetEdgeL(int edgeId)
        {
            IEdgeInterpolate edgeInterpolate = Interpolate as IEdgeInterpolate;
            double[] ret = edgeInterpolate.GetEdgeL(edgeId);
            return ret;
        }

        // ξ([-1, 1])からL1,L2に変換
        public static double[] GetLFromXi(double xi)
        {
            return new double[2] { (1.0 - xi) / 2.0, (1.0 + xi) / 2.0 };
        }

        public static IntegrationPoints GetIntegrationPoints(LineIntegrationPointCount integrationPointCount)
        {
            foreach (var ip in IntegrationPoints.LineIntegrationPoints)
            {
                if (ip.PointCount == (int)integrationPointCount)
                {
                    return ip;
                }
            }
            System.Diagnostics.Debug.Assert(false);
            return null;
        }

        /// <summary>
        /// N
        /// </summary>
        /// <param name="L"></param>
        /// <returns></returns>
        public double[] CalcN(double[] L)
        {
            double[] ret = Interpolate.CalcN(L);
            return ret;
        }

        /// <summary>
        /// dN/du
        /// </summary>
        /// <returns></returns>
        public double[][] CalcNu(double[] L)
        {
            double[][] ret = Interpolate.CalcNu(L);
            return ret;
        }

        /// <summary>
        /// S{N}Tdx
        /// </summary>
        /// <returns></returns>
        public double[] CalcSN()
        {
            double[] ret = Interpolate.CalcSN();
            return ret;
        }

        /// <summary>
        /// S{N}{N}Tdx
        /// </summary>
        /// <returns></returns>
        public double[,] CalcSNN()
        {
            double[,] ret = Interpolate.CalcSNN();
            return ret;
        }

        /// <summary>
        /// S{Nx}{Nx}Tdx
        /// </summary>
        /// <returns></returns>
        public double[,] CalcSNxNx()
        {
            double[,][,] sNuNv = Interpolate.CalcSNuNv();

            double[,] ret = sNuNv[0, 0];
            return ret;
        }

        /// <summary>
        /// vecN
        /// </summary>
        /// <param name="L"></param>
        /// <returns></returns>
        public double[][] CalcEdgeN(double[] L)
        {
            IEdgeInterpolate edgeInterpolate = Interpolate as IEdgeInterpolate;
            double[][] ret = edgeInterpolate.CalcEdgeN(L);
            return ret;
        }

        /// <summary>
        /// [rot(vecN)]z
        /// </summary>
        /// <param name="L"></param>
        /// <returns></returns>
        public double[] CalcRotEdgeN(double[] L)
        {
            IEdgeInterpolate edgeInterpolate = Interpolate as IEdgeInterpolate;
            double[] ret = edgeInterpolate.CalcRotEdgeN(L);
            return ret;
        }

        /// <summary>
        /// d(vecN)/du (u=x,y)
        /// </summary>
        /// <param name="L"></param>
        /// <returns></returns>
        public double[][][] CalcEdgeNu(double[] L)
        {
            IEdgeInterpolate edgeInterpolate = Interpolate as IEdgeInterpolate;
            double[][][] ret = edgeInterpolate.CalcEdgeNu(L);
            return ret;
        }
    }
}
