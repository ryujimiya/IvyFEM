using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class TriangleFE : FE
    {
        public TriangleFE() : base()
        {
            Type = ElementType.Tri;
            Order = 1;
            FEType = FiniteElementType.ScalarLagrange;
            VertexCount = 3;
            CreateInterpolate();
            NodeCount = GetNodeCount();
            EdgeCount = GetEdgeCount();
        }

        public TriangleFE(int order, FiniteElementType feType) : base()
        {
            Type = ElementType.Tri;
            Order = order;
            FEType = feType;
            VertexCount = 3;
            CreateInterpolate();
            NodeCount = GetNodeCount();
            EdgeCount = GetEdgeCount();
        }

        public TriangleFE(TriangleFE src)
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

        public double[] GetNodeL(int nodeId)
        {
            double[] ret = Interpolate.GetNodeL(nodeId);
            return ret;
        }

        public double GetArea()
        {
            double[] co1 = World.GetVertexCoord(VertexCoordIds[0]);
            double[] co2 = World.GetVertexCoord(VertexCoordIds[1]);
            double[] co3 = World.GetVertexCoord(VertexCoordIds[2]);
            double area = CadUtils.TriArea(co1, co2, co3);
            return area;
        }

        public double[] GetEdgeLengths()
        {
            double[] a;
            double[] b;
            double[] c;
            CalcTransMatrix(out a, out b, out c);
            // 1/(2A)の係数をなくしたa, b, cを求める
            double A = GetArea();
            for (int i = 0; i < 3; i++)
            {
                a[i] *= (2.0 * A);
                b[i] *= (2.0 * A);
                c[i] *= (2.0 * A);
            }

            double[] edgeLengths = new double[3];
            for (int i = 0; i < 3; i++)
            {
                edgeLengths[i] = Math.Sqrt(b[i] * b[i] + c[i] * c[i]);
            }
            return edgeLengths;
        }

        public void CalcTransMatrix(out double[] a, out double[] b, out double[] c)
        {
            a = new double[3];
            b = new double[3];
            c = new double[3];
            double[] co2D1;
            double[] co2D2;
            double[] co2D3;
            if (World.Dimension == 2)
            {
                co2D1 = World.GetVertexCoord(VertexCoordIds[0]);
                co2D2 = World.GetVertexCoord(VertexCoordIds[1]);
                co2D3 = World.GetVertexCoord(VertexCoordIds[2]);
            }
            else if (World.Dimension == 3)
            {
                OpenTK.Vector2d[] projected = ProjectVertexsFrom3D();
                co2D1 = new double[2] { projected[0].X, projected[0].Y };
                co2D2 = new double[2] { projected[1].X, projected[1].Y };
                co2D3 = new double[2] { projected[2].X, projected[2].Y };
            }
            else
            {
                throw new NotImplementedException();
            }
            OpenTK.Vector2d v1 = new OpenTK.Vector2d(co2D1[0], co2D1[1]);
            OpenTK.Vector2d v2 = new OpenTK.Vector2d(co2D2[0], co2D2[1]);
            OpenTK.Vector2d v3 = new OpenTK.Vector2d(co2D3[0], co2D3[1]);
            double A = CadUtils2D.TriArea(v1, v2, v3);
            OpenTK.Vector2d[] v = { v1, v2, v3 };
            for (int k = 0; k < 3; k++)
            {
                int l = (k + 1) % 3;
                int m = (k + 2) % 3;
                a[k] = (1.0 / (2.0 * A)) * (v[l].X * v[m].Y - v[m].X * v[l].Y);
                b[k] = (1.0 / (2.0 * A)) * (v[l].Y - v[m].Y);
                c[k] = (1.0 / (2.0 * A)) * (v[m].X - v[l].X);
            }            
        }

        public double[] L2Coord(double[] L)
        {
            uint dim = World.Dimension;
            double[] pt = new double[dim];
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
                for (int iDof = 0; iDof < dim; iDof++)
                {
                    pt[iDof] += N[iNode] * ptValue[iNode][iDof];
                }
            }
            return pt;
        }

        public double[] Coord2L(double[] pt)
        {
            double[] pt2D;
            if (World.Dimension == 2)
            {
                pt2D = pt; 
            }
            else if (World.Dimension == 3)
            {
                OpenTK.Vector2d projected = ProjectFrom3D(new OpenTK.Vector3d(pt[0], pt[1], pt[2]));
                pt2D = new double[] { projected.X, projected.Y };
            }
            else
            {
                throw new NotImplementedException();
            }

            double x = pt2D[0];
            double y = pt2D[1];
            double[] a;
            double[] b;
            double[] c;
            CalcTransMatrix(out a, out b, out c);
            double L1 = a[0] + b[0] * x + c[0] * y;
            double L2 = a[1] + b[1] * x + c[1] * y;
            double L3 = a[2] + b[2] * x + c[2] * y;
            double[] L = { L1, L2, L3 };
            return L;
        }

        public double[] GetEdgeL(int edgeId)
        {
            IEdgeInterpolate edgeInterpolate = Interpolate as IEdgeInterpolate;
            double[] ret = edgeInterpolate.GetEdgeL(edgeId);
            return ret;
        }

        public static IntegrationPoints GetIntegrationPoints(TriangleIntegrationPointCount integrationPointCount)
        {
            foreach (var ip in IntegrationPoints.TriangleIntegrationPoints)
            {
                if (ip.PointCount == (int)integrationPointCount)
                {
                    return ip;
                }
            }
            System.Diagnostics.Debug.Assert(false);
            return null;
        }

        public double GetDetJacobian(double[] L)
        {
            double A = GetArea();
            return 2.0 * A;
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
        /// d^2 N/(dudv)
        /// </summary>
        /// <param name="L"></param>
        /// <returns></returns>
        public double[,][] CalcNuv(double[] L)
        {
            double[,][] ret = Interpolate.CalcNuv(L);
            return ret;
        }

        /// <summary>
        /// SNdx
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
        /// S{Nu}{Nv}Tdx, u,v = x, y
        /// </summary>
        /// <returns></returns>
        public double[,][,] CalcSNuNv()
        {
            double[,][,] ret = Interpolate.CalcSNuNv();
            return ret;
        }

        /// <summary>
        /// S{Nu}{N}Tdx, u = x, y
        /// </summary>
        /// <returns></returns>
        public double[][,] CalcSNuN()
        {
            double[][,] ret = Interpolate.CalcSNuN();
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
