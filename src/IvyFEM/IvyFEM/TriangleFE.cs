using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TriangleFE : FE
    {
        public TriangleFE() : base()
        {
            Type = ElementType.Tri;
            Order = 1;
            FEType = FiniteElementType.ScalarLagrange;
            VertexCount = 3;
            CreateInterpolate();
            NodeCount = GetNodeCount();
        }

        public TriangleFE(int order, FiniteElementType feType) : base()
        {
            Type = ElementType.Tri;
            Order = order;
            FEType = feType;
            VertexCount = 3;
            CreateInterpolate();
            NodeCount = GetNodeCount();
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
            OpenTK.Vector2d v1 = new OpenTK.Vector2d(co1[0], co1[1]);
            OpenTK.Vector2d v2 = new OpenTK.Vector2d(co2[0], co2[1]);
            OpenTK.Vector2d v3 = new OpenTK.Vector2d(co3[0], co3[1]);
            double area = CadUtils.TriArea(v1, v2, v3);
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
            double[] co1 = World.GetVertexCoord(VertexCoordIds[0]);
            double[] co2 = World.GetVertexCoord(VertexCoordIds[1]);
            double[] co3 = World.GetVertexCoord(VertexCoordIds[2]);
            OpenTK.Vector2d v1 = new OpenTK.Vector2d(co1[0], co1[1]);
            OpenTK.Vector2d v2 = new OpenTK.Vector2d(co2[0], co2[1]);
            OpenTK.Vector2d v3 = new OpenTK.Vector2d(co3[0], co3[1]);
            double A = CadUtils.TriArea(v1, v2, v3);
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

            double x = pt[0];
            double y = pt[1];
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
    }
}
