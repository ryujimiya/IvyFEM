using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class LineFE : FE
    {
        public LineFE() : base()
        {
            Type = ElementType.Line;
            Order = 1;
            VertexCount = 2;
            NodeCount = GetNodeCount();
        }

        public LineFE(int order) : base()
        {
            Type = ElementType.Line;
            Order = order;
            VertexCount = 2;
            NodeCount = GetNodeCount();
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
            uint nodeCnt = 0;
            if (Order == 1)
            {
                nodeCnt = 2;
            }
            else if (Order == 2)
            {
                nodeCnt = 3;
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return nodeCnt;
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
            /*
            OpenTK.Vector2d v1 = new OpenTK.Vector2d(co1[0], co1[1]);
            OpenTK.Vector2d v2 = new OpenTK.Vector2d(co2[0], co2[1]);
            var t = v2 - v1;
            t = CadUtils.Normalize(t);
            // n = t x e3
            double[] normal = { t.Y, -t.X};
            */
            double[] normal = IvyFEM.CadUtils.GetNormal2D(co1, co2);
            return normal;
        }

        public double[] GetNodeL(int nodeId)
        {
            double[] ret = null;
            if (Order == 1)
            {
                ret = Get1stNodeL(nodeId);
            }
            else if (Order == 2)
            {
                ret = Get2ndNodeL(nodeId);
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
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
            double[] ret = null;
            if (Order == 1)
            {
                ret = Calc1stN(L);
            }
            else if (Order == 2)
            {
                ret = Calc2ndN(L);
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return ret;
        }

        /// <summary>
        /// dN/du
        /// </summary>
        /// <returns></returns>
        public double[][] CalcNu(double[] L)
        {
            double[][] ret = null;
            if (Order == 1)
            {
                ret = Calc1stNu(L);
            }
            else if (Order == 2)
            {
                ret = Calc2ndNu(L);
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return ret;
        }

        /// <summary>
        /// S{N}Tdx
        /// </summary>
        /// <returns></returns>
        public double[] CalcSN()
        {
            double[] ret = null;
            if (Order == 1)
            {
                ret = Calc1stSN();
            }
            else if (Order == 2)
            {
                ret = Calc2ndSN();
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return ret;
        }

        /// <summary>
        /// S{N}{N}Tdx
        /// </summary>
        /// <returns></returns>
        public double[,] CalcSNN()
        {
            double[,] ret = null;
            if (Order == 1)
            {
                ret = Calc1stSNN();
            }
            else if (Order == 2)
            {
                ret = Calc2ndSNN();
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return ret;
        }

        /// <summary>
        /// S{Nx}{Nx}Tdx
        /// </summary>
        /// <returns></returns>
        public double[,] CalcSNxNx()
        {
            double[,] ret = null;
            if (Order == 1)
            {
                ret = Calc1stSNxNx();
            }
            else if (Order == 2)
            {
                ret = Calc2ndSNxNx();
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return ret;
        }

    }
}
