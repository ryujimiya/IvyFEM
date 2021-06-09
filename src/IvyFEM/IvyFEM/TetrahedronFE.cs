using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TetrahedronFE : FE
    {
        public TetrahedronFE() : base()
        {
            Type = ElementType.Tet;
            Order = 1;
            FEType = FiniteElementType.ScalarLagrange;
            VertexCount = 4;
            CreateInterpolate();
            NodeCount = GetNodeCount();
            EdgeCount = GetEdgeCount();
        }

        public TetrahedronFE(int order, FiniteElementType feType) : base()
        {
            Type = ElementType.Tet;
            Order = order;
            FEType = feType;
            VertexCount = 4;
            CreateInterpolate();
            NodeCount = GetNodeCount();
            EdgeCount = GetEdgeCount();
        }

        public TetrahedronFE(TriangleFE src)
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
            if (!(Interpolate is IEdgeInterpolate3D))
            {
                return 0;
            }
            IEdgeInterpolate3D edgeInterpolate = Interpolate as IEdgeInterpolate3D;
            return edgeInterpolate.GetEdgeCount();
        }

        public double GetVolume()
        {
            double[] co1 = World.GetVertexCoord(VertexCoordIds[0]);
            double[] co2 = World.GetVertexCoord(VertexCoordIds[1]);
            double[] co3 = World.GetVertexCoord(VertexCoordIds[2]);
            double[] co4 = World.GetVertexCoord(VertexCoordIds[3]);
            co1 = AddDisplacement(0, co1);
            co2 = AddDisplacement(1, co2);
            co3 = AddDisplacement(2, co3);
            co4 = AddDisplacement(3, co4);
            double vol = CadUtils.TetVolume(co1, co2, co3, co4);
            return vol;
        }

        public double[] GetEdgeLengths()
        {
            double[] edgeLengths = new double[6];
            int[][] edgePointId = new int[6][]
            {
                new int[]{ 0, 1 },
                new int[]{ 1, 2 },
                new int[]{ 0, 2 },
                new int[]{ 0, 3 },
                new int[]{ 3, 1 },
                new int[]{ 2, 3 }
            };

            for (int eIndex = 0; eIndex < 6; eIndex++)
            {
                int[] edgePointId1 = edgePointId[eIndex];
                int iPtId = edgePointId1[0];
                int jPtId = edgePointId1[1];
                int coId1 = NodeCoordIds[iPtId];
                int coId2 = NodeCoordIds[jPtId];
                double[] co1 = World.GetCoord((uint)QuantityId, coId1);
                double[] co2 = World.GetCoord((uint)QuantityId, coId2);
                edgeLengths[eIndex] = Math.Sqrt(
                    (co2[0] - co1[0]) * (co2[0] - co1[0]) +
                    (co2[1] - co1[1]) * (co2[1] - co1[1]) +
                    (co2[2] - co1[2]) * (co2[2] - co1[2]));
            }
            return edgeLengths;
        }

        public void CalcTransMatrix(out double[] a, out double[] b, out double[] c, out double[] d)
        {
            a = new double[4];
            b = new double[4];
            c = new double[4];
            d = new double[4];
            double[] co1 = World.GetVertexCoord(VertexCoordIds[0]);
            double[] co2 = World.GetVertexCoord(VertexCoordIds[1]);
            double[] co3 = World.GetVertexCoord(VertexCoordIds[2]);
            double[] co4 = World.GetVertexCoord(VertexCoordIds[3]);
            co1 = AddDisplacement(0, co1);
            co2 = AddDisplacement(1, co2);
            co3 = AddDisplacement(2, co3);
            co4 = AddDisplacement(3, co4);
            OpenTK.Vector3d v1 = new OpenTK.Vector3d(co1[0], co1[1], co1[2]);
            OpenTK.Vector3d v2 = new OpenTK.Vector3d(co2[0], co2[1], co2[2]);
            OpenTK.Vector3d v3 = new OpenTK.Vector3d(co3[0], co3[1], co3[2]);
            OpenTK.Vector3d v4 = new OpenTK.Vector3d(co4[0], co4[1], co4[2]);
            double V = CadUtils3D.TetVolume(v1, v2, v3, v4);
            OpenTK.Vector3d[] p = { v1, v2, v3, v4 };
            int[] e = { 1, -1, 1, -1 };
            for (int k = 0; k < 4; k++)
            {
                int l = (k + 1) % 4;
                int m = (k + 2) % 4;
                int n = (k + 3) % 4;
                a[k] = (1.0 / (6.0 * V)) * e[k] * (
                    p[l].X * (p[m].Y * p[n].Z - p[n].Y * p[m].Z) +
                    p[m].X * (p[n].Y * p[l].Z - p[l].Y * p[n].Z) +
                    p[n].X * (p[l].Y * p[m].Z - p[m].Y * p[l].Z));
                b[k] = (1.0 / (6.0 * V)) * e[k] * (
                    p[l].Y * (p[n].Z - p[m].Z) +
                    p[m].Y * (p[l].Z - p[n].Z) +
                    p[n].Y * (p[m].Z - p[l].Z));
                c[k] = (1.0 / (6.0 * V)) * e[k] * (
                    p[l].Z * (p[n].X - p[m].X) +
                    p[m].Z * (p[l].X - p[n].X) +
                    p[n].Z * (p[m].X - p[l].X));
                d[k] = (1.0 / (6.0 * V)) * e[k] * (
                    p[l].X * (p[n].Y - p[m].Y) +
                    p[m].X * (p[l].Y - p[n].Y) +
                    p[n].X * (p[m].Y - p[l].Y));
            }
        }

        public double[] GetNodeL(int nodeId)
        {
            double[] ret = Interpolate.GetNodeL(nodeId);
            return ret;
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
                ptValue[iNode] = AddDisplacement(iNode, ptValue[iNode]);
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
            double x = pt[0];
            double y = pt[1];
            double z = pt[2];
            double[] a;
            double[] b;
            double[] c;
            double[] d;
            CalcTransMatrix(out a, out b, out c, out d);
            double L1 = a[0] + b[0] * x + c[0] * y + d[0] * z;
            double L2 = a[1] + b[1] * x + c[1] * y + d[1] * z;
            double L3 = a[2] + b[2] * x + c[2] * y + d[2] * z;
            double L4 = a[3] + b[3] * x + c[3] * y + d[3] * z;
            double[] L = { L1, L2, L3, L4 };
            return L;
        }

        public double[] GetEdgeL(int edgeId)
        {
            IEdgeInterpolate3D edgeInterpolate = Interpolate as IEdgeInterpolate3D;
            double[] ret = edgeInterpolate.GetEdgeL(edgeId);
            return ret;
        }

        public static IntegrationPoints GetIntegrationPoints(TetrahedronIntegrationPointCount integrationPointCount)
        {
            foreach (var ip in IntegrationPoints.TetrahedronIntegrationPoints)
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
            double V = GetVolume();
            return 6.0 * V;
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
            IEdgeInterpolate3D edgeInterpolate = Interpolate as IEdgeInterpolate3D;
            double[][] ret = edgeInterpolate.CalcEdgeN(L);
            return ret;
        }

        /// <summary>
        /// [rot(vecN)]z
        /// </summary>
        /// <param name="L"></param>
        /// <returns></returns>
        public double[][] CalcRotEdgeN(double[] L)
        {
            IEdgeInterpolate3D edgeInterpolate = Interpolate as IEdgeInterpolate3D;
            double[][] ret = edgeInterpolate.CalcRotEdgeN(L);
            return ret;
        }
    }
}
