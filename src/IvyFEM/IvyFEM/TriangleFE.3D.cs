using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    partial class TriangleFE
    {
        public OpenTK.Vector2d ProjectFrom3D(OpenTK.Vector3d p)
        {
            System.Diagnostics.Debug.Assert(World.Dimension == 3);

            OpenTK.Vector3d origin;
            OpenTK.Vector3d normal;
            OpenTK.Vector3d xdir;
            if (FixedOrigin3D != null &&
                FixedNormal3D != null &&
                FixedXDir3D != null)
            {
                origin = new OpenTK.Vector3d(FixedOrigin3D[0], FixedOrigin3D[1], FixedOrigin3D[2]);
                normal = new OpenTK.Vector3d(FixedNormal3D[0], FixedNormal3D[1], FixedNormal3D[2]);
                xdir = new OpenTK.Vector3d(FixedXDir3D[0], FixedXDir3D[1], FixedXDir3D[2]);
            }
            else
            {
                double[] co1 = World.GetVertexCoord(VertexCoordIds[0]);
                double[] co2 = World.GetVertexCoord(VertexCoordIds[1]);
                double[] co3 = World.GetVertexCoord(VertexCoordIds[2]);
                co1 = AddDisplacement(0, co1);
                co2 = AddDisplacement(1, co2);
                co3 = AddDisplacement(2, co3);
                OpenTK.Vector3d v1 = new OpenTK.Vector3d(co1[0], co1[1], co1[2]);
                OpenTK.Vector3d v2 = new OpenTK.Vector3d(co2[0], co2[1], co2[2]);
                OpenTK.Vector3d v3 = new OpenTK.Vector3d(co3[0], co3[1], co3[2]);

                var normal0 = CadUtils3D.TriNormal(v1, v2, v3);

                // X平面固定
                origin = new OpenTK.Vector3d(0.0, 0.0, 0.0);
                normal = new OpenTK.Vector3d(0.0, 0.0, 1.0);
                xdir = new OpenTK.Vector3d(1.0, 0.0, 0.0);
                System.Diagnostics.Debug.Assert(Math.Abs(Math.Abs(OpenTK.Vector3d.Dot(normal, normal0)) - 1.0) < 1.0e-12);
            }

            OpenTK.Vector2d pt2D = CadUtils3D.ProjectToPlane(p, origin, normal, xdir);
            return pt2D;
        }

        public OpenTK.Vector2d[] ProjectVertexsFrom3D()
        {
            System.Diagnostics.Debug.Assert(World.Dimension == 3);

            OpenTK.Vector3d origin;
            OpenTK.Vector3d normal;
            OpenTK.Vector3d xdir;
            if (FixedOrigin3D != null &&
                FixedNormal3D != null &&
                FixedXDir3D != null)
            {
                origin = new OpenTK.Vector3d(FixedOrigin3D[0], FixedOrigin3D[1], FixedOrigin3D[2]);
                normal = new OpenTK.Vector3d(FixedNormal3D[0], FixedNormal3D[1], FixedNormal3D[2]);
                xdir = new OpenTK.Vector3d(FixedXDir3D[0], FixedXDir3D[1], FixedXDir3D[2]);
            }
            else
            {
                double[] co1 = World.GetVertexCoord(VertexCoordIds[0]);
                double[] co2 = World.GetVertexCoord(VertexCoordIds[1]);
                double[] co3 = World.GetVertexCoord(VertexCoordIds[2]);
                co1 = AddDisplacement(0, co1);
                co2 = AddDisplacement(1, co2);
                co3 = AddDisplacement(2, co3);
                OpenTK.Vector3d v1 = new OpenTK.Vector3d(co1[0], co1[1], co1[2]);
                OpenTK.Vector3d v2 = new OpenTK.Vector3d(co2[0], co2[1], co2[2]);
                OpenTK.Vector3d v3 = new OpenTK.Vector3d(co3[0], co3[1], co3[2]);

                var normal0 = CadUtils3D.TriNormal(v1, v2, v3);

                // X平面固定
                origin = new OpenTK.Vector3d(0.0, 0.0, 0.0);
                normal = new OpenTK.Vector3d(0.0, 0.0, 1.0);
                xdir = new OpenTK.Vector3d(1.0, 0.0, 0.0);
                System.Diagnostics.Debug.Assert(Math.Abs(Math.Abs(OpenTK.Vector3d.Dot(normal, normal0)) - 1.0) < 1.0e-12);
            }

            OpenTK.Vector2d[] pt2Ds = new OpenTK.Vector2d[VertexCount];
            for (int iVertex = 0; iVertex < VertexCount; iVertex++)
            {
                int coId = VertexCoordIds[iVertex];
                double[] coord = World.GetVertexCoord(coId);
                coord = AddDisplacement(iVertex, coord);
                OpenTK.Vector3d p = new OpenTK.Vector3d(coord[0], coord[1], coord[2]);
                OpenTK.Vector2d pt2D = CadUtils3D.ProjectToPlane(p, origin, normal, xdir);
                pt2Ds[iVertex] = pt2D;

            }
            return pt2Ds;
        }

        public OpenTK.Vector2d[] ProjectNodesFrom3D()
        {
            System.Diagnostics.Debug.Assert(World.Dimension == 3);

            OpenTK.Vector3d origin;
            OpenTK.Vector3d normal;
            OpenTK.Vector3d xdir;
            if (FixedOrigin3D != null &&
                FixedNormal3D != null &&
                FixedXDir3D != null)
            {
                origin = new OpenTK.Vector3d(FixedOrigin3D[0], FixedOrigin3D[1], FixedOrigin3D[2]);
                normal = new OpenTK.Vector3d(FixedNormal3D[0], FixedNormal3D[1], FixedNormal3D[2]);
                xdir = new OpenTK.Vector3d(FixedXDir3D[0], FixedXDir3D[1], FixedXDir3D[2]);
            }
            else
            {
                double[] co1 = World.GetVertexCoord(VertexCoordIds[0]);
                double[] co2 = World.GetVertexCoord(VertexCoordIds[1]);
                double[] co3 = World.GetVertexCoord(VertexCoordIds[2]);
                co1 = AddDisplacement(0, co1);
                co2 = AddDisplacement(1, co2);
                co3 = AddDisplacement(2, co3);
                OpenTK.Vector3d v1 = new OpenTK.Vector3d(co1[0], co1[1], co1[2]);
                OpenTK.Vector3d v2 = new OpenTK.Vector3d(co2[0], co2[1], co2[2]);
                OpenTK.Vector3d v3 = new OpenTK.Vector3d(co3[0], co3[1], co3[2]);

                var normal0 = CadUtils3D.TriNormal(v1, v2, v3);

                // X平面固定
                origin = new OpenTK.Vector3d(0.0, 0.0, 0.0);
                normal = new OpenTK.Vector3d(0.0, 0.0, 1.0);
                xdir = new OpenTK.Vector3d(1.0, 0.0, 0.0);
                System.Diagnostics.Debug.Assert(Math.Abs(Math.Abs(OpenTK.Vector3d.Dot(normal, normal0)) - 1.0) < 1.0e-12);
            }

            OpenTK.Vector2d[] pt2Ds = new OpenTK.Vector2d[NodeCount];
            for (int iNode = 0; iNode < NodeCount; iNode++)
            {
                int coId = NodeCoordIds[iNode];
                double[] coord = World.GetCoord((uint)QuantityId, coId);
                coord = AddDisplacement(iNode, coord);
                OpenTK.Vector3d p = new OpenTK.Vector3d(coord[0], coord[1], coord[2]);
                OpenTK.Vector2d pt2D = CadUtils3D.ProjectToPlane(p, origin, normal, xdir);
                pt2Ds[iNode] = pt2D;
            }
            return pt2Ds;
        }
    }
}
