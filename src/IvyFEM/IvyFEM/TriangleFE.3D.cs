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

            double[] co1 = World.GetVertexCoord(VertexCoordIds[0]);
            double[] co2 = World.GetVertexCoord(VertexCoordIds[1]);
            double[] co3 = World.GetVertexCoord(VertexCoordIds[2]);
            OpenTK.Vector3d v1 = new OpenTK.Vector3d(co1[0], co1[1], co1[2]);
            OpenTK.Vector3d v2 = new OpenTK.Vector3d(co2[0], co2[1], co2[2]);
            OpenTK.Vector3d v3 = new OpenTK.Vector3d(co3[0], co3[1], co3[2]);

            OpenTK.Vector3d origin = new OpenTK.Vector3d(co1[0], co1[1], co1[2]);
            OpenTK.Vector3d normal = CadUtils3D.TriNormal(v1, v2, v3);
            OpenTK.Vector3d xdir = CadUtils3D.GetDirection(v1, v2);

            double x = OpenTK.Vector3d.Dot((p - origin), xdir);
            double y = OpenTK.Vector3d.Dot((p - origin), OpenTK.Vector3d.Cross(normal, xdir));
            return new OpenTK.Vector2d(x, y);
        }

        public OpenTK.Vector2d[] ProjectVertexsFrom3D()
        {
            System.Diagnostics.Debug.Assert(World.Dimension == 3);

            double[] co1 = World.GetVertexCoord(VertexCoordIds[0]);
            double[] co2 = World.GetVertexCoord(VertexCoordIds[1]);
            double[] co3 = World.GetVertexCoord(VertexCoordIds[2]);
            OpenTK.Vector3d v1 = new OpenTK.Vector3d(co1[0], co1[1], co1[2]);
            OpenTK.Vector3d v2 = new OpenTK.Vector3d(co2[0], co2[1], co2[2]);
            OpenTK.Vector3d v3 = new OpenTK.Vector3d(co3[0], co3[1], co3[2]);

            OpenTK.Vector3d origin = new OpenTK.Vector3d(co1[0], co1[1], co1[2]);
            OpenTK.Vector3d normal = CadUtils3D.TriNormal(v1, v2, v3);
            OpenTK.Vector3d xdir = CadUtils3D.GetDirection(v1, v2);

            OpenTK.Vector2d[] pt2Ds = new OpenTK.Vector2d[3];
            for (int iVertex = 0; iVertex < VertexCount; iVertex++)
            {
                int coId = VertexCoordIds[iVertex];
                double[] coord = World.GetVertexCoord(coId);
                OpenTK.Vector3d p = new OpenTK.Vector3d(coord[0], coord[1], coord[2]);
                double x = OpenTK.Vector3d.Dot((p - origin), xdir);
                double y = OpenTK.Vector3d.Dot((p - origin), OpenTK.Vector3d.Cross(normal, xdir));
                pt2Ds[iVertex] = new OpenTK.Vector2d(x, y);

            }
            return pt2Ds;
        }

        public OpenTK.Vector2d[] ProjectNodesFrom3D()
        {
            System.Diagnostics.Debug.Assert(World.Dimension == 3);

            double[] co1 = World.GetVertexCoord(VertexCoordIds[0]);
            double[] co2 = World.GetVertexCoord(VertexCoordIds[1]);
            double[] co3 = World.GetVertexCoord(VertexCoordIds[2]);
            OpenTK.Vector3d v1 = new OpenTK.Vector3d(co1[0], co1[1], co1[2]);
            OpenTK.Vector3d v2 = new OpenTK.Vector3d(co2[0], co2[1], co2[2]);
            OpenTK.Vector3d v3 = new OpenTK.Vector3d(co3[0], co3[1], co3[2]);

            OpenTK.Vector3d origin = new OpenTK.Vector3d(co1[0], co1[1], co1[2]);
            OpenTK.Vector3d normal = CadUtils3D.TriNormal(v1, v2, v3);
            OpenTK.Vector3d xdir = CadUtils3D.GetDirection(v1, v2);

            OpenTK.Vector2d[] pt2Ds = new OpenTK.Vector2d[3];
            for (int iNode = 0; iNode < NodeCount; iNode++)
            {
                int coId = NodeCoordIds[iNode];
                double[] coord = World.GetCoord((uint)QuantityId, coId);
                OpenTK.Vector3d p = new OpenTK.Vector3d(coord[0], coord[1], coord[2]);
                double x = OpenTK.Vector3d.Dot((p - origin), xdir);
                double y = OpenTK.Vector3d.Dot((p - origin), OpenTK.Vector3d.Cross(normal, xdir));
                pt2Ds[iNode] = new OpenTK.Vector2d(x, y);
            }
            return pt2Ds;
        }
    }
}
