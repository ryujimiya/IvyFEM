using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Vertex3D : IObject
    {
        public OpenTK.Vector3d Point { get; set; }
        public double[] Color { get; set; } = new double[3];

        public Vertex3D()
        {
            Point = new OpenTK.Vector3d();
            Color = new double[3] { 0.0, 0.0, 0.0 };
        }

        public Vertex3D(OpenTK.Vector3d point)
        {
            Point = new OpenTK.Vector3d(point.X, point.Y, point.Z);
            Color = new double[3] { 0.0, 0.0, 0.0 };
        }

        public Vertex3D(Vertex3D src)
        {
            Copy(src);
        }

        public void Copy(IObject src)
        {
            Vertex3D srcVertex = src as Vertex3D;

            Point = new OpenTK.Vector3d(srcVertex.Point.X, srcVertex.Point.Y, srcVertex.Point.Z);
            Color = new double[srcVertex.Color.Length];
            srcVertex.Color.CopyTo(Color, 0);
        }
    }
}
