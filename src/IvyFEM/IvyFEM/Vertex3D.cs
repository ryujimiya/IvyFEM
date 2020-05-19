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

        public Vertex3D()
        {
            Point = new OpenTK.Vector3d();
        }

        public Vertex3D(OpenTK.Vector3d point)
        {
            Point = new OpenTK.Vector3d(point.X, point.Y, point.Z);
        }

        public Vertex3D(Vertex3D src)
        {
            Copy(src);
        }

        public void Copy(IObject src)
        {
            Vertex3D srcVertex = src as Vertex3D;

            Point = new OpenTK.Vector3d(srcVertex.Point.X, srcVertex.Point.Y, srcVertex.Point.Z);
        }
    }
}
