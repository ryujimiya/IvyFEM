using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshPoint2D
    {
        public OpenTK.Vector2d Point { get; set; } = new OpenTK.Vector2d();
        public int Element { get; set; } = -1;
        public uint Node { get; set; } = 0;

        public MeshPoint2D()
        {

        }

        public MeshPoint2D(double x, double y, int element, uint node)
        {
            Point = new OpenTK.Vector2d(x, y);
            Element = element;
            Node = node;
        }

        public MeshPoint2D(MeshPoint2D src)
        {
            Point = src.Point;
            Element = src.Element;
            Node = src.Node;
        }
    }

}
