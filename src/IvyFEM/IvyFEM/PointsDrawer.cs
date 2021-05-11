using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    // DrawerArray, FieldDrawerArray兼用
    public class PointsDrawer : IDrawer, IFieldDrawer
    {
        public double[] Color { get; set; } = new double[] { 0.5, 0.5, 0.5 };

        public RotMode SutableRotMode => RotMode.RotMode2D;

        public bool IsAntiAliasing { get; set; } = false;

        // 位置と半径のリスト
        public IList<KeyValuePair<OpenTK.Vector2d, double>> PointRadiuss { get; set; } =
            new List<KeyValuePair<OpenTK.Vector2d, double>>();

        public PointsDrawer()
        {

        }

        public BoundingBox3D GetBoundingBox(OpenTK.Matrix3d rot)
        {
            if (PointRadiuss.Count == 0)
            {
                return new BoundingBox3D();
            }

            IList<double> coords = new List<double>();
            foreach (var pair in PointRadiuss)
            {
                OpenTK.Vector2d pt1 = pair.Key;
                coords.Add(pt1.X);
                coords.Add(pt1.Y);
            }

            int pointCnt = coords.Count / 2;

            BoundingBox3D bb;
            {
                double x1 = coords[0];
                double y1 = coords[1];
                double z1 = 0.0;
                bb = new BoundingBox3D(x1, x1, y1, y1, z1, z1);
            }
            for (int iPt = 1; iPt < pointCnt; iPt++)
            {
                double x1 = coords[iPt * 2];
                double y1 = coords[iPt * 2 + 1];
                double z1 = 0.0;
                bb.MaxX = (x1 > bb.MaxX) ? x1 : bb.MaxX; bb.MinX = (x1 < bb.MinX) ? x1 : bb.MinX;
                bb.MaxY = (y1 > bb.MaxY) ? y1 : bb.MaxY; bb.MinY = (y1 < bb.MinY) ? y1 : bb.MinY;
                bb.MaxZ = (z1 > bb.MaxZ) ? z1 : bb.MaxZ; bb.MinZ = (z1 < bb.MinZ) ? z1 : bb.MinZ;
            }
            return bb;
        }

        public void Draw()
        {
            foreach (var pair in PointRadiuss)
            {
                var pt = pair.Key;
                var r = pair.Value;
                int cnt = 36;
                OpenTK.Vector2d[] pts = new OpenTK.Vector2d[cnt];
                for (int i = 0; i < cnt; i++)
                {
                    double theta = i * (2.0 * Math.PI / cnt);
                    double x = pt.X + r * Math.Cos(theta);
                    double y = pt.Y + r * Math.Sin(theta);
                    pts[i] = new OpenTK.Vector2d(x, y);
                }
                GL.LineWidth(4);
                GL.Color3(Color);
                GL.Begin(PrimitiveType.Lines);
                for (int i = 0; i < cnt; i++)
                {
                    var pt1 = pts[i];
                    var pt2 = pts[(i + 1) % cnt];
                    GL.Vertex2(pt1.X, pt1.Y);
                    GL.Vertex2(pt2.X, pt2.Y);
                }
                GL.End();
            }
        }

        public void DrawSelection(uint idraw)
        {

        }

        public void AddSelected(int[] selectFlag)
        {

        }

        public void ClearSelected()
        {

        }

        public void Update(FEWorld world)
        {

        }
    }
}
