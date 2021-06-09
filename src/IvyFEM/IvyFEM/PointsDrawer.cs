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
        public uint Dimension { get; private set; } = 0;
        public double[] Color { get; set; } = new double[] { 0.5, 0.5, 0.5 };

        public RotMode SutableRotMode { get; private set; } = RotMode.RotModeNotSet;

        public bool IsAntiAliasing { get; set; } = false;

        // 位置と半径のリスト
        public IList<KeyValuePair<double[], double>> PointRadiuss { get; set; } =
            new List<KeyValuePair<double[], double>>();

        public PointsDrawer()
        {

        }

        public BoundingBox3D GetBoundingBox(OpenTK.Matrix3d rot)
        {
            if (PointRadiuss.Count == 0)
            {
                return new BoundingBox3D();
            }
            Dimension = (uint)PointRadiuss[0].Key.Length;
            SutableRotMode = Dimension == 3 ? RotMode.RotMode3D : RotMode.RotMode2D;

            IList<double> coords = new List<double>();
            foreach (var pair in PointRadiuss)
            {
                double[] pt1 = pair.Key;
                for (int idim = 0; idim < Dimension; idim++)
                {
                    coords.Add(pt1[idim]);
                }
            }

            int pointCnt = coords.Count / (int)Dimension;

            BoundingBox3D bb;
            {
                double x1;
                double y1;
                double z1;
                if (Dimension == 3)
                {
                    x1 = coords[0];
                    y1 = coords[1];
                    z1 = coords[2];
                }
                else
                {
                    x1 = coords[0];
                    y1 = coords[1];
                    z1 = 0.0;
                }
                bb = new BoundingBox3D(x1, x1, y1, y1, z1, z1);
            }
            for (int iPt = 1; iPt < pointCnt; iPt++)
            {
                double x1;
                double y1;
                double z1;
                if (Dimension == 3)
                {
                    x1 = coords[iPt * 3];
                    y1 = coords[iPt * 3 + 1];
                    z1 = coords[iPt * 3 + 2];
                }
                else
                {
                    x1 = coords[iPt * 2];
                    y1 = coords[iPt * 2 + 1];
                    z1 = 0.0;
                }
                bb.MaxX = (x1 > bb.MaxX) ? x1 : bb.MaxX; bb.MinX = (x1 < bb.MinX) ? x1 : bb.MinX;
                bb.MaxY = (y1 > bb.MaxY) ? y1 : bb.MaxY; bb.MinY = (y1 < bb.MinY) ? y1 : bb.MinY;
                bb.MaxZ = (z1 > bb.MaxZ) ? z1 : bb.MaxZ; bb.MinZ = (z1 < bb.MinZ) ? z1 : bb.MinZ;
            }
            return bb;
        }

        public void Draw()
        {
            if (PointRadiuss.Count == 0)
            {
                return;
            }
            Dimension = (uint)PointRadiuss[0].Key.Length;
            SutableRotMode = Dimension == 3 ? RotMode.RotMode3D : RotMode.RotMode2D;

            if (Dimension == 3)
            {
                Draw3D();
            }
            else
            {
                Draw2D();
            }
        }

        private void Draw2D()
        {
            foreach (var pair in PointRadiuss)
            {
                double[] pt = pair.Key;
                double r = pair.Value;
                int cnt = 36;
                OpenTK.Vector2d[] pts = new OpenTK.Vector2d[cnt];
                for (int i = 0; i < cnt; i++)
                {
                    double theta = i * (2.0 * Math.PI / cnt);
                    double x = pt[0] + r * Math.Cos(theta);
                    double y = pt[1] + r * Math.Sin(theta);
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

        private void Draw3D()
        {
            foreach (var pair in PointRadiuss)
            {
                double[] pt = pair.Key;
                double r = pair.Value;
                int cnt = 36;
                OpenTK.Vector3d[] pts = new OpenTK.Vector3d[cnt];
                for (int i = 0; i < cnt; i++)
                {
                    double theta = i * (2.0 * Math.PI / cnt);
                    double x = pt[0] + r * Math.Cos(theta);
                    double y = pt[1] + r * Math.Sin(theta);
                    double z = pt[2];
                    pts[i] = new OpenTK.Vector3d(x, y, z);
                }
                GL.LineWidth(4);
                GL.Color3(Color);
                GL.Begin(PrimitiveType.Lines);
                for (int i = 0; i < cnt; i++)
                {
                    var pt1 = pts[i];
                    var pt2 = pts[(i + 1) % cnt];
                    GL.Vertex3(pt1.X, pt1.Y, pt1.Z);
                    GL.Vertex3(pt2.X, pt2.Y, pt2.Z);
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
