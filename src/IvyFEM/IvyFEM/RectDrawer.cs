using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    public class RectDrawer : IDrawer
    {
        private double BeginX;
        private double BeginY;
        private double EndX;
        private double EndY;

        public RotMode SutableRotMode => RotMode.RotMode2D;

        public bool IsAntiAliasing { get; set; }

        public RectDrawer()
        {
            BeginX = 0;
            BeginY = 0;
        }

        public RectDrawer(double x, double y)
        {
            EndX = x;
            EndY = y;
            BeginX = x;
            BeginY = y;
        }

        public void SetInitialPositionMode(double x0, double y0)
        {
            BeginX = x0;
            BeginY = y0;
        }

        public void SetPosition(double x1, double y1)
        {
            EndX = x1;
            EndY = y1;
        }

        public void GetPosition(out double x0, out double y0, out double x1, out double y1)
        {
            x0 = BeginX;
            y0 = BeginY;
            x1 = EndX;
            y1 = EndY;
        }

        public void GetCenterSize(out double centX, out double centY, out double sizeX, out double sizeY)
        {
            centX = (BeginX + EndX) * 0.5;
            centY = (BeginY + EndY) * 0.5;
            sizeX = Math.Abs(BeginX - EndX);
            sizeY = Math.Abs(BeginY - EndY);
        }

        public void AddSelected(int[] selectFlag)
        {

        }

        public void ClearSelected()
        {

        }

        public void Draw()
        {
            GL.Color3(0.0, 0.0, 0.0);
            GL.LineWidth(1);
            double offset = 0.1;

            GL.Begin(PrimitiveType.Lines);

            GL.Vertex3(BeginX, BeginY, offset);
            GL.Vertex3(BeginX, EndY, offset);
            GL.Vertex3(BeginX, EndY, offset);
            GL.Vertex3(EndX, EndY, offset);
            GL.Vertex3(EndX, EndY, offset);
            GL.Vertex3(EndX, BeginY, offset);
            GL.Vertex3(EndX, BeginY, offset);
            GL.Vertex3(BeginX, BeginY, offset);

            GL.End();
        }

        public void DrawSelection(uint idraw)
        {

        }

        public BoundingBox3D GetBoundingBox(OpenTK.Matrix3d rot)
        {
            return new BoundingBox3D();
        }
    }
}
