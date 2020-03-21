using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    // constraintはDrawerArray, FieldDrawerArray兼用
    public class ConstraintDrawer : IDrawer, IFieldDrawer
    {
        private Constraint Constraint = null;

        public RotMode SutableRotMode { get; private set; } = RotMode.RotMode2D;
        public bool IsAntiAliasing { get; set; } = false;

        public ConstraintDrawer(Constraint constraint)
        {
            Constraint = constraint;
        }

        public BoundingBox3D GetBoundingBox(Matrix3d rot)
        {
            return Constraint.GetBoundingBox(rot);
        }

        public void Draw()
        {
            if (Constraint is LineConstraint)
            {
                DrawLineConstraint();
            }
            else if (Constraint is CircleConstraint)
            {
                DrawCircleConstraint();
            }
        }

        private void DrawLineConstraint()
        {
            LineConstraint lineConstraint = Constraint as LineConstraint;
            var pt = lineConstraint.Point;
            var normal = lineConstraint.Normal;
            var horizontal = new OpenTK.Vector2d(-normal.Y, normal.X);
            var pt1 = pt - horizontal * 1.0e+6;
            var pt2 = pt + horizontal * 1.0e+6;
            GL.LineWidth(2);
            GL.Color3(0.5, 0.5, 0.5);
            GL.Begin(PrimitiveType.Lines);
            GL.Vertex2(pt1.X, pt1.Y);
            GL.Vertex2(pt2.X, pt2.Y);
            GL.End();
        }

        private void DrawCircleConstraint()
        {
            CircleConstraint circleConstraint = Constraint as CircleConstraint;
            var pt = circleConstraint.Point;
            var r = circleConstraint.R;
            int cnt = 36;
            OpenTK.Vector2d[] pts = new OpenTK.Vector2d[cnt];
            for (int i = 0; i < cnt; i++)
            {
                double theta = i * (2.0 * Math.PI / cnt);
                double x = pt.X + r * Math.Cos(theta);
                double y = pt.Y + r * Math.Sin(theta);
                pts[i] = new OpenTK.Vector2d(x, y);
            }
            GL.LineWidth(2);
            GL.Color3(0.5, 0.5, 0.5);
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
