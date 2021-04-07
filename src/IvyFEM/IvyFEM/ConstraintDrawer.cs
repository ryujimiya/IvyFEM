using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    // constraintはDrawerArray, FieldDrawerArray兼用
    public class ConstraintDrawer : IDrawer, IFieldDrawer
    {
        private Constraint Constraint = null;

        public RotMode SutableRotMode { get; private set; } = RotMode.RotModeNotSet;
        public bool IsAntiAliasing { get; set; } = false;

        public ConstraintDrawer(Constraint constraint)
        {
            Constraint = constraint;
            SutableRotMode = constraint.SutableRotMode;
        }

        public BoundingBox3D GetBoundingBox(OpenTK.Matrix3d rot)
        {
            return Constraint.GetBoundingBox(rot);
        }

        public void Draw()
        {
            if (Constraint is LineConstraint)
            {
                DrawLineConstraint();
            }
            else if (Constraint is LineSegConstraint)
            {
                DrawLineSegConstraint();
            }
            else if (Constraint is CircleConstraint)
            {
                DrawCircleConstraint();
            }
            else if (Constraint is PlaneConstraint)
            {
                DrawPlaneConstraint();
            }
            else if (Constraint is SphereConstraint)
            {
                DrawSphereConstraint();
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

        private void DrawLineSegConstraint()
        {
            LineSegConstraint lineSegConstraint = Constraint as LineSegConstraint;
            var pt1 = lineSegConstraint.Point1;
            var pt2 = lineSegConstraint.Point2;
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

        private void DrawPlaneConstraint()
        {
            PlaneConstraint planeConstraint = Constraint as PlaneConstraint;
            var pt = planeConstraint.Point;
            var normal = planeConstraint.Normal;
            double bbLen = planeConstraint.BoundingBoxLen;
            var horizontalX = CadUtils3D.GetVerticalUnitVector(normal);
            var horizontalY = OpenTK.Vector3d.Cross(normal, horizontalX);
            var pt1 = pt - horizontalX * bbLen - horizontalY * bbLen;
            var pt2 = pt + horizontalX * bbLen - horizontalY * bbLen;
            var pt3 = pt + horizontalX * bbLen + horizontalY * bbLen;
            var pt4 = pt - horizontalX * bbLen + horizontalY * bbLen;
            GL.LineWidth(2);
            GL.Color3(0.5, 0.5, 0.5);
            GL.Begin(PrimitiveType.Polygon);
            GL.Vertex3(pt1.X, pt1.Y, pt1.Z);
            GL.Vertex3(pt2.X, pt2.Y, pt2.Z);
            GL.Vertex3(pt3.X, pt3.Y, pt3.Z);
            GL.Vertex3(pt4.X, pt4.Y, pt4.Z);
            GL.End();
        }

        private void DrawSphereConstraint()
        {
            SphereConstraint sphereConstraint = Constraint as SphereConstraint;
            var cPt = sphereConstraint.Point;
            var r = sphereConstraint.R;

            OpenTK.Vector3d ptU = new OpenTK.Vector3d(cPt.X, cPt.Y, cPt.Z + r);
            OpenTK.Vector3d ptD = new OpenTK.Vector3d(cPt.X, cPt.Y, cPt.Z - r);
            int divCnt = 36;
            IList<IList<OpenTK.Vector3d>> ptss = new List<IList<OpenTK.Vector3d>>();
            int divCntXY = divCnt;
            int divCntZY = divCnt / 2;
            for (int iXY = 0; iXY < divCntXY; iXY++)
            {
                double thetaXY = (iXY / (double)divCntXY) * 360.0 * (Math.PI / 180.0);
                IList<OpenTK.Vector3d> pts = new List<OpenTK.Vector3d>();
                ptss.Add(pts);

                pts.Add(ptU);
                for (int iZY = 1; iZY <= (divCntZY - 1); iZY++)
                {
                    double thetaZY = (iZY / (double)divCntZY) * 180.0 * (Math.PI / 180);
                    double x = cPt.X + r * Math.Sin(thetaZY) * Math.Cos(thetaXY);
                    double y = cPt.Y + r * Math.Sin(thetaZY) * Math.Sin(thetaXY);
                    double z = cPt.Z + r * Math.Cos(thetaZY);
                    OpenTK.Vector3d pt = new OpenTK.Vector3d(x, y, z);
                    pts.Add(pt);
                }
                pts.Add(ptD);
            }

            GL.LineWidth(1);
            GL.Color3(0.5, 0.5, 0.5);
            GL.Begin(PrimitiveType.Lines);
            for (int iXY = 0; iXY < ptss.Count; iXY++)
            {
                IList<OpenTK.Vector3d> pts = ptss[iXY];
                for (int iZ = 0; iZ < pts.Count - 1; iZ++)
                {
                    var pt1 = pts[iZ];
                    var pt2 = pts[iZ + 1];
                    GL.Vertex3(pt1.X, pt1.Y, pt1.Z);
                    GL.Vertex3(pt2.X, pt2.Y, pt2.Z);
                }
            }
            for (int iXY = 0; iXY < ptss.Count - 1; iXY++)
            {
                IList<OpenTK.Vector3d> pts = ptss[iXY];
                IList<OpenTK.Vector3d> nextPts = ptss[iXY + 1];
                for (int iZ = 1; iZ < pts.Count - 1; iZ++)
                {
                    var pt1 = pts[iZ];
                    var pt2 = nextPts[iZ];
                    GL.Vertex3(pt1.X, pt1.Y, pt1.Z);
                    GL.Vertex3(pt2.X, pt2.Y, pt2.Z);
                }
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
