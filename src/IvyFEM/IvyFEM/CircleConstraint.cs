using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class CircleConstraint : Constraint
    {
        public OpenTK.Vector2d Point { get; set; } = new OpenTK.Vector2d();
        public double R { get; set; } = 0;

        public CircleConstraint()
        {

        }

        public CircleConstraint(OpenTK.Vector2d point, double r,
            EqualityType equality = EqualityType.Eq)
        {
            Equality = equality;
            Point = new OpenTK.Vector2d(point.X, point.Y);
            R = r;
        }

        public override double GetValue(double[] x)
        {
            System.Diagnostics.Debug.Assert(x.Length == 2);
            OpenTK.Vector2d xVec = new OpenTK.Vector2d(x[0], x[1]);
            OpenTK.Vector2d rVec = xVec - Point;
            double r = rVec.Length;
            double value = R - r;
            if (Equality == EqualityType.LessEq)
            {
                value *= -1.0;
            }
            return value;
        }

        public override double GetDerivative(int iDof, double[] x)
        {
            System.Diagnostics.Debug.Assert(x.Length == 2);
            if (iDof >= 2)
            {
                throw new InvalidOperationException();
            }
            OpenTK.Vector2d xVec = new OpenTK.Vector2d(x[0], x[1]);
            OpenTK.Vector2d rVec = xVec - Point;
            double r = rVec.Length;
            double[] rVec2 = { rVec.X, rVec.Y };
            double value = -rVec2[iDof] / r;
            if (Equality == EqualityType.LessEq)
            {
                value *= -1.0;
            }
            return value;
        }

        public override double Get2ndDerivative(int iDof, int jDof, double[] x)
        {
            System.Diagnostics.Debug.Assert(x.Length == 2);
            if (iDof >= 2 || jDof >= 2)
            {
                throw new InvalidOperationException();
            }
            OpenTK.Vector2d xVec = new OpenTK.Vector2d(x[0], x[1]);
            OpenTK.Vector2d rVec = xVec - Point;
            double r = rVec.Length;
            double[] rVec2 = { rVec.X, rVec.Y };
            double value = 0;
            if (iDof  == jDof)
            {
                value += -1.0 / r;
            }
            value += rVec2[iDof] * rVec2[jDof] / (r * r * r);
            if (Equality == EqualityType.LessEq)
            {
                value *= -1.0;
            }
            return value;
        }

        public override BoundingBox3D GetBoundingBox(OpenTK.Matrix3d rot)
        {
            // 2D
            IList<double> coords = new List<double>();
            coords.Add(Point.X + R);
            coords.Add(0.0);
            coords.Add(0.0);
            coords.Add(Point.Y + R);
            coords.Add(Point.X - R);
            coords.Add(0.0);
            coords.Add(0.0);
            coords.Add(Point.Y - R);
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
    }
}
