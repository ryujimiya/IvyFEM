using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class LineConstraint : Constraint
    {
        public OpenTK.Vector2d Point { get; set; } = new OpenTK.Vector2d();
        public OpenTK.Vector2d Normal { get; set; } = new OpenTK.Vector2d();
        public double BoundingBoxLen { get; set; } = 1.0;

        public LineConstraint()
        {

        }

        public LineConstraint(OpenTK.Vector2d point, OpenTK.Vector2d normal,
            EqualityType equality = EqualityType.Eq)
        {
            Equality = equality;
            Point = new OpenTK.Vector2d(point.X, point.Y);
            Normal = new OpenTK.Vector2d(normal.X, normal.Y);
        }

        public override double GetValue(double[] x)
        {
            System.Diagnostics.Debug.Assert(x.Length == 2);
            OpenTK.Vector2d xVec = new OpenTK.Vector2d(x[0], x[1]);
            OpenTK.Vector2d lineVec = xVec - Point;
            double value = -OpenTK.Vector2d.Dot(lineVec, Normal);
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
            double value = -Normal[iDof];
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
            double value = 0;
            if (Equality == EqualityType.LessEq)
            {
                value *= -1.0;
            }
            return value;
        }

        public override BoundingBox3D GetBoundingBox(OpenTK.Matrix3d rot)
        {
            // 2D
            var pt = Point;
            var normal = Normal;
            var horizontal = new OpenTK.Vector2d(-normal.Y, normal.X);
            var pt1 = pt - horizontal * BoundingBoxLen;
            var pt2 = pt + horizontal * BoundingBoxLen;

            IList<double> coords = new List<double>();
            coords.Add(pt1.X);
            coords.Add(pt1.Y);
            coords.Add(pt2.X);
            coords.Add(pt2.Y);
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
