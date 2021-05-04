﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class LineSegConstraint : Constraint
    {
        public OpenTK.Vector2d Point1 { get; set; } = new OpenTK.Vector2d();
        public OpenTK.Vector2d Point2 { get; set; } = new OpenTK.Vector2d();
        public OpenTK.Vector2d Normal { get; set; } = new OpenTK.Vector2d();
        public override RotMode SutableRotMode => RotMode.RotMode2D;

        public LineSegConstraint()
        {

        }

        public LineSegConstraint(OpenTK.Vector2d point1, OpenTK.Vector2d point2, OpenTK.Vector2d normal,
            EqualityType equality = EqualityType.Eq)
        {
            Equality = equality;
            Point1 = new OpenTK.Vector2d(point1.X, point1.Y);
            Point2 = new OpenTK.Vector2d(point2.X, point2.Y);
            Normal = new OpenTK.Vector2d(normal.X, normal.Y);
        }

        private bool IsInsideLineSeg(double[] x)
        {
            // 線分上にある？
            /* FIXME:
            System.Diagnostics.Debug.Assert(x.Length == 2);
            OpenTK.Vector2d xVec = new OpenTK.Vector2d(x[0], x[1]);
            OpenTK.Vector2d segDir = OpenTK.Vector2d.Normalize(Point2 - Point1);
            double scale = OpenTK.Vector2d.Dot(xVec, segDir);
            if (scale < 0.0 || scale > segDir.Length)
            {
                return false;
            }
            return true;
            */
            return true;
        }

        public override double GetValue(double[] x)
        {
            if (!IsInsideLineSeg(x)) // FIXME:
            {
                return 0.0;
            }
            System.Diagnostics.Debug.Assert(x.Length == 2);
            OpenTK.Vector2d xVec = new OpenTK.Vector2d(x[0], x[1]);
            OpenTK.Vector2d lineVec = xVec - Point1;
            double value = -OpenTK.Vector2d.Dot(lineVec, Normal);
            if (Equality == EqualityType.LessEq)
            {
                value *= -1.0;
            }
            return value;
        }

        public override double GetDerivative(int iDof, double[] x)
        {
            if (!IsInsideLineSeg(x)) // FIXME:
            {
                return 0.0;
            }
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
            if (!IsInsideLineSeg(x)) // FIXME:
            {
                return 0.0;
            }
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
            var pt1 = Point1;
            var pt2 = Point2;

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