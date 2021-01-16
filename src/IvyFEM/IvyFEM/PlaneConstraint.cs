using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class PlaneConstraint : Constraint
    {
        public OpenTK.Vector3d Point { get; set; } = new OpenTK.Vector3d();
        public OpenTK.Vector3d Normal { get; set; } = new OpenTK.Vector3d();
        public double BoundingBoxLen { get; set; } = 1.0;
        public override RotMode SutableRotMode => RotMode.RotMode3D;

        public PlaneConstraint()
        {

        }

        public PlaneConstraint(OpenTK.Vector3d point, OpenTK.Vector3d normal,
            EqualityType equality = EqualityType.Eq)
        {
            Equality = equality;
            Point = new OpenTK.Vector3d(point);
            Normal = new OpenTK.Vector3d(normal);
        }

        public override double GetValue(double[] x)
        {
            System.Diagnostics.Debug.Assert(x.Length == 3);
            OpenTK.Vector3d xVec = new OpenTK.Vector3d(x[0], x[1], x[2]);
            OpenTK.Vector3d lineVec = xVec - Point;
            double value = -OpenTK.Vector3d.Dot(lineVec, Normal);
            if (Equality == EqualityType.LessEq)
            {
                value *= -1.0;
            }
            return value;
        }

        public override double GetDerivative(int iDof, double[] x)
        {
            System.Diagnostics.Debug.Assert(x.Length == 3);
            if (iDof >= 3)
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
            System.Diagnostics.Debug.Assert(x.Length == 3);
            if (iDof >= 3 || jDof >= 3)
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
            var pt = Point;
            var normal = Normal;
            var horizontalX = CadUtils3D.GetVerticalUnitVector(normal);
            var horizontalY = OpenTK.Vector3d.Cross(normal, horizontalX);
            OpenTK.Vector3d[] pts = new OpenTK.Vector3d[4];
            pts[0] = pt - horizontalX * BoundingBoxLen - horizontalY * BoundingBoxLen;
            pts[1] = pt + horizontalX * BoundingBoxLen - horizontalY * BoundingBoxLen;
            pts[2] = pt + horizontalX * BoundingBoxLen + horizontalY * BoundingBoxLen;
            pts[3] = pt - horizontalX * BoundingBoxLen + horizontalY * BoundingBoxLen;

            IList<double> coords = new List<double>();
            foreach(OpenTK.Vector3d tmpPt in pts)
            {
                coords.Add(tmpPt.X);
                coords.Add(tmpPt.Y);
                coords.Add(tmpPt.Z);
            }
            int pointCnt = coords.Count / 3;

            BoundingBox3D bb;
            {
                double x1 = coords[0];
                double y1 = coords[1];
                double z1 = coords[2];
                bb = new BoundingBox3D(x1, x1, y1, y1, z1, z1);
            }
            for (int iPt = 1; iPt < pointCnt; iPt++)
            {
                double x1 = coords[iPt * 3];
                double y1 = coords[iPt * 3 + 1];
                double z1 = coords[iPt * 3 + 2];
                bb.MaxX = (x1 > bb.MaxX) ? x1 : bb.MaxX; bb.MinX = (x1 < bb.MinX) ? x1 : bb.MinX;
                bb.MaxY = (y1 > bb.MaxY) ? y1 : bb.MaxY; bb.MinY = (y1 < bb.MinY) ? y1 : bb.MinY;
                bb.MaxZ = (z1 > bb.MaxZ) ? z1 : bb.MaxZ; bb.MinZ = (z1 < bb.MinZ) ? z1 : bb.MinZ;
            }
            return bb;
        }
    }
}
