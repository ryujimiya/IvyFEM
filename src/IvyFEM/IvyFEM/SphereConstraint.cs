using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class SphereConstraint : Constraint
    {
        public OpenTK.Vector3d Point { get; set; } = new OpenTK.Vector3d();
        public double R { get; set; } = 0;
        public override RotMode SutableRotMode => RotMode.RotMode2D;

        public SphereConstraint()
        {

        }

        public SphereConstraint(OpenTK.Vector3d point, double r,
            EqualityType equality = EqualityType.Eq)
        {
            Equality = equality;
            Point = new OpenTK.Vector3d(point);
            R = r;
        }

        public override double GetValue(double[] x)
        {
            System.Diagnostics.Debug.Assert(x.Length == 3);
            OpenTK.Vector3d xVec = new OpenTK.Vector3d(x[0], x[1], x[2]);
            OpenTK.Vector3d rVec = xVec - Point;
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
            System.Diagnostics.Debug.Assert(x.Length == 3);
            if (iDof >= 3)
            {
                throw new InvalidOperationException();
            }
            OpenTK.Vector3d xVec = new OpenTK.Vector3d(x[0], x[1], x[2]);
            OpenTK.Vector3d rVec = xVec - Point;
            double r = rVec.Length;
            double[] rVecXYZ = { rVec.X, rVec.Y, rVec.Z };
            double value = -rVecXYZ[iDof] / r;
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
            OpenTK.Vector3d xVec = new OpenTK.Vector3d(x[0], x[1], x[2]);
            OpenTK.Vector3d rVec = xVec - Point;
            double r = rVec.Length;
            double[] rVecXYZ = { rVec.X, rVec.Y, rVec.Z };
            double value = 0;
            if (iDof  == jDof)
            {
                value += -1.0 / r;
            }
            value += rVecXYZ[iDof] * rVecXYZ[jDof] / (r * r * r);
            if (Equality == EqualityType.LessEq)
            {
                value *= -1.0;
            }
            return value;
        }

        public override BoundingBox3D GetBoundingBox(OpenTK.Matrix3d rot)
        {
            IList<double> coords = new List<double>();
            coords.Add(Point.X + R);
            coords.Add(Point.Y);
            coords.Add(Point.Z);

            coords.Add(Point.X);
            coords.Add(Point.Y + R);
            coords.Add(Point.Z);

            coords.Add(Point.X - R);
            coords.Add(Point.Y);
            coords.Add(Point.Z);

            coords.Add(Point.X);
            coords.Add(Point.Y - R);
            coords.Add(Point.Z);

            coords.Add(Point.X);
            coords.Add(Point.Y);
            coords.Add(Point.Z + R);

            coords.Add(Point.X);
            coords.Add(Point.Y - R);
            coords.Add(Point.Z);

            coords.Add(Point.X);
            coords.Add(Point.Y);
            coords.Add(Point.Z - R);

            coords.Add(Point.X);
            coords.Add(Point.Y + R);
            coords.Add(Point.Z);

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
