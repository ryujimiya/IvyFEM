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
    }
}
