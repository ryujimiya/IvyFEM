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
    }
}
