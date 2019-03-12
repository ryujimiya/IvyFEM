using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    abstract public class Constraint
    {
        public EqualityType Equality { get; set; } = EqualityType.Eq;
        public abstract double GetValue(double[] x);
        public abstract double GetDerivative(int iDof, double[] x);
        public abstract double Get2ndDerivative(int iDof, int jDof, double[] x);
    }
}
