using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public interface IEquationSolver
    {
        bool DoubleSolve(out double[] X, IvyFEM.Linear.DoubleSparseMatrix A, double[] B);
        bool ComplexSolve(
            out System.Numerics.Complex[] X, IvyFEM.Linear.ComplexSparseMatrix A, System.Numerics.Complex[] B);
    }
}
