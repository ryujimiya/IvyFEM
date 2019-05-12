using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public partial class IvyFEMEquationSolver : IEquationSolver
    {
        public IvyFEMEquationSolverMethod Method { get; set; } = IvyFEMEquationSolverMethod.Default;
        public double ConvRatioTolerance { get; set; } = IvyFEM.Linear.Constants.ConvRatioTolerance;
        public int ILUFillinLevel { get; set; } = 0;

        public bool DoubleSolve(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            bool success = false;
            X = null;

            switch (Method)
            {
                case IvyFEMEquationSolverMethod.Default:
                case IvyFEMEquationSolverMethod.NoPreconCG:
                    success = DoubleSolveNoPreconCG(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.CG:
                    success = DoubleSolveCG(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.CGWithPivoting:
                    success = DoubleSolveCGWithPivoting(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.ICCG:
                    success = DoubleSolveICCG(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.NoPreconBiCGSTAB:
                    success = DoubleSolveNoPreconBiCGSTAB(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.BiCGSTAB:
                    success = DoubleSolveBiCGSTAB(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.BiCGSTABWithPivoting:
                    success = DoubleSolveBiCGSTABWithPivoting(out X, A, B);
                    break;

                default:
                    throw new NotImplementedException();
                    //break;
            }

            return success;
        }

        public bool ComplexSolve(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            bool success = false;
            X = null;

            switch (Method)
            {
                case IvyFEMEquationSolverMethod.CG:
                    // エルミート行列の場合はCG
                    // 必要なケースがでてくれば実装する
                    throw new NotImplementedException();
                    //break;

                case IvyFEMEquationSolverMethod.Default:
                case IvyFEMEquationSolverMethod.NoPreconCOCG:
                    success = ComplexSolveNoPreconCOCG(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.COCG:
                    // 複素対称行列の場合はCOCG
                    success = ComplexSolveCOCG(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.ICCOCG:
                    success = ComplexSolveICCOCG(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.NoPreconBiCGSTAB:
                    success = ComplexSolveNoPreconBiCGSTAB(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.BiCGSTAB:
                    success = ComplexSolveBiCGSTAB(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.BiCGSTABWithPivoting:
                    success = ComplexSolveBiCGSTABWithPivoting(out X, A, B);
                    break;

                default:
                    throw new NotImplementedException();
                    //break;
            }

            return success;
        }
    }
}
