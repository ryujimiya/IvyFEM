using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public class LisEquationSolver : IEquationSolver
    {
        public LisEquationSolverMethod Method { get; set; } = LisEquationSolverMethod.Default;

        public bool DoubleSolve(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            bool success = false;
            using (IvyFEM.Lis.LisInitializer LisInitializer = new IvyFEM.Lis.LisInitializer())
            using (IvyFEM.Lis.LisMatrix lisA = new IvyFEM.Lis.LisMatrix(A))
            using (IvyFEM.Lis.LisVector lisB = new IvyFEM.Lis.LisVector(B))
            using (IvyFEM.Lis.LisVector lisX = new IvyFEM.Lis.LisVector())
            using (IvyFEM.Lis.LisSolver lisSolver = new IvyFEM.Lis.LisSolver())
            {
                int ret;
                int n = B.Length;
                lisX.SetSize(0, n);
                A = null;
                ret = lisSolver.SetOption(GetOption());
                System.Diagnostics.Debug.Assert(ret == 0);
                ret = lisSolver.SetOptionC();
                System.Diagnostics.Debug.Assert(ret == 0);
                ret = lisSolver.Solve(lisA, lisB, lisX);
                System.Diagnostics.Debug.Assert(ret == 0);
                success = (ret == 0);
                int iter;
                ret = lisSolver.GetIter(out iter);
                System.Diagnostics.Debug.Assert(ret == 0);
                System.Diagnostics.Debug.WriteLine("Lis Solve iter = " + iter);
                System.Numerics.Complex[] complexX = new System.Numerics.Complex[n];
                ret = lisX.GetValues(0, n, complexX);
                System.Diagnostics.Debug.Assert(ret == 0);
                X = new double[n];
                for (int i = 0; i < n; i++)
                {
                    X[i] = complexX[i].Real;
                }
            }
            return success;
        }

        public bool ComplexSolve(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            bool success = false;
            using (IvyFEM.Lis.LisInitializer LisInitializer = new IvyFEM.Lis.LisInitializer())
            using (IvyFEM.Lis.LisMatrix lisA = new IvyFEM.Lis.LisMatrix(A))
            using (IvyFEM.Lis.LisVector lisB = new IvyFEM.Lis.LisVector(B))
            using (IvyFEM.Lis.LisVector lisX = new IvyFEM.Lis.LisVector())
            using (IvyFEM.Lis.LisSolver lisSolver = new IvyFEM.Lis.LisSolver())
            {
                int ret;
                int n = B.Length;
                lisX.SetSize(0, n);
                A = null;
                ret = lisSolver.SetOption(GetOption());
                System.Diagnostics.Debug.Assert(ret == 0);
                ret = lisSolver.SetOptionC();
                System.Diagnostics.Debug.Assert(ret == 0);
                ret = lisSolver.Solve(lisA, lisB, lisX);
                System.Diagnostics.Debug.Assert(ret == 0);
                success = (ret == 0);
                int iter;
                ret = lisSolver.GetIter(out iter);
                System.Diagnostics.Debug.Assert(ret == 0);
                System.Diagnostics.Debug.WriteLine("Lis Solve iter = " + iter);
                X = new System.Numerics.Complex[n];
                ret = lisX.GetValues(0, n, X);
                System.Diagnostics.Debug.Assert(ret == 0);
            }
            return success;
        }

        private string GetOption()
        {
            string option = "-print mem";
            string methodOption = "";
            switch (Method)
            {
                case LisEquationSolverMethod.Default:
                    break;

                case LisEquationSolverMethod.CG:
                    methodOption = "-i cg";
                    break;

                case LisEquationSolverMethod.BiCG:
                    methodOption = "-i bicg";
                    break;

                case LisEquationSolverMethod.CGS:
                    methodOption = "-i cgs";
                    break;

                case LisEquationSolverMethod.BiCGSTAB:
                    methodOption = "-i bicgstab";
                    break;

                case LisEquationSolverMethod.BiCGSTABl:
                    methodOption = "-i bicgstabl";
                    break;

                case LisEquationSolverMethod.GPBiCG:
                    methodOption = "-i gpbicg";
                    break;

                case LisEquationSolverMethod.TFQMR:
                    methodOption = "-i tfqmr";
                    break;

                case LisEquationSolverMethod.Orthominm:
                    methodOption = "-i orthomin";
                    break;

                case LisEquationSolverMethod.GMRESm:
                    methodOption = "-i gmres";
                    break;

                case LisEquationSolverMethod.Jacobi:
                    methodOption = "-i jacobi";
                    break;

                case LisEquationSolverMethod.GaussSeidel:
                    methodOption = "-i gs";
                    break;

                case LisEquationSolverMethod.SOR:
                    methodOption = "-i sor";
                    break;

                case LisEquationSolverMethod.BiCGSafe:
                    methodOption = "-i bicgsafe";
                    break;

                case LisEquationSolverMethod.CR:
                    methodOption = "-i cr";
                    break;

                case LisEquationSolverMethod.BiCR:
                    methodOption = "-i bicr";
                    break;

                case LisEquationSolverMethod.CRS:
                    methodOption = "-i crs";
                    break;

                case LisEquationSolverMethod.BiCRSTAB:
                    methodOption = "-i bicrstab";
                    break;

                case LisEquationSolverMethod.GPBiCR:
                    methodOption = "-i gpbicr";
                    break;

                case LisEquationSolverMethod.BiCRSafe:
                    methodOption = "-i bicrsafe";
                    break;

                case LisEquationSolverMethod.FGMRESm:
                    methodOption = "-i fgmres";
                    break;

                case LisEquationSolverMethod.IDRs:
                    methodOption = "-i idrs";
                    break;

                case LisEquationSolverMethod.IDR1:
                    methodOption = "-i idr1";
                    break;

                case LisEquationSolverMethod.MINRES:
                    methodOption = "-i minres";
                    break;

                case LisEquationSolverMethod.COCG:
                    methodOption = "-i cocg";
                    break;

                case LisEquationSolverMethod.COCR:
                    methodOption = "-i cocr";
                    break;
            }
            if (methodOption != "")
            {
                option += " " + methodOption;
            }
            return option;
        }
    }
}
