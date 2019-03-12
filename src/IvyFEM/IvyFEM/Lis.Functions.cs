using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using LisScalar = System.Numerics.Complex;

namespace IvyFEM.Lis
{
    public class Functions
    {
        /////////////////////////////////////////////////////////////////////////
        // Utilities
        public static int Initialize()
        {
            string[] args = System.Environment.GetCommandLineArgs();
            int argc = args.Length;

            int ret = 0;
            unsafe
            {
                byte*[] argv = new byte*[argc];
                for (int i = 0; i < argc; i++)
                {
                    byte[] v = System.Text.Encoding.UTF8.GetBytes(args[i]);
                    fixed (byte* vP = &v[0])
                    {
                        argv[i] = vP;
                    }
                }
                fixed (byte** argvP = &argv[0])
                {
                    ret = IvyFEM.Lis.ImportedFunctions.lis_initialize(&argc, argvP);
                }
            }
            return ret;
        }

        public static int Finalize()
        {
            int ret = IvyFEM.Lis.ImportedFunctions.lis_finalize();
            return ret;
        }

        /////////////////////////////////////////////////////////////////////////
        // Matrix Operations
        public static int MatrixCreate(int comm, LisMatrix A)
        {
            int ret = 0;
            unsafe
            {
                A.Native = null;
                fixed (NativeLisMatrix** nativeAPP = &A.Native)
                {
                    ret = IvyFEM.Lis.ImportedFunctions.lis_matrix_create(comm, nativeAPP);
                }
            }
            return ret;
        }

        public static int MatrixDestroy(LisMatrix A)
        {
            int ret = 0;
            unsafe
            {
                ret = IvyFEM.Lis.ImportedFunctions.lis_matrix_destroy(A.Native);
                A.Native = null;
            }
            return ret;
        }

        public static int MatrixAssemble(LisMatrix A)
        {
            int ret = 0;
            unsafe
            {
                ret = IvyFEM.Lis.ImportedFunctions.lis_matrix_assemble(A.Native);
            }
            return ret;
        }

        public static int MatrixSetSize(LisMatrix A, int localN, int globalN)
        {
            int ret = 0;
            unsafe
            {
                ret = IvyFEM.Lis.ImportedFunctions.lis_matrix_set_size(A.Native, localN, globalN);
            }
            return ret;
        }

        public static int MatrixGetSize(LisMatrix A, out int localN, out int globalN)
        {
            int ret = 0;

            localN = 0;
            globalN = 0;
            unsafe
            {
                fixed (int* localNP = &localN)
                fixed (int* globalNP = &globalN)
                {
                    ret = IvyFEM.Lis.ImportedFunctions.lis_matrix_get_size(A.Native, localNP, globalNP);

                }
            }
            return ret;
        }

        public static int MatrixGetRange(LisMatrix A, out int @is, out int ie)
        {
            int ret = 0;

            @is = 0;
            ie = 0;
            unsafe
            {
                fixed (int* isP = &@is)
                fixed (int* ieP = &ie)
                {
                    ret = IvyFEM.Lis.ImportedFunctions.lis_matrix_get_range(A.Native, isP, ieP);
                }
            }
            return ret;
        }

        public static int MatrixSetType(LisMatrix A, MatrixType matrixType)
        {
            int ret = 0;
            unsafe
            {
                IvyFEM.Lis.ImportedFunctions.lis_matrix_set_type(A.Native, matrixType);
            }
            return ret;
        }

        public static int MatrixGetType(LisMatrix A, out MatrixType matrixType)
        {
            int ret = 0;
            unsafe
            {
                fixed (MatrixType* matrixTypeP = &matrixType)
                {
                    ret = IvyFEM.Lis.ImportedFunctions.lis_matrix_get_type(A.Native, matrixTypeP);

                }
            }
            return ret;
        }

        public static int MatrixSetValue(
            SetValueFlag flag, int i, int j, LisScalar value, LisMatrix A)
        {
            int ret = 0;
            unsafe
            {
                ret = IvyFEM.Lis.ImportedFunctions.lis_matrix_set_value(
                        flag, i, j, value, A.Native);
            }
            return ret;
        }

        public static int MatrixSetValueNew(
            SetValueFlag flag, int i, int j, LisScalar value, LisMatrix A)
        {
            int ret = 0;
            unsafe
            {
                ret = IvyFEM.Lis.ImportedFunctions.lis_matrix_set_value_new(
                    flag, i, j, value, A.Native);
            }
            return ret;
        }


        public static int MatrixSetValues(
            SetValueFlag flag, int n, LisScalar[] values, LisMatrix A)
        {
            int ret = 0;
            unsafe
            {
            	fixed (LisScalar* valuesP = &values[0])
            	{
                    ret = IvyFEM.Lis.ImportedFunctions.lis_matrix_set_values(
                        flag, n, valuesP, A.Native);
                }
            }
            return ret;
        }

        /////////////////////////////////////////////////////////////////////////
        // Vector Operations
        public static int VectorCreate(int comm, LisVector v)
        {
            int ret = 0;
            unsafe
            {
                v.Native = null;
                fixed (NativeLisVector ** nativeVPP = &v.Native)
                {
                    ret = IvyFEM.Lis.ImportedFunctions.lis_vector_create(comm, nativeVPP);
                }

            }
            return ret;
        }

        public static int VectorSetSize(LisVector v, int localN, int globalN)
        {
            int ret = 0;
            unsafe
            {
                ret = IvyFEM.Lis.ImportedFunctions.lis_vector_set_size(v.Native, localN, globalN);
            }
            return ret;
        }

        public static int VectorDestroy(LisVector v)
        {
            int ret = 0;
            unsafe
            {
                ret = IvyFEM.Lis.ImportedFunctions.lis_vector_destroy(v.Native);
            }
            return ret;
        }

        public static int VectorGetSize(LisVector v, out int localN, out int globalN)
        {
            int ret = 0;
            localN = 0;
            globalN = 0;
            unsafe
            {
                fixed (int* localNP = &localN)
                fixed (int* globalNP = &globalN)
                {
                    ret = IvyFEM.Lis.ImportedFunctions.lis_vector_get_size(v.Native, localNP, globalNP);
                }
            }
            return ret;
        }

        public static int VectorGetRange(LisVector v, out int @is, out int ie)
        {
            int ret = 0;
            @is = 0;
            ie = 0;
            unsafe
            {
                fixed (int* isP = &@is)
                fixed (int* ieP = &ie)
                {
                    ret = IvyFEM.Lis.ImportedFunctions.lis_vector_get_range(v.Native, isP, ieP);
                }
            }
            return ret;
        }

        public static int VectorGetValue(LisVector v, int i, out LisScalar value)
        {
            int ret = 0;
            value = 0;
            unsafe
            {
                fixed (LisScalar* valueP = &value)
                {
                    ret = IvyFEM.Lis.ImportedFunctions.lis_vector_get_value(v.Native, i, valueP);
                }
            }
            return ret;
        }

        public static int VectorGetValues(LisVector v, int start, int count, LisScalar[] values)
        {
            int ret = 0;
            unsafe
            {
            	fixed (LisScalar* valuesP = &values[0])
            	{
                    ret = IvyFEM.Lis.ImportedFunctions.lis_vector_get_values(v.Native, start, count, valuesP);
                }

            }
            return ret;
        }

        public static int VectorSetValue(SetValueFlag flag, int i, LisScalar value, LisVector v)
        {
            int ret = 0;
            unsafe
            {
                ret = IvyFEM.Lis.ImportedFunctions.lis_vector_set_value(flag, i, value, v.Native);
            }
            return ret;
        }

        public static int VectorSetValues(SetValueFlag flag, int count, int[] index, LisScalar[] values, LisVector v)
        {
            int ret = 0;
            unsafe
            {
            	fixed (LisScalar* valuesP = &values[0])
            	{
                    ret = IvyFEM.Lis.ImportedFunctions.lis_vector_set_values(
                        flag, count, index, valuesP, v.Native);
                }
            }
            return ret;
        }

        public static int VectorSetValues2(SetValueFlag flag, int start, int count, LisScalar[] values, LisVector v)
        {
            int ret = 0;
            unsafe
            {
            	fixed (LisScalar* valuesP = &values[0])
            	{
                    ret = IvyFEM.Lis.ImportedFunctions.lis_vector_set_values2(
                        flag, start, count, valuesP, v.Native);
                }
            }
            return ret;
        }

        public static int VectorSetAll(LisScalar alpha, LisVector v)
        {
            int ret = 0;
            unsafe
            {
                ret = IvyFEM.Lis.ImportedFunctions.lis_vector_set_all(alpha, v.Native);
            }
            return ret;
        }

        public static int VectorConjugate(LisVector v)
        {
            int ret = 0;
            unsafe
            {
                ret = IvyFEM.Lis.ImportedFunctions.lis_vector_conjugate(v.Native);
            }
            return ret;
        }

        /////////////////////////////////////////////////////////////////////////
        // Matrix-Vector Operations
        public static int Matvec(LisMatrix A, LisVector x, LisVector y)
        {
            int ret = 0;

            unsafe
            {
                ret = IvyFEM.Lis.ImportedFunctions.lis_matvec(A.Native, x.Native, y.Native);
            }
            return ret;
        }

        /////////////////////////////////////////////////////////////////////////
        // Linear Solvers
        public static int SolverCreate(LisSolver solver)
        {
            int ret = 0;
            unsafe
            {
                solver.Native = null;
                fixed (NativeLisSolver** nativePP = &solver.Native)
                ret = IvyFEM.Lis.ImportedFunctions.lis_solver_create(nativePP);
            }
            return ret;
        }

        public static int SolverDestroy(LisSolver solver)
        {
            int ret = 0;
            unsafe
            {
                ret = IvyFEM.Lis.ImportedFunctions.lis_solver_destroy(solver.Native);
                solver.Native = null;
            }
            return ret;
        }

        public static int SolverGetIter(LisSolver solver, out int iter)
        {
            int ret = 0;
            iter = 0;
            unsafe
            {
                fixed (int* iterP = &iter)
                {
                    ret = IvyFEM.Lis.ImportedFunctions.lis_solver_get_iter(solver.Native, iterP);
                }
            }
            return ret;
        }

        public static int SolverSetOption(string text, LisSolver solver)
        {
            int ret = 0;
            byte[] textBytes = System.Text.Encoding.UTF8.GetBytes(text);
            unsafe
            {
                ret = IvyFEM.Lis.ImportedFunctions.lis_solver_set_option(textBytes, solver.Native);
            }
            return ret;
        }

        public static int SolverSetOptionC(LisSolver solver)
        {
            int ret = 0;
            unsafe
            {
                ret = IvyFEM.Lis.ImportedFunctions.lis_solver_set_optionC(solver.Native);
            }
            return ret;
        }

        public static int Solve(LisMatrix A, LisVector b, LisVector x, LisSolver solver)
        {
            int ret = 0;
            unsafe
            {
                ret = IvyFEM.Lis.ImportedFunctions.lis_solve(A.Native, b.Native, x.Native, solver.Native);
            }
            return ret;
        }

    }
}
