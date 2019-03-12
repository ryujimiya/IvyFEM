using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; // DllImport
using LisScalar = System.Numerics.Complex;

namespace IvyFEM.Lis
{
    public class ImportedFunctions
    {
        /////////////////////////////////////////////////////////////////////////
        // Utilities
        [DllImport("liblis.dll")]
        public static extern unsafe int lis_initialize(int* argc, byte** argv);

        [DllImport("liblis.dll")]
        public static extern int lis_finalize();

        /////////////////////////////////////////////////////////////////////////
        // Matrix Operations
        [DllImport("liblis.dll")]
        public static extern unsafe int lis_matrix_create(int comm, NativeLisMatrix** Amat);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_matrix_destroy(NativeLisMatrix* Amat);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_matrix_assemble(NativeLisMatrix* A);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_matrix_set_size(NativeLisMatrix* A, int local_n, int global_n);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_matrix_get_size(NativeLisMatrix* A, int* local_n, int* global_n);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_matrix_get_range(NativeLisMatrix* A, int*@is, int* ie);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_matrix_set_type(NativeLisMatrix* A, MatrixType matrix_type);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_matrix_get_type(NativeLisMatrix* A, MatrixType* matrix_type);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_matrix_set_value(
            SetValueFlag flag, int i, int j, LisScalar value, NativeLisMatrix* A);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_matrix_set_value_new(
            SetValueFlag flag, int i, int j, LisScalar value, NativeLisMatrix* A);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_matrix_set_values(
            SetValueFlag flag, int n, LisScalar* value, NativeLisMatrix* A);

        /////////////////////////////////////////////////////////////////////////
        // Vector Operations
        [DllImport("liblis.dll")]
        public static extern unsafe int lis_vector_create(int comm, NativeLisVector** vec);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_vector_set_size(NativeLisVector* vec, int local_n, int global_n);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_vector_destroy(NativeLisVector* vec);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_vector_get_size(NativeLisVector* v, int* local_n, int* global_n);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_vector_get_range(NativeLisVector* v, int* @is, int* ie);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_vector_get_value(NativeLisVector* v, int i, LisScalar* value);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_vector_get_values(
            NativeLisVector* v, int start, int count, LisScalar* value);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_vector_set_value(
            SetValueFlag flag, int i, LisScalar value, NativeLisVector* v);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_vector_set_values(
            SetValueFlag flag, int count, int[] index, LisScalar* value, NativeLisVector* v);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_vector_set_values2(
            SetValueFlag flag, int start, int count, LisScalar* value, NativeLisVector* v);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_vector_set_all(LisScalar alpha, NativeLisVector* vx);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_vector_conjugate(NativeLisVector* vx);

        /////////////////////////////////////////////////////////////////////////
        // Matrix-Vector Operations
        [DllImport("liblis.dll")]
        public static extern unsafe int lis_matvec(NativeLisMatrix* A, NativeLisVector* x, NativeLisVector* y);

        /////////////////////////////////////////////////////////////////////////
        // Linear Solvers
        [DllImport("liblis.dll")]
        public static extern unsafe int lis_solver_create(NativeLisSolver** solver);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_solver_destroy(NativeLisSolver* solver);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_solver_get_iter(NativeLisSolver* solver, int* iter);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_solver_set_option(byte[] text, NativeLisSolver* solver);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_solver_set_optionC(NativeLisSolver* solver);

        [DllImport("liblis.dll")]
        public static extern unsafe int lis_solve(
            NativeLisMatrix* A, NativeLisVector* b, NativeLisVector* x, NativeLisSolver* solver);
    }
}
