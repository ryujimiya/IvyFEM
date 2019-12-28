using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; // DllImport

namespace IvyFEM.Lapack
{
    /// <summary>
    /// liblapacke.dll
    /// プロジェクトターゲットはAnyCPU でなくx64を指定する
    /// （AnyCPUを指定すると呼び出し時System.BadImageFormatExceptionが発生してしまう）
    /// </summary>
    public class ImportedFunctions
    {
        ////////////////////////////////////////////////////////////////
        // BLAS
        ////////////////////////////////////////////////////////////////
        // double
        [DllImport("libblas.dll")]
        public static extern unsafe void daxpy_(int* n, double* da, double[] dx,int* incx,
            double[] dy, int* incy);

        [DllImport("libblas.dll")]
        public static extern unsafe double ddot_(int* n, double[] dx, int* incx, double[] dy, int* incy);

        [DllImport("libblas.dll")]
        public static extern unsafe void dgbmv_(
            byte* trans, int* m, int* n, int* kl, int* ku,
            double* alpha, double[] a, int* lda,
            double[] x, int* incx, 
            double* beta, double[] y, int* incy);

        [DllImport("libblas.dll")]
        public static extern unsafe void dgemm_(
            byte* transa, byte* transb,
            int* m, int* n, int* k,
            double* alpha, double[] a, int* lda,
            double[] b, int* ldb,
            double* beta, double[] c, int* ldc);

        [DllImport("libblas.dll")]
        public static extern unsafe void dgemv_(
            byte* trans, int* m, int* n,
            double* alpha, double[] a, int* lda,
            double[] x, int* incx,
            double* beta, double[] y, int* incy);

        [DllImport("libblas.dll")]
        public static extern unsafe double dnrm2_(int* n, double[] x, int* incx);

        [DllImport("libblas.dll")]
        public static extern unsafe void dscal_(int* n, double* da, double[] dx, int* incx);

        // complex
        [DllImport("libblas.dll")]
        public static extern unsafe void zaxpy_(
            int* n, System.Numerics.Complex* za, System.Numerics.Complex* zx, int* incx,
            System.Numerics.Complex* zy, int* incy);

        // X^H * Y
        [DllImport("libblas.dll")]
        public static extern unsafe System.Numerics.Complex zdotc_(
            int* n, System.Numerics.Complex* zx, int* incx,
            System.Numerics.Complex* zy, int* incy);

        // X^T * Y
        [DllImport("libblas.dll")]
        public static extern unsafe System.Numerics.Complex zdotu_(
            int* n, System.Numerics.Complex* zx, int* incx,
            System.Numerics.Complex* zy, int* incy);

        [DllImport("libblas.dll")]
        public static extern unsafe void zgbmv_(
            byte* trans, int* m, int* n, int* kl, int* ku,
            System.Numerics.Complex* alpha, System.Numerics.Complex* a, int* lda,
            System.Numerics.Complex* x, int* incx, 
            System.Numerics.Complex* beta, System.Numerics.Complex* y, int* incy);

        [DllImport("libblas.dll")]
        public static extern unsafe void zgemm_(
            byte* transa, byte* transb,
            int* m, int* n, int* k,
            System.Numerics.Complex* alpha, System.Numerics.Complex* a, int* lda,
            System.Numerics.Complex* b, int* ldb,
            System.Numerics.Complex* beta, System.Numerics.Complex* c, int* ldc);

        [DllImport("libblas.dll")]
        public static extern unsafe void zgemv_(
            byte* trans, int* m, int* n,
            System.Numerics.Complex* alpha, System.Numerics.Complex* a, int* lda,
            System.Numerics.Complex* x, int* incx, 
            System.Numerics.Complex* beta, System.Numerics.Complex* y, int* incy);

        [DllImport("libblas.dll")]
        public static extern unsafe void zscal_(int* n, System.Numerics.Complex* za,
            System.Numerics.Complex* zx, int* incx);


        ////////////////////////////////////////////////////////////////
        // LAPACK
        ////////////////////////////////////////////////////////////////
        // alternate LAPACKE version exists.
        [DllImport("liblapack.dll")]
        public static extern unsafe void dgesv_(
            int* n, int* nrhs,
            double[] a, int* lda, int[] ipiv,
            double[] b, int* ldb,
            int* info);

        // alternate LAPACKE version exists.
        [DllImport("liblapack.dll")]
        public static extern unsafe void zlacgv_(int* n, System.Numerics.Complex* x, int* incx);

        ////////////////////////////////////////////////////////////////
        // LAPACKE
        ////////////////////////////////////////////////////////////////
        // double
        [DllImport("liblapacke.dll")]
        public static extern int LAPACKE_dgetrf(
            MatrixLayout matrix_layout,
            int m, int n, double[] a,
            int lda, int[] ipiv);

        [DllImport("liblapacke.dll")]
        public static extern int LAPACKE_dgetri(
            MatrixLayout matrix_layout,
            int n, double[] a,
            int lda, int[] ipiv );

        [DllImport("liblapacke.dll")]
        public static extern int LAPACKE_dgbsv(
            MatrixLayout matrix_layout,
            int n, int kl, int ku, int nrhs,
            double[] ab, int ldab, int[] ipiv,
            double[] b, int ldb);

        [DllImport("liblapacke.dll")]
        public static extern int LAPACKE_dgeev(
            MatrixLayout matrix_layout,
            byte jobvl, byte jobvr,
            int n, double[] a, int lda,
            double[] wr, double[] wi,
            double[] vl, int ldvl,
            double[] vr, int ldvr);


        [DllImport("liblapacke.dll")]
        public static extern int LAPACKE_dgesv(
            MatrixLayout matrix_layout,
            int n, int nrhs,
            double[] a, int lda, int[] ipiv,
            double[] b, int ldb);

        [DllImport("liblapacke.dll")]
        public static extern int LAPACKE_dggev(
            MatrixLayout matrix_layout,
            byte jobvl, byte jobvr,
            int n, double[] a, int lda, 
            double[] b, int ldb,
            double[] alphar, double[] alphai, double[] beta,
            double[] vl, int ldvl,
            double[] vr, int ldvr);

        [DllImport("liblapacke.dll")]
        public static extern int LAPACKE_dpbsv(
            MatrixLayout matrix_layout,
            byte uplo, int n, int kd, int nrhs, double[] ab, int ldab,
            double[] b, int ldb);

        [DllImport("liblapacke.dll")]
        public static extern int LAPACKE_dsbev(
            MatrixLayout matrix_layout,
            byte jobz, byte uplo, int n, int kd, double[] ab, int ldab,
            double[] w, double[] z, int ldz);

        [DllImport("liblapacke.dll")]
        public static extern int LAPACKE_dsbgv(
            MatrixLayout matrix_layout,
            byte jobz, byte uplo,
            int n, int ka, int kb,
            double[] ab, int ldab, double[] bb, int ldbb,
            double[] w, double[] z, int ldz);

        [DllImport("liblapacke.dll")]
        public static extern int LAPACKE_dgbtrs(
            MatrixLayout matrix_layout,
            byte trans,
            int n, int kl, int ku, int nrhs,
            double[] ab, int ldab, int[] ipiv, double[] b, int ldb);

        // complex
        [DllImport("liblapacke.dll")]
        public static extern unsafe int LAPACKE_zgetrf(
            MatrixLayout matrix_layout,
            int m, int n, System.Numerics.Complex* a,
            int lda, int[] ipiv);

        [DllImport("liblapacke.dll")]
        public static extern unsafe int LAPACKE_zgetri(
            MatrixLayout matrix_layout,
            int n, System.Numerics.Complex* a,
            int lda, int[] ipiv );

        [DllImport("liblapacke.dll")]
        public static extern unsafe int LAPACKE_zgbsv(
            MatrixLayout matrix_layout,
            int n, int kl, int ku, int nrhs,
            System.Numerics.Complex* ab, int ldab, int[] ipiv,
            System.Numerics.Complex* b, int ldb);

        [DllImport("liblapacke.dll")]
        public static extern unsafe int LAPACKE_zgesv(
            MatrixLayout matrix_layout,
            int n, int nrhs,
            System.Numerics.Complex* a, int lda, int[] ipiv,
            System.Numerics.Complex* b, int ldb);

        [DllImport("liblapacke.dll")]
        public static extern unsafe int LAPACKE_zhbev(
            MatrixLayout matrix_layout,
            byte jobz, byte uplo,
            int n, int kd, System.Numerics.Complex* ab, int ldab,
            double[] w, System.Numerics.Complex* z, int ldz);

        [DllImport("liblapacke.dll")]
        public static extern unsafe int LAPACKE_zhbgv(
            MatrixLayout matrix_layout,
            byte jobz, byte uplo,
            int n, int ka, int kb,
            System.Numerics.Complex* ab, int ldab, System.Numerics.Complex* bb, int ldbb,
            double[] w, System.Numerics.Complex* z, int ldz);

        [DllImport("liblapacke.dll")]
        public static extern unsafe int LAPACKE_zggev(
            MatrixLayout matrix_layout,
            byte jobvl, byte jobvr,
            int n, System.Numerics.Complex* a, int lda,
            System.Numerics.Complex* b, int ldb,
            System.Numerics.Complex* alpha, System.Numerics.Complex* beta,
            System.Numerics.Complex* vl, int ldvl,
            System.Numerics.Complex* vr, int ldvr);

        [DllImport("liblapacke.dll")]
        public static extern unsafe int LAPACKE_zlacgv(int n, System.Numerics.Complex* x, int incx);

        [DllImport("liblapacke.dll")]
        public static extern unsafe int LAPACKE_zpbsv(
            MatrixLayout matrix_layout,
            byte uplo, int n, int kd, int nrhs,
            System.Numerics.Complex* ab, int ldab,
            System.Numerics.Complex* b, int ldb);

        [DllImport("liblapacke.dll")]
        public static extern unsafe int LAPACKE_zgbtrs(
            MatrixLayout matrix_layout,
            byte trans, 
            int n, int kl, int ku, int nrhs,
            System.Numerics.Complex* ab, int ldab, int[] ipiv, System.Numerics.Complex* b, int ldb);
    }
}
