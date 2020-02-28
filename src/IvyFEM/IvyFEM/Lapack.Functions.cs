using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Lapack
{
    public class Functions
    {
        ////////////////////////////////////////////////////////////////
        // Utils
        private static void DoubleBandMatrixLapack2Blas(out double[] blasAB, out int lda,
            double[] AB, int m, int n, int kl, int ku)
        {
            lda = kl + ku + 1;
            blasAB = new double[lda * n];
            for (int col = 0; col < n; col++)
            {
                int k = ku - col;
                for (int row = Math.Max(0, col - ku); row <= Math.Min(m - 1, col + kl); row++)
                {
                    if (!(row >= col - ku && row <= col + kl))
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    blasAB[(k + row) + col * lda] = AB[(row - col) + kl + ku + col * (kl * 2 + ku + 1)];//matrix(I, J)
                }
            }
        }

        private static void ComplexBandMatrixLapack2Blas(out System.Numerics.Complex[] blasAB, out int lda,
            System.Numerics.Complex[] AB, int m, int n, int kl, int ku)
        {
            lda = kl + ku + 1;
            blasAB = new System.Numerics.Complex[lda * n];
            for (int col = 0; col < n; col++)
            {
                int k = ku - col;
                for (int row = Math.Max(0, col - ku); row <= Math.Min(m - 1, col + kl); row++)
                {
                    if (!(row >= col - ku && row <= col + kl))
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    blasAB[(k + row) + col * lda] = AB[(row - col) + kl + ku + col * (kl * 2 + ku + 1)];//matrix(I, J)
                }
            }
        }

        ////////////////////////////////////////////////////////////////
        // BLAS
        ////////////////////////////////////////////////////////////////
        public static double[] daxpy(double a, double[] X, double[] Y)
        {
            double[] Z;

            if (X.Length != Y.Length)
            {
                throw new ArgumentException("Mismatched size: v1(" + X.Length + "), v2(" + Y.Length + ")");
            }

            int n = X.Length;
            int incX = 1;
            int incY = 1;

            Z = new double[n];
            Y.CopyTo(Z, 0);
            unsafe
            {
                IvyFEM.Lapack.ImportedFunctions.daxpy_(&n, &a, X, &incX, Z, &incY);
            }

            return Z;
        }

        public static double ddot(double[] X, double[] Y)
        {
            if (X.Length != Y.Length)
            {
                throw new ArgumentException("Mismatched size: v1(" + X.Length + "), v2(" + Y.Length + ")");
            }

            int n = X.Length;
            int incX = 1;
            int incY = 1;

            double ret = 0;
            unsafe
            {
                ret = IvyFEM.Lapack.ImportedFunctions.ddot_(&n, X, &incX, Y, &incY);
            }
            return ret;
        }

        public static void dgbmvAX(out double[] Y,
            double[] AB, int aRow, int aCol,
            int subdia, int superdia, TransposeType aTransposeType,
            double[] X)
        {
            if (aCol != X.Length)
            {
                throw new ArgumentException("Mismatched size: aCol != X.Length(" + aCol + " != " + X.Length + ")");
            }

            Y = new double[X.Length];

            byte trans = Trans.FromTransposeType(aTransposeType);
            int m = aRow;
            int n = aCol;
            int kl = subdia;
            int ku = superdia;
            double alpha = 1.0;
            //int lda = kl + ku + 1;
            // LDA must be at least (kl + ku + 1).
            int lda;
            int incX = 1;
            double beta = 0.0;
            int incY = 1;

            // BlasのバンドマトリックスとLapackのバンドマトリックスは違う
            double[] blasAB;
            DoubleBandMatrixLapack2Blas(out blasAB, out lda, AB, m, n, kl, ku);

            unsafe
            {
                IvyFEM.Lapack.ImportedFunctions.dgbmv_(
                    &trans, &m, &n, &kl, &ku,
                    &alpha, blasAB, &lda,
                    X, &incX, 
                    &beta, Y, &incY);
            }
        }

        public static void dgemmAB(out double[] C, out int cRow, out int cCol,
            double[]  A, int aRow, int aCol, TransposeType aTransposeType,
            double[]  B, int bRow, int bCol, TransposeType bTransposeType)
        {
            if (aCol != bRow)
            {
                throw new ArgumentException("Mismatched size: aCol != bRow(" + aCol + " != " + bRow + ")");
            }

            cRow = aRow;
            cCol = bCol;

            C = new double[cRow * cCol];

            byte transa = Trans.FromTransposeType(aTransposeType);
            byte transb = Trans.FromTransposeType(bTransposeType);

            int m = aRow;
            int n = bCol;
            int k = aCol;

            double alpha = 1.0;

            int lda = aRow;
            int ldb = bRow;

            double beta = 0.0;
            int ldc = aRow;

            unsafe
            {
                IvyFEM.Lapack.ImportedFunctions.dgemm_(
                    &transa, &transb,
                    &m, &n, &k,
                    &alpha, A, &lda,
                    B, &ldb,
                    &beta, C, &ldc);
            }
        }

        public static void dgemvAX(out double[] Y,
            double[] A, int aRow, int aCol, TransposeType aTransposeType,
            double[] X)
        {
            if (aCol != X.Length)
            {
                throw new ArgumentException("Mismatched size: aCol != X.Length(" + aCol + " != " + X.Length + ")");
            }

            Y = new double[aRow];

            byte trans = Trans.FromTransposeType(aTransposeType);
            int m = aRow;
            int n = aCol;
            double alpha = 1.0;
            int lda = aRow;
            int incX = 1;
            double beta = 0.0;
            int incY = 1;

            unsafe
            {
                IvyFEM.Lapack.ImportedFunctions.dgemv_(&trans, &m, &n, &alpha, A, &lda,
                                        X, &incX, &beta, Y, &incY);
            }
        }

        public static double dnrm2(double[] X)
        {
            int n = X.Length;
            int incX = 1;

            double ret = 0;
            unsafe
            {
                ret = IvyFEM.Lapack.ImportedFunctions.dnrm2_(&n, X, &incX);
            }
            return ret;
        }

        public static double[] dscal(double[] X, double a)
        {
            int n = X.Length;
            int incX = 1;

            double[] Y = new double[n];
            X.CopyTo(Y, 0);

            unsafe
            {
                IvyFEM.Lapack.ImportedFunctions.dscal_(&n, &a, Y, &incX);
            }
            return Y;
        }

        public static System.Numerics.Complex[] zaxpy(
            System.Numerics.Complex a, System.Numerics.Complex[] X, System.Numerics.Complex[] Y)
        {
            System.Numerics.Complex[] Z;
            if (X.Length != Y.Length)
            {
                throw new ArgumentException("Mismatched size: v1(" + X.Length + "), v2(" + Y.Length + ")");
            }

            int n = X.Length;
            int incX = 1;
            int incY = 1;

            Z = new System.Numerics.Complex[n];
            Y.CopyTo(Z, 0);
            unsafe
            {
                fixed (System.Numerics.Complex *XP = &X[0])
                fixed (System.Numerics.Complex* ZP = &Z[0])
                {
                    IvyFEM.Lapack.ImportedFunctions.zaxpy_(&n, &a, XP, &incX, ZP, &incY);
                }
            }
            return Z;
        }

        // X^H * Y
        public static System.Numerics.Complex zdotc(System.Numerics.Complex[] X, System.Numerics.Complex[] Y)
        {
            if (X.Length != Y.Length)
            {
                throw new ArgumentException("Mismatched size: v1(" + X.Length + "), v2(" + Y.Length + ")");
            }

            int n = X.Length;
            int incX = 1;
            int incY = 1;

            System.Numerics.Complex ret = 0;
            unsafe
            {
                fixed (System.Numerics.Complex* XP = &X[0])
                fixed (System.Numerics.Complex* YP = &Y[0])
                {
                    ret = IvyFEM.Lapack.ImportedFunctions.zdotc_(&n, XP, &incX, YP, &incY);
                }
            }
            return ret;
        }

        // X^T * Y
        public static System.Numerics.Complex zdotu(System.Numerics.Complex[] X, System.Numerics.Complex[] Y)
        {
            if (X.Length != Y.Length)
            {
                throw new ArgumentException("Mismatched size: v1(" + X.Length + "), v2(" + Y.Length + ")");
            }

            int n = X.Length;
            int incX = 1;
            int incY = 1;

            System.Numerics.Complex ret = 0;
            unsafe
            {
                fixed (System.Numerics.Complex* XP = &X[0])
                fixed (System.Numerics.Complex* YP = &Y[0])
                {
                    ret = IvyFEM.Lapack.ImportedFunctions.zdotu_(&n, XP, &incX, YP, &incY);
                }
            }
            return ret;
        }

        public static void zgbmvAX(out System.Numerics.Complex[] Y,
            System.Numerics.Complex[] AB, int aRow, int aCol,
            int subdia, int superdia, TransposeType aTransposeType,
            System.Numerics.Complex[] X)
        {
            if (aCol != X.Length)
            {
                throw new ArgumentException("Mismatched size: aCol != X.Length(" + aCol + " != " + X.Length + ")");
            }

            Y = new System.Numerics.Complex[X.Length];

            byte trans = Trans.FromTransposeType(aTransposeType);
            int m = aRow;
            int n = aCol;
            int kl = subdia;
            int ku = superdia;
            System.Numerics.Complex alpha = 1.0;
            //int lda = kl + ku + 1;
            // LDA must be at least (kl + ku + 1).
            int lda;
            int incX = 1;
            System.Numerics.Complex beta = 0.0;
            int incY = 1;

            // BlasのバンドマトリックスとLapackのバンドマトリックスは違う
            System.Numerics.Complex[] blasAB;
            ComplexBandMatrixLapack2Blas(out blasAB, out lda, AB, m, n, kl, ku);

            unsafe
            {
                fixed (System.Numerics.Complex* blasABP = &blasAB[0])
                fixed (System.Numerics.Complex* XP = &X[0])
                fixed (System.Numerics.Complex* YP = &Y[0])
                {
                    IvyFEM.Lapack.ImportedFunctions.zgbmv_(
                    &trans, &m, &n, &kl, &ku,
                    &alpha, blasABP, &lda,
                    XP, &incX, &beta,
                    YP, &incY);
                }
            }
        }

        public static int zgemmAB(out System.Numerics.Complex[] C, out int cRow, out int cCol,
            System.Numerics.Complex[] A, int aRow, int aCol, TransposeType aTransposeType,
            System.Numerics.Complex[] B, int bRow, int bCol, TransposeType bTransposeType)
        {
            if (aCol != bRow)
            {
                throw new ArgumentException("Mismatched size: aCol != bRow(" + aCol + " != " + bRow + ")");
            }

            cRow = aRow;
            cCol = bCol;

            C = new System.Numerics.Complex[cRow * cCol];

            byte transa = Trans.FromTransposeType(aTransposeType);
            byte transb = Trans.FromTransposeType(bTransposeType);
            int m = aRow;
            int n = bCol;
            int k = aCol;

            System.Numerics.Complex alpha = 1.0;

            int lda = aRow;
            int ldb = bRow;

            System.Numerics.Complex beta = 0.0;
            int ldc = aRow;

            unsafe
            {
                fixed (System.Numerics.Complex* AP = &A[0])
                fixed (System.Numerics.Complex* BP = &B[0])
                fixed (System.Numerics.Complex* CP = &C[0])
                {
                    IvyFEM.Lapack.ImportedFunctions.zgemm_(
                        &transa, &transb,
                        &m, &n, &k,
                        &alpha, AP, &lda,
                        BP, &ldb,
                        &beta, CP, &ldc);
                }
            }
            return 0;
        }

        public static void zgemvAX(out System.Numerics.Complex[] Y,
            System.Numerics.Complex[] A, int aRow, int aCol, TransposeType aTransposeType,
            System.Numerics.Complex[] X)
        {
            if (aCol != X.Length)
            {
                throw new ArgumentException("Mismatched size: aCol != X.Length(" + aCol + " != " + X.Length + ")");
            }

            Y = new System.Numerics.Complex[aRow];

            byte trans = Trans.FromTransposeType(aTransposeType);
            int m = aRow;
            int n = aCol;
            System.Numerics.Complex alpha = 1.0;
            int lda = aRow;
            int incX = 1;
            System.Numerics.Complex beta = 0.0;
            int incY = 1;

            unsafe
            {
                fixed (System.Numerics.Complex* AP = &A[0])
                fixed (System.Numerics.Complex* XP = &X[0])
                fixed (System.Numerics.Complex* YP = &Y[0])
                {
                    IvyFEM.Lapack.ImportedFunctions.zgemv_(&trans, &m, &n, &alpha, AP, &lda,
                                        XP, &incX, &beta, YP, &incY);
                }
            }
        }

        public static System.Numerics.Complex[] zscal(System.Numerics.Complex[] X, System.Numerics.Complex a)
        {
            int n = X.Length;
            int incX = 1;

            System.Numerics.Complex[] Y = new System.Numerics.Complex[n];
            X.CopyTo(Y, 0);

            unsafe
            {
                fixed (System.Numerics.Complex* YP = &Y[0])
                {
                    IvyFEM.Lapack.ImportedFunctions.zscal_(&n, &a, YP, &incX);
                }
            }
            return Y;
        }

        ////////////////////////////////////////////////////////////////
        // LAPACK
        ////////////////////////////////////////////////////////////////


        ////////////////////////////////////////////////////////////////
        // LAPACKE
        ////////////////////////////////////////////////////////////////
        public static int dgetrf(out double[] X, out int[] piv, double[] A, int aRow, int aCol)
        {
            int m = aRow;
            int n = aCol;

            X = new double[A.Length];
            A.CopyTo(X, 0);
            // A = P*L*U; the unit diagonal elements of L are not stored.

            int lda = m;
            int[] ipiv = new int[Math.Min(m, n)];
            // The pivot indices; for 1 <= i <= min(M,N), 
            // row i of the matrix was interchanged with row IPIV(i).

            int ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_dgetrf(
                MatrixLayout.ColMajor,
                m, n, X, lda, ipiv);
            if (ret != 0)
            {
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            // ipiv: 1 <= i <= min(M,N), row i
            piv = new int[ipiv.Length];
            for (int i = 0; i < piv.Length; i++)
            {
                piv[i] = ipiv[i] - 1;
            }
            return ret;
        }

        public static int dgetri(out double[] X, int[] piv, double[] A, int aRow, int aCol)
        {
            if (aRow != aCol)
            {
                throw new ArgumentException();
            }
            int n = aRow;

            X = new double[A.Length];
            A.CopyTo(X, 0);
            // A = P*L*U; the unit diagonal elements of L are not stored.

            int lda = n;
            int[] ipiv = new int[piv.Length];
            // The pivot indices; for 1 <= i <= min(M,N), 
            // row i of the matrix was interchanged with row IPIV(i).
            // ipiv: 1 <= i <= min(M,N), row i
            for (int i = 0; i < piv.Length; i++)
            {
                ipiv[i] = piv[i] + 1;
            }

            int ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_dgetri(
                MatrixLayout.ColMajor,
                n, X, lda, ipiv);
            if (ret != 0)
            {
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            return ret;
        }

        public static int dgbsv(out double[] X, out int xRow, out int xCol,
            double[]  A, int aRow, int aCol, int subdia, int superdia,
            double[]  B, int bRow, int bCol)
        {
            double[] dirtyA = new double[A.Length];
            A.CopyTo(dirtyA, 0);

            int[] ipiv;
            int ret = dgbsv_dirty(out X, out xRow, out xCol,
                dirtyA, aRow, aCol, subdia, superdia,
                B, bRow, bCol, out ipiv);
            return ret;
        }

        public static int dgbsv_dirty(out double[] X, out int xRow, out int xCol,
            double[] dirtyA, int aRow, int aCol, int subdia, int superdia,
            double[] B, int bRow, int bCol, out int[] ipiv)
        {
            // DOUBLE PRECISION array, dimension (LDAB,N)
            // On entry, the matrix A in band storage, in rows KL+1 to 2*KL+KU+1; rows 1 to KL of the array need not be set.
            // The j-th column of A is stored in the j-th column of the array AB as follows:
            // AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)

            // DOUBLE PRECISION array, dimension (LDB,NRHS)
            // On entry, the N-by-NRHS right hand side matrix B.

            int n = aRow;
            // The number of linear equations, i.e., the order of the matrix A.  N >= 0.

            int kl = subdia;
            // The number of subdiagonals within the band of A.  KL >= 0.

            int ku = superdia;
            // The number of superdiagonals within the band of A.  KU >= 0.

            int nrhs = bCol;
            // The number of right hand sides, i.e., the number of columns

            int lda = 2 * subdia + superdia + 1;
            // The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.

            ipiv = new int[n];

            int ldb = bRow;

            double[] C = new double[B.Length];
            B.CopyTo(C, 0);

            int ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_dgbsv(
                MatrixLayout.ColMajor,
                n, kl, ku, nrhs,
                dirtyA, lda, ipiv,
                C, ldb);
            if (ret != 0)
            {
                X = null;
                xRow = 0;
                xCol = 0;
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            X = C;
            xRow = bRow;
            xCol = bCol;
            return ret;
        }

        public static int dgeev(double[] A, int aRow, int aCol,
            out System.Numerics.Complex[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            double[] dirtyA = new double[A.Length];
            A.CopyTo(dirtyA, 0);
            int ret = dgeev_dirty(dirtyA, aRow, aCol, out eVals, out eVecs);
            return ret;
        }

        public static int dgeev_dirty(double[] dirtyA, int aRow, int aCol,
            out System.Numerics.Complex[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            byte jobvl = Job.DontCompute;
            byte jobvr = Job.Compute;
            int n = aCol;
            int lda = n;

            double[] wr = new double[n];
            double[] wi = new double[n];

            int ldvl = 1;
            double[] vl = null;

            int ldvr = n;
            double[] vr = new double[ldvr * n];

            int ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_dgeev(
                MatrixLayout.ColMajor,
                jobvl, jobvr,
                n, dirtyA, lda,
                wr, wi,
                vl, ldvl,
                vr, ldvr);

            if (ret != 0)
            {
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            // 固有値を格納
            eVals = new System.Numerics.Complex[n];
            for (int i = 0; i < n; i++)
            {
                eVals[i] = new System.Numerics.Complex(wr[i], wi[i]);
            }

            // 固有ベクトルを格納
            eVecs = new System.Numerics.Complex[n][];
            for (int i = 0; i < n; i++)
            {
                if (Math.Abs(wi[i]) < IvyFEM.Constants.PrecisionLowerLimit)
                {
                    // 実数の固有ベクトル
                    eVecs[i] = new System.Numerics.Complex[ldvr];
                    for (int j = 0; j < ldvr; j++)
                    {
                        eVecs[i][j] = new System.Numerics.Complex(vr[i * ldvr + j], 0.0);
                    }
                }
                else
                {
                    // 複素数（複素共役）の固有ベクトル
                    var vec1 = new System.Numerics.Complex[ldvr];
                    var vec2 = new System.Numerics.Complex[ldvr];

                    for (int j = 0; j < ldvr; j++)
                    {
                        vec1[j] = new System.Numerics.Complex(vr[i * ldvr + j], vr[(i + 1) * ldvr + j]);
                        vec2[j] = new System.Numerics.Complex(vr[i * ldvr + j], -vr[(i + 1) * ldvr + j]);
                    }
                    eVecs[i] = vec1;
                    eVecs[i + 1] = vec2;

                    // 2列参照するので
                    i++;
                }
            }

            return ret;
        }

        public static int dgesv(out double[] X, out int xRow, out int xCol,
            double[] A, int aRow, int aCol,
            double[] B, int bRow, int bCol)
        {
            X = null;
            xRow = 0;
            xCol = 0;

            int n = aRow;
            int nrhs = bCol;
            int lda = n;
            int[] ipiv = new int[n];
            int ldb = bRow;
            double[] C = new double[B.Length];
            B.CopyTo(C, 0);

            double[] dirtyA = new double[A.Length];
            A.CopyTo(dirtyA, 0);

            // LAPACKE
            int ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_dgesv(
                MatrixLayout.ColMajor, n, nrhs, dirtyA, lda, ipiv, C, ldb);
            /*
            // LAPACK
            int ret = 0;
            unsafe
            {
                IvyFEM.Lapack.ImportedFunctions.dgesv_(
                    &n, &nrhs,
                    dirtyA, &lda, ipiv,
                    B, &ldb,
                    &ret);
            }
            */
            if (ret != 0)
            {
                X = null;
                xRow = 0;
                xCol = 0;
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            X = C;
            xRow = bRow;
            xCol = bCol;
            return ret;
        }

        public static int dggev(double[] A, int aRow, int aCol,
            double[] B, int bRow, int bCol,
            out System.Numerics.Complex[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            double[] dirtyA = new double[A.Length];
            A.CopyTo(dirtyA, 0);
            double[] dirtyB = new double[B.Length];
            B.CopyTo(dirtyB, 0);
            int ret = dggev_dirty(dirtyA, aRow, aCol, dirtyB, bRow, bCol, out eVals, out eVecs);
            return ret;
        }

        // A,Bは書き換わる
        public static int dggev_dirty(double[] dirtyA, int aRow, int aCol,
            double[] dirtyB, int bRow, int bCol,
            out System.Numerics.Complex[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            byte jobvl = Job.DontCompute;
            byte jobvr = Job.Compute;
            int n = aCol;
            int lda = n;
            int ldb = n;

            // 計算された固有値の実部の分子部分が入る．
            double[] alphar = new double[n];
            // 計算された固有値の虚部の分子部分が入る．
            // 複素共役対の場合は，alphai[j]=(正値)，alphai[j+1]=(負値) の順に入る．
            double[] alphai = new double[n];
            // 計算された固有値の分母部分が入る
            //  (alphar[j] + alphai[j] * i) / beta[j] (i:虚数単位)が一般化固有値となる
            double[] beta = new double[n];

            int ldvl = 1;
            double[] vl = null;

            int ldvr = n;
            double[] vr = new double[ldvr * n];

            int ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_dggev(
                MatrixLayout.ColMajor, jobvl, jobvr,
                n, dirtyA, lda, dirtyB, ldb,
                alphar, alphai, beta, vl, ldvl, vr, ldvr);
            if (ret != 0)
            {
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            eVals = new System.Numerics.Complex[n];
            for (int i = 0; i < n; i++)
            {
                if (Math.Abs(beta[i]) < IvyFEM.Constants.PrecisionLowerLimit)
                {
                    double real = ((Math.Abs(alphar[i]) < IvyFEM.Constants.PrecisionLowerLimit) ?
                        (double.NaN)
                        : ((Math.Abs(alphar[i]) > 0) ? (double.PositiveInfinity) : (double.NegativeInfinity)));
                    double imag = ((Math.Abs(alphai[i]) < IvyFEM.Constants.PrecisionLowerLimit) ?
                        (double.NaN)
                        : ((Math.Abs(alphai[i]) > 0) ? (double.PositiveInfinity) : (double.NegativeInfinity)));
                    eVals[i] = new System.Numerics.Complex(real, imag);
                }
                else
                {
                    double real = ((Math.Abs(alphar[i]) < IvyFEM.Constants.PrecisionLowerLimit) ?
                        (0.0) : (alphar[i] / beta[i]));
                    double imag = ((Math.Abs(alphai[i]) < IvyFEM.Constants.PrecisionLowerLimit) ?
                        (0.0) : (alphai[i] / beta[i]));
                    eVals[i] = new System.Numerics.Complex(real, imag);
                }
            }

            eVecs = new System.Numerics.Complex[n][];
            for (int i = 0; i < n; i++)
            {
                if (Math.Abs(alphai[i]) < IvyFEM.Constants.PrecisionLowerLimit)
                {
                    // 実数の固有ベクトル
                    eVecs[i] = new System.Numerics.Complex[ldvr];

                    for (int j = 0; j < ldvr; j++)
                    {
                        eVecs[i][j] = new System.Numerics.Complex(vr[i * ldvr + j], 0.0);
                    }
                }
                else
                {
                    // 複素数（複素共役）の固有ベクトル
                    var vec1 = new System.Numerics.Complex[ldvr];
                    var vec2 = new System.Numerics.Complex[ldvr];

                    for (int j = 0; j < ldvr; j++)
                    {
                        vec1[j] = new System.Numerics.Complex(vr[i * ldvr + j], vr[(i + 1) * ldvr + j]);
                        vec2[j] = new System.Numerics.Complex(vr[i * ldvr + j], -vr[(i + 1) * ldvr + j]);
                    }
                    eVecs[i] = vec1;
                    eVecs[i + 1] = vec2;

                    i++;
                }
            }

            return ret;
        }

        public static int dpbsv(out double[] X, out int xRow, out int xCol,
            double[] A, int aRow, int aCol, int superdia,
            double[] B, int bRow, int bCol)
        {
            //if UPLO = 'U', AB(KD + 1 + i - j, j) = A(i, j) for max(1, j - KD) <= i <= j;
            byte uplo = UpperLower.Upper;
        
            int n = aRow;

            int kd = superdia;
            //The number of superdiagonals of the matrix A if UPLO = 'U'

            int nrhs = bCol;
            // The number of right hand sides, i.e., the number of columns

            int lda = kd + 1;
            // The leading dimension of the array AB.LDAB >= KD + 1.

            int ldb = bRow;

            double[] C = new double[B.Length];
            B.CopyTo(C, 0);

            double[] dirtyA = new double[A.Length];
            A.CopyTo(dirtyA, 0);

            int ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_dpbsv(
                MatrixLayout.ColMajor,
                uplo, n, kd, nrhs, dirtyA, lda,
                C, ldb);
            if (ret != 0)
            {
                X = null;
                xRow = 0;
                xCol = 0;
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            X = C;
            xRow = bRow;
            xCol = bCol;
            return ret;
        }

        public static int dsbev(double[] A, int aRow, int aCol, int superDia,
            out double[] eVals, out double[][] eVecs)
        {
            byte jobz = Job.Compute;
            byte uplo = UpperLower.Upper;
            int n = aCol;
            int kd = superDia;
            int ldab = kd + 1;
            // The leading dimension of the array AB.  LDAB >= KD + 1.

            double[] w = new double[n];
            int ldz = n;
            double[] z = new double[ldz * n];

            double[] dirtyA = new double[A.Length];
            A.CopyTo(dirtyA, 0);

            int ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_dsbev(
                MatrixLayout.ColMajor,
                jobz, uplo, n, kd, dirtyA, ldab,
                w, z, ldz);
            if (ret != 0)
            {
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            eVals = new double[n];
            for (int i = 0; i < n; i++)
            {
                eVals[i] = w[i];
            }

            eVecs = new double[n][];
            for (int i = 0; i < n; i++)
            {
                eVecs[i] = new double[ldz];
                for (int j = 0; j < ldz; j++)
                {
                    eVecs[i][j] = z[i * ldz + j];
                }
            }
            return ret;
        }

        public static int dsbgv(double[] A, int aRow, int aCol, int aSuperdia,
            double[] B, int bRow, int bCol, int bSuperdia,
            out double[] eVals, out double[][]  eVecs)
        {
            byte jobz = Job.Compute;
            byte uplo = UpperLower.Upper;
            int n = aCol;
            int ka = aSuperdia;
            int kb = bSuperdia;
            int ldab = ka + 1;
            // the leading dimension of the array AB. ldab >= ka + 1.
            int ldbb = kb + 1;
            // the leading dimension of the array BB. ldbb >= kb + 1.

            double[] w = new double[n];
            int ldz = n;
            // if jobz == 'V' then N <= ldz
            double[] z = new double[ldz * n];

            double[] dirtyA = new double[A.Length];
            A.CopyTo(dirtyA, 0);
            double[] dirtyB = new double[B.Length];
            B.CopyTo(dirtyB, 0);

            int ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_dsbgv(
                MatrixLayout.ColMajor, jobz, uplo,
                n, ka, kb, dirtyA, ldab, dirtyB, ldbb,
                w, z, ldz);
            if (ret != 0)
            {
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }
            eVals = new double[n];
            for (int i = 0; i < n; i++)
            {
                eVals[i] = w[i];
            }

            eVecs = new double[n][];
            for (int i = 0; i < n; i++)
            {
                eVecs[i] = new double[ldz];
                for (int j = 0; j < ldz; j++)
                {
                    eVecs[i][j] = z[i * ldz + j];
                }
            }
            return ret;
        }

        public static int dgbtrs(out double[] X, out int xRow, out int xCol,
            int[] ipiv,
            double[] A, int aRow, int aCol, int subdia, int superdia,
            double[] B, int bRow, int bCol)
        {
            byte trans = Trans.Nop;
            int n = aRow;
            int kl = subdia;
            int ku = superdia;
            int nrhs = bCol;
            // The number of right hand sides, i.e., the number of columns
            // of the matrix B
            int ldab = 2 * kl + ku + 1;
            int ldb = n;
            double[] C = new double[B.Length];
            B.CopyTo(C, 0);

            int ret = ImportedFunctions.LAPACKE_dgbtrs(
                MatrixLayout.ColMajor,
                trans,
                n, kl, ku, nrhs,
                A, ldab, ipiv, C, ldb);
            if (ret != 0)
            {
                X = null;
                xRow = 0;
                xCol = 0;
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            X = C;
            xRow = bRow;
            xCol = bCol;
            return ret;
        }

        public static int zgetrf(out System.Numerics.Complex[] X, out int[] piv,
            System.Numerics.Complex[] A, int aRow, int aCol)
        {
            int m = aRow;
            int n = aCol;

            X = new System.Numerics.Complex[A.Length];
            A.CopyTo(X, 0);
            // A = P*L*U; the unit diagonal elements of L are not stored.

            int lda = m;
            int[] ipiv = new int[Math.Min(m, n)];
            // The pivot indices; for 1 <= i <= min(M,N), 
            // row i of the matrix was interchanged with row IPIV(i).

            int ret = 0;
            unsafe
            {
                fixed (System.Numerics.Complex* XP = &X[0])
                {
                    ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_zgetrf(
                        MatrixLayout.ColMajor,
                        m, n, XP, lda, ipiv);
                }
            }
            if (ret != 0)
            {
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            // ipiv: 1 <= i <= min(M,N), row i
            piv = new int[ipiv.Length];
            for (int i = 0; i < piv.Length; i++)
            {
                piv[i] = ipiv[i] - 1;
            }
            return ret;
        }

        public static int zgetri(
            out System.Numerics.Complex[] X, int[] piv, System.Numerics.Complex[] A, int aRow, int aCol)
        {
            if (aRow != aCol)
            {
                throw new ArgumentException();
            }
            int n = aRow;

            X = new System.Numerics.Complex[A.Length];
            A.CopyTo(X, 0);
            // A = P*L*U; the unit diagonal elements of L are not stored.

            int lda = n;
            int[] ipiv = new int[piv.Length];
            // The pivot indices; for 1 <= i <= min(M,N), 
            // row i of the matrix was interchanged with row IPIV(i).
            // ipiv: 1 <= i <= min(M,N), row i
            for (int i = 0; i < piv.Length; i++)
            {
                ipiv[i] = piv[i] + 1;
            }

            int ret = 0;
            unsafe
            {
                fixed (System.Numerics.Complex* XP = &X[0])
                {
                    ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_zgetri(
                        MatrixLayout.ColMajor,
                        n, XP, lda, ipiv);
                }
            }
            if (ret != 0)
            {
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            return ret;
        }

        public static int zgbsv(out System.Numerics.Complex[] X, out int xRow, out int xCol,
            System.Numerics.Complex[] A, int aRow, int aCol, int subdia, int superdia,
            System.Numerics.Complex[] B, int bRow, int bCol)
        {
            System.Numerics.Complex[] dirtyA = new System.Numerics.Complex[A.Length];
            A.CopyTo(dirtyA, 0);

            int[] ipiv;
            int ret = zgbsv_dirty(out X, out xRow, out xCol,
                dirtyA, aRow, aCol, subdia, superdia,
                B, bRow, bCol, out ipiv);
            return ret;
        }

        public static int zgbsv_dirty(out System.Numerics.Complex[] X, out int xRow, out int xCol,
            System.Numerics.Complex[] dirtyA, int aRow, int aCol, int subdia, int superdia,
            System.Numerics.Complex[] B, int bRow, int bCol, out int[] ipiv)
        {
            // COMPLEX*16 array, dimension (LDAB,N)
            // On entry, the matrix A in band storage, in rows KL+1 to 2*KL+KU+1; rows 1 to KL of the array need not be set.
            // The j-th column of A is stored in the j-th column of the array AB as follows:
            // AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)

            // COMPLEX*16 array, dimension (LDB,NRHS)
            // On entry, the N-by-NRHS right hand side matrix B.

            int n = aRow;
            // The number of linear equations, i.e., the order of the matrix A.  N >= 0.

            int kl = subdia;
            // The number of subdiagonals within the band of A.  KL >= 0.

            int ku = superdia;
            // The number of superdiagonals within the band of A.  KU >= 0.

            int nrhs = bCol;
            // The number of right hand sides, i.e., the number of columns

            int lda = 2 * subdia + superdia + 1;
            // The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.

            ipiv = new int[n];

            int ldb = bRow;

            System.Numerics.Complex[] C = new System.Numerics.Complex[B.Length];
            B.CopyTo(C, 0);

            int ret;
            unsafe
            {
                fixed (System.Numerics.Complex* dirtyAP = &dirtyA[0])
                fixed (System.Numerics.Complex* CP = &C[0])
                {
                    ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_zgbsv(
                        MatrixLayout.ColMajor,
                        n, kl, ku, nrhs, dirtyAP, lda, ipiv, CP, ldb);

                }
            }
            if (ret != 0)
            {
                X = null;
                xRow = 0;
                xCol = 0;
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            X = C;
            xRow = bRow;
            xCol = bCol;

            return ret;
        }

        public static int zgesv(out System.Numerics.Complex[] X, out int xRow, out int xCol,
            System.Numerics.Complex[] A, int aRow, int aCol,
            System.Numerics.Complex[] B, int bRow, int bCol)
        {
            int n = aRow;
            int nrhs = bCol;
            int lda = n;
            int[] ipiv = new int[n];
            int ldb = bRow;
            System.Numerics.Complex[] C = new System.Numerics.Complex[B.Length];
            B.CopyTo(C, 0);

            System.Numerics.Complex[] dirtyA = new System.Numerics.Complex[A.Length];
            A.CopyTo(dirtyA, 0);

            int ret = 0;
            unsafe
            {
                fixed (System.Numerics.Complex* dirtyAP = &dirtyA[0])
                fixed (System.Numerics.Complex* CP = &C[0])
                {
                    ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_zgesv(
                        MatrixLayout.ColMajor,
                        n, nrhs, dirtyAP, lda, ipiv, CP, ldb);
                }
            }

            if (ret != 0)
            {
                X = null;
                xRow = 0;
                xCol = 0;
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            X = C;
            xRow = bRow;
            xCol = bCol;

            return ret;
        }

        public static int zhbev(System.Numerics.Complex[] A, int aRow, int aCol, int superDia,
            out double[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            byte jobz = Job.Compute;
            byte uplo = UpperLower.Upper;
            int n = aCol;
            int kd = superDia;
            int ldab = kd + 1;
            // The leading dimension of the array AB.  LDAB >= KD + 1.

            double[] w = new double[n];
            int ldz = n;
            System.Numerics.Complex[] z = new System.Numerics.Complex[ldz * n];

            System.Numerics.Complex[] dirtyA = new System.Numerics.Complex[A.Length];
            A.CopyTo(dirtyA, 0);

            int ret = 0;
            unsafe
            {
                fixed (System.Numerics.Complex* dirtyAP = &dirtyA[0])
                fixed (System.Numerics.Complex* zP = &z[0])
                {
                    ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_zhbev(
                    MatrixLayout.ColMajor,
                    jobz, uplo, n, kd, dirtyAP, ldab,
                    w, zP, ldz);
                }
            }
            if (ret != 0)
            {
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            eVals = new double[n];
            for (int i = 0; i < n; i++)
            {
                eVals[i] = w[i];
            }

            eVecs = new System.Numerics.Complex[n][];
            for (int i = 0; i < n; i++)
            {
                eVecs[i] = new System.Numerics.Complex[ldz];
                for (int j = 0; j < ldz; j++)
                {
                    eVecs[i][j] = z[i * ldz + j];
                }
            }
            return ret;
        }

        public static int zhbgv(System.Numerics.Complex[] A, int aRow, int aCol, int aSuperdia,
            System.Numerics.Complex[] B, int bRow, int bCol, int bSuperdia,
            out double[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            System.Numerics.Complex[] dirtyA = new System.Numerics.Complex[A.Length];
            A.CopyTo(dirtyA, 0);
            System.Numerics.Complex[] dirtyB = new System.Numerics.Complex[B.Length];
            B.CopyTo(dirtyB, 0);
            int ret = zhbgv_dirty(dirtyA, aRow, aCol, aSuperdia, dirtyB, bRow, bCol, bSuperdia, out eVals, out eVecs);
            return ret;
        }

        // A,Bは書き換わる
        public static int zhbgv_dirty(System.Numerics.Complex[] dirtyA, int aRow, int aCol, int aSuperdia,
            System.Numerics.Complex[] dirtyB, int bRow, int bCol, int bSuperdia,
            out double[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            byte jobz = Job.Compute;
            byte uplo = UpperLower.Upper;
            int n = aCol;
            int ka = aSuperdia;
            int kb = bSuperdia;
            int ldab = ka + 1;
            // the leading dimension of the array AB. ldab >= ka + 1.
            int ldbb = kb + 1;
            // the leading dimension of the array BB. ldbb >= kb + 1.

            double[] w = new double[n];
            int ldz = n;
            // if jobz == 'V' then N <= ldz
            System.Numerics.Complex[] z = new System.Numerics.Complex[ldz * n];

            int ret;
            unsafe
            {
                fixed (System.Numerics.Complex* dirtyAP = &dirtyA[0])
                fixed (System.Numerics.Complex* dirtyBP = &dirtyB[0])
                fixed (System.Numerics.Complex* zP = &z[0])
                {
                    ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_zhbgv(
                        MatrixLayout.ColMajor, jobz, uplo,
                        n, ka, kb, dirtyAP, ldab, dirtyBP, ldbb,
                        w, zP, ldz);

                }
            }
            if (ret != 0)
            {
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }
            eVals = new double[n];
            for (int i = 0; i < n; i++)
            {
                eVals[i] = w[i];
            }

            eVecs = new System.Numerics.Complex[n][];
            for (int i = 0; i < n; i++)
            {
                eVecs[i] = new System.Numerics.Complex[ldz];
                for (int j = 0; j < ldz; j++)
                {
                    eVecs[i][j] = z[i * ldz + j];
                }
            }
            return ret;
        }

        public static int zgeev(System.Numerics.Complex[] A, int aRow, int aCol,
            out System.Numerics.Complex[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            System.Numerics.Complex[] dirtyA = new System.Numerics.Complex[A.Length];
            A.CopyTo(dirtyA, 0);
            int ret = zgeev_dirty(dirtyA, aRow, aCol, out eVals, out eVecs);
            return ret;
        }

        public static int zgeev_dirty(System.Numerics.Complex[] dirtyA, int aRow, int aCol,
            out System.Numerics.Complex[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            byte jobvl = Job.DontCompute;
            byte jobvr = Job.Compute;
            int n = aCol;
            int lda = n;

            System.Numerics.Complex[] w = new System.Numerics.Complex[n];

            int ldvl = 1;
            System.Numerics.Complex[] vl = new System.Numerics.Complex[ldvl * n];

            int ldvr = n;
            System.Numerics.Complex[] vr = new System.Numerics.Complex[ldvr * n];

            int ret = -1;
            unsafe
            {
                fixed (System.Numerics.Complex* dirtyAP = &dirtyA[0])
                fixed (System.Numerics.Complex* wP = &w[0])
                fixed (System.Numerics.Complex* vlP = &vl[0])
                fixed (System.Numerics.Complex* vrP = &vr[0])
                {
                    ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_zgeev(
                        MatrixLayout.ColMajor,
                        jobvl, jobvr,
                        n, dirtyAP, lda,
                        wP,
                        vlP, ldvl,
                        vrP, ldvr);
                }
            }

            if (ret != 0)
            {
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            // 固有値を格納
            eVals = new System.Numerics.Complex[n];
            for (int i = 0; i < n; i++)
            {
                eVals[i] = w[i];
            }

            // 固有ベクトルを格納
            eVecs = new System.Numerics.Complex[n][];
            for (int i = 0; i < n; i++)
            {
                eVecs[i] = new System.Numerics.Complex[ldvr];
                for (int j = 0; j < ldvr; j++)
                {
                    eVecs[i][j] = vr[i * ldvr + j];
                }
            }

            return ret;
        }

        public static int zggev(System.Numerics.Complex[] A, int aRow, int aCol,
            System.Numerics.Complex[] B, int bRow, int bCol,
            out System.Numerics.Complex[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            System.Numerics.Complex[] dirtyA = new System.Numerics.Complex[A.Length];
            A.CopyTo(dirtyA, 0);
            System.Numerics.Complex[] dirtyB = new System.Numerics.Complex[B.Length];
            B.CopyTo(dirtyB, 0);
            int ret = zggev_dirty(
                dirtyA, aRow, aCol, dirtyB, bRow, bCol,
                out eVals, out eVecs);
            return ret;
        }

        public static int zggev_dirty(System.Numerics.Complex[] dirtyA, int aRow, int aCol,
            System.Numerics.Complex[] dirtyB, int bRow, int bCol,
            out System.Numerics.Complex[] eVals, out System.Numerics.Complex[][] eVecs)
        {
            byte jobvl = Job.DontCompute;
            byte jobvr = Job.Compute;
            int n = aCol;
            int lda = n;
            int ldb = n;
            System.Numerics.Complex[] alpha = new System.Numerics.Complex[n];
            System.Numerics.Complex[] beta = new System.Numerics.Complex[n];
            int ldvl = 1;
            System.Numerics.Complex[] vl = new System.Numerics.Complex[ldvl * n];
            int ldvr = n;
            System.Numerics.Complex[] vr = new System.Numerics.Complex[ldvr * n];

            int ret;
            unsafe
            {
                fixed (System.Numerics.Complex* dirtyAP = &dirtyA[0])
                fixed (System.Numerics.Complex* dirtyBP = &dirtyB[0])
                fixed (System.Numerics.Complex* alphaP = &alpha[0])
                fixed (System.Numerics.Complex* betaP = &beta[0])
                fixed (System.Numerics.Complex* vlP = &vl[0])
                fixed (System.Numerics.Complex* vrP = &vr[0])
                {
                    ret =  ImportedFunctions.LAPACKE_zggev(
                        MatrixLayout.ColMajor,
                        jobvl, jobvr,
                        n, dirtyAP, lda,
                        dirtyBP, ldb,
                        alphaP, betaP,
                        vlP, ldvl,
                        vrP, ldvr);
                }
            }
            if (ret != 0)
            {
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }
            eVals = new System.Numerics.Complex[n];
            for (int i = 0; i < n; i++)
            {
                if (beta[i].Magnitude < 1.0e-50)
                {
                    eVals[i] = double.NaN;
                }
                else
                {
                    eVals[i] = alpha[i] / beta[i];
                }
            }

            eVecs = new System.Numerics.Complex[n][];
            for (int i = 0; i < n; i++)
            {
                eVecs[i] = new System.Numerics.Complex[ldvr];
                for (int j = 0; j < ldvr; j++)
                {
                    eVecs[i][j] = vr[i * ldvr + j];
                }
            }
            return ret;
        }

        // conjugate
        public static System.Numerics.Complex[] zlacgv(System.Numerics.Complex[] X)
        {
            int n = X.Length;
            int incX = 1;

            System.Numerics.Complex[] Y = new System.Numerics.Complex[n];
            X.CopyTo(Y, 0);

            unsafe
            {
                fixed (System.Numerics.Complex* YP = &Y[0])
                {
                    // LAPACK
                    //IvyFEM.Lapack.ImportedFunctions.zlacgv_(&n, YP, &incX);

                    // LAPACKE
                    IvyFEM.Lapack.ImportedFunctions.LAPACKE_zlacgv(n, YP, incX);
                }
            }
            return Y;
        }

        public static int zpbsv(out System.Numerics.Complex[] X, out int xRow, out int xCol,
            System.Numerics.Complex[] A, int aRow, int aCol, int superdia,
            System.Numerics.Complex[] B, int bRow, int bCol)
        {
            //if UPLO = 'U', AB(KD + 1 + i - j, j) = A(i, j) for max(1, j - KD) <= i <= j;
            byte uplo = UpperLower.Upper;

            int n = aRow;

            int kd = superdia;
            //The number of superdiagonals of the matrix A if UPLO = 'U'

            int nrhs = bCol;
            // The number of right hand sides, i.e., the number of columns

            int lda = kd + 1;
            // The leading dimension of the array AB.LDAB >= KD + 1.

            int ldb = bRow;

            System.Numerics.Complex[] C = new System.Numerics.Complex[B.Length];
            B.CopyTo(C, 0);

            System.Numerics.Complex[] dirtyA = new System.Numerics.Complex[A.Length];
            A.CopyTo(dirtyA, 0);

            int ret = 0;
            unsafe
            {
                fixed (System.Numerics.Complex* dirtyAP = &dirtyA[0])
                fixed (System.Numerics.Complex* CP = &C[0])
                {
                    ret = IvyFEM.Lapack.ImportedFunctions.LAPACKE_zpbsv(
                        MatrixLayout.ColMajor,
                        uplo, n, kd, nrhs, dirtyAP, lda,
                        CP, ldb);
                }

            }
            if (ret != 0)
            {
                X = null;
                xRow = 0;
                xCol = 0;
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            X = C;
            xRow = bRow;
            xCol = bCol;
            return ret;
        }

        public static int zgbtrs(out System.Numerics.Complex[] X, out int xRow, out int xCol,
            int[] ipiv,
            System.Numerics.Complex[] A, int aRow, int aCol, int subdia, int superdia,
            System.Numerics.Complex[] B, int bRow, int bCol)
        {
            byte trans = Trans.Nop;
            int n = aRow;
            int kl = subdia;
            int ku = superdia;
            int nrhs = bCol;
            // The number of right hand sides, i.e., the number of columns
            // of the matrix B
            int ldab = 2 * kl + ku + 1;
            int ldb = n;
            System.Numerics.Complex[] C = new System.Numerics.Complex[B.Length];
            B.CopyTo(C, 0);

            int ret = -1;
            unsafe
            {
                fixed (System.Numerics.Complex* AP = &A[0])
                fixed (System.Numerics.Complex* CP = &C[0])
                {
                    ret = ImportedFunctions.LAPACKE_zgbtrs(
                        MatrixLayout.ColMajor,
                        trans,
                        n, kl, ku, nrhs,
                        AP, ldab, ipiv, CP, ldb);
                }
            }
            if (ret != 0)
            {
                X = null;
                xRow = 0;
                xCol = 0;
                throw new InvalidOperationException("Error occurred: ret = " + ret);
                //return ret;
            }

            X = C;
            xRow = bRow;
            xCol = bCol;
            return ret;
        }
    }
}
