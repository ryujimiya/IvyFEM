using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Lapack
{
    public class Utils
    {
        public static System.Numerics.Complex[] Conjugate(System.Numerics.Complex[] A)
        {
            return IvyFEM.Lapack.Functions.zlacgv(A);
        }

        public static double[] NormalizeDoubleVector(double[] X)
        {
            double norm = IvyFEM.Lapack.Functions.dnrm2(X);
            if (Math.Abs(norm) < IvyFEM.Constants.PrecisionLowerLimit)
            {
                // 0ベクトル
                return X;
            }
            double[] normalized = IvyFEM.Lapack.Functions.dscal(X, 1.0 / norm);
            return normalized;
        }

        public static System.Numerics.Complex[] NormalizeComplexVector(System.Numerics.Complex[] X)
        {
            System.Numerics.Complex squareNorm = IvyFEM.Lapack.Functions.zdotc(X, X);
            System.Numerics.Complex norm = System.Numerics.Complex.Sqrt(squareNorm);
            if (norm.Magnitude < IvyFEM.Constants.PrecisionLowerLimit)
            {
                // 0ベクトル
                return X;
            }
            System.Numerics.Complex[] normalized = IvyFEM.Lapack.Functions.zscal(X, 1.0 / norm);
            return normalized;
        }

        public static DoubleMatrix DoubleInverse1(DoubleMatrix A)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;
            DoubleMatrix workA = new DoubleMatrix(A);
            DoubleMatrix workB = new DoubleMatrix(n, n);
            workB.Identity(); // 単位行列
            double[] a = workA.Buffer;
            double[] b = workB.Buffer;
            // [A][X] = [B]
            //  [B]の内容が書き換えられるので、matXを新たに生成せず、Bを出力に指定している
            int xRow = 0;
            int xCol = 0;
            IvyFEM.Lapack.Functions.dgesv(out b, out xRow, out xCol, a, n, n, b, n, n);

            bool alloc = false; // 指定したバッファを使用する
            DoubleMatrix X = new DoubleMatrix(b, xRow, xCol, alloc);
            return X;
        }

        public static ComplexMatrix ComplexInverse1(ComplexMatrix A)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;
            ComplexMatrix workA = new ComplexMatrix(A);
            ComplexMatrix workB = new ComplexMatrix(n, n);
            workB.Identity(); // 単位行列
            System.Numerics.Complex[] a = workA.Buffer;
            System.Numerics.Complex[] b = workB.Buffer;
            // [A][X] = [B]
            //  [B]の内容が書き換えられるので、matXを新たに生成せず、Bを出力に指定している
            int xRow = 0;
            int xCol = 0;
            IvyFEM.Lapack.Functions.zgesv(out b, out xRow, out xCol, a, n, n, b, n, n);

            bool alloc = false; // 指定したバッファを使用する
            ComplexMatrix X = new ComplexMatrix(b, xRow, xCol, alloc);
            return X;
        }

        public static DoubleMatrix DoubleInverse2(DoubleMatrix A)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;
            double[] LU = null;
            int[] piv = null;
            // LU分解
            IvyFEM.Lapack.Functions.dgetrf(out LU, out piv, A.Buffer, n, n);
            // 逆行列を求める
            double[] _X = null;
            IvyFEM.Lapack.Functions.dgetri(out _X, piv, LU, n, n);

            bool alloc = false; // 指定したバッファを使用する
            DoubleMatrix X = new DoubleMatrix(_X, n, n, alloc);
            return X;
        }

        public static ComplexMatrix ComplexInverse2(ComplexMatrix A)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;
            System.Numerics.Complex[] LU = null;
            int[] piv = null;
            // LU分解
            IvyFEM.Lapack.Functions.zgetrf(out LU, out piv, A.Buffer, n, n);
            // 逆行列を求める
            System.Numerics.Complex[] _X = null;
            IvyFEM.Lapack.Functions.zgetri(out _X, piv, LU, n, n);

            bool alloc = false; // 指定したバッファを使用する
            ComplexMatrix X = new ComplexMatrix(_X, n, n, alloc);
            return X;
        }
    }
}
