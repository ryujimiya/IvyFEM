using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public class Utils
    {
        public static bool DoubleMatrixIsSymmetric(double[,] A)
        {
            bool isSymmetric = true;
            System.Diagnostics.Debug.Assert(A.GetLength(0) == A.GetLength(1));
            int n = A.GetLength(0);
            for (int row = 0; row < n; row++)
            {
                for (int col = 0; col < row; col++)
                {
                    double diff = A[row, col] - A[col, row];
                    if (Math.Abs(diff) >= IvyFEM.Constants.PrecisionLowerLimit)
                    {
                        isSymmetric = false;
                        break;
                    }
                }
                if (!isSymmetric)
                {
                    break;
                }
            }
            return isSymmetric;
        }

        private static void GetBandMatrixSubDiaSuperDia(
            bool[,] matPattern,
            out int rowcolLength,
            out int subdiaLength,
            out int superdiaLength)
        {
            rowcolLength = matPattern.GetLength(0);

            // subdiaサイズ、superdiaサイズを取得する
            subdiaLength = 0;
            superdiaLength = 0;
            // Note: c == rowcolLength - 1は除く
            for (int c = 0; c < rowcolLength - 1; c++)
            {
                if (subdiaLength >= (rowcolLength - 1 - c))
                {
                    break;
                }
                int cnt = 0;
                for (int r = rowcolLength - 1; r >= c + 1; r--)
                {
                    // 非０要素が見つかったら抜ける
                    if (matPattern[r, c])
                    {
                        cnt = r - c;
                        break;
                    }
                }
                if (cnt > subdiaLength)
                {
                    subdiaLength = cnt;
                }
            }
            // Note: c == 0は除く
            for (int c = rowcolLength - 1; c >= 1; c--)
            {
                if (superdiaLength >= c)
                {
                    break;
                }
                int cnt = 0;
                for (int r = 0; r <= c - 1; r++)
                {
                    // 非０要素が見つかったら抜ける
                    if (matPattern[r, c])
                    {
                        cnt = c - r;
                        break;
                    }
                }
                if (cnt > superdiaLength)
                {
                    superdiaLength = cnt;
                }
            }
            //System.Diagnostics.Debug.WriteLine("rowcolLength: {0} subdiaLength: {1} superdiaLength: {2}", rowcolLength, subdiaLength, superdiaLength);
        }

        private static bool[,] GetDoubleMatrixNonzeroPattern(DoubleSparseMatrix A)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;
            bool[,] nonzeroPattern = new bool[n, n];
            for (int row = 0; row < n; row++)
            {
                foreach (var pair in A.RowColIndexValues[row])
                {
                    int col = pair.Key;
                    nonzeroPattern[row, col] = true;
                }
            }
            return nonzeroPattern;
        }

        public static bool OrderToDoubleBandMatrix(
            out DoubleSparseMatrix orderedA, out double[] orderedB, out int[] indexs,
            DoubleSparseMatrix A, double[] B)
        {
            orderedA = null;
            orderedB = null;
            indexs = null;

            // バンド幅を縮小する
            // 非０要素のパターンを取得
            bool[,] matPattern = GetDoubleMatrixNonzeroPattern(A);
            int n = matPattern.GetLength(0);
            // subdiagonal、superdiagonalのサイズを取得する
            int iniSubdiaLength = 0;
            int iniSuperdiaLength = 0;
            {
                int rowcolLength;
                int subdiaLength;
                int superdiaLength;
                GetBandMatrixSubDiaSuperDia(matPattern, out rowcolLength, out subdiaLength, out superdiaLength);
                System.Diagnostics.Debug.WriteLine(
                    "Initial rowcolLength: {0} subdiaLength: {1} superdiaLength: {2}",
                    rowcolLength, subdiaLength, superdiaLength);

                iniSubdiaLength = subdiaLength;
                iniSuperdiaLength = superdiaLength;
            }

            // 非０要素出現順に節点番号を格納
            IList<int> newIndexs = new List<int>();
            Queue<int> check = new Queue<int>();
            int[] remain = new int[matPattern.GetLength(0)];
            for (int i = 0; i < matPattern.GetLength(0); i++)
            {
                remain[i] = i;
            }
            while (newIndexs.Count < n)
            {
                for (int iRemain = 0; iRemain < remain.Length; iRemain++)
                {
                    int i = remain[iRemain];
                    if (i == -1) continue;
                    check.Enqueue(i);
                    remain[iRemain] = -1;
                    break;
                }
                while (check.Count > 0)
                {
                    int i = check.Dequeue();
                    newIndexs.Add(i);
                    for (int iRemain = 0; iRemain < remain.Length; iRemain++)
                    {
                        int j = remain[iRemain];
                        if (j == -1) continue;
                        if (matPattern[i, j])
                        {
                            check.Enqueue(j);
                            remain[iRemain] = -1;
                        }
                    }
                }
            }
            System.Diagnostics.Debug.Assert(newIndexs.Count == n);
            DoubleSparseMatrix newA = new DoubleSparseMatrix(n, n);
            // 遅い
            //for (int row = 0; row < n; row++)
            //{
            //    for (int col = 0; col < n; col++)
            //    {
            //        newA[row, col] = A[newIndexs[row], newIndexs[col]];
            //    }
            //}
            int[] oldIndexs = new int[n];
            for (int i = 0; i < n; i++)
            {
                oldIndexs[newIndexs[i]] = i;
            }

            for (int row = 0; row < n; row++)
            {
                foreach (var pair in A.RowColIndexValues[row])
                {
                    int col = pair.Key;
                    double value = pair.Value;
                    newA[oldIndexs[row], oldIndexs[col]] = value;
                }
            }

            // 改善できないこともあるのでチェックする
            bool improved = false;
            // 非０パターンを取得
            bool[,] newMatPattern = GetDoubleMatrixNonzeroPattern(newA);
            // check
            {
                int rowcolLength;
                int subdiaLength;
                int superdiaLength;
                GetBandMatrixSubDiaSuperDia(newMatPattern, out rowcolLength, out subdiaLength, out superdiaLength);
                System.Diagnostics.Debug.WriteLine(
                    "Ordered rowcolLength: {0} subdiaLength: {1} superdiaLength: {2}",
                    rowcolLength, subdiaLength, superdiaLength);

                if (subdiaLength <= iniSubdiaLength && superdiaLength <= iniSuperdiaLength)
                {
                    improved = true;
                }
            }
            if (improved)
            {
                // 置き換え
                orderedA = newA;
                indexs = newIndexs.ToArray();

                orderedB = new double[n];
                for (int row = 0; row < n; row++)
                {
                    orderedB[row] = B[newIndexs[row]];
                }
            }
            else
            {
                System.Diagnostics.Debug.WriteLine("band with not optimized!");
            }
            return improved;
        }

        public static void GetDoubleSparseMatrixSubDiaSuperDia(
            IvyFEM.Linear.DoubleSparseMatrix A,
            out int rowcolLength,
            out int subdiaLength,
            out int superdiaLength)
        {
            // 非０要素のパターンを取得
            bool[,] matPattern = GetDoubleMatrixNonzeroPattern(A);
            int n = matPattern.GetLength(0);
            // subdiagonal、superdiagonalのサイズを取得する
            GetBandMatrixSubDiaSuperDia(matPattern, out rowcolLength, out subdiaLength, out superdiaLength);
        }

        private static bool[,] GetComplexMatrixNonzeroPattern(ComplexSparseMatrix A)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;
            bool[,] nonzeroPattern = new bool[n, n];
            for (int row = 0; row < n; row++)
            {
                foreach (var pair in A.RowColIndexValues[row])
                {
                    int col = pair.Key;
                    nonzeroPattern[row, col] = true;
                }
            }
            return nonzeroPattern;
        }

        public static bool OrderToComplexBandMatrix(
            out ComplexSparseMatrix orderedA, out System.Numerics.Complex[] orderedB, out int[] indexs,
            ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            orderedA = null;
            orderedB = null;
            indexs = null;

            // バンド幅を縮小する
            // 非０要素のパターンを取得
            bool[,] matPattern = GetComplexMatrixNonzeroPattern(A);
            int n = matPattern.GetLength(0);
            // subdiagonal、superdiagonalのサイズを取得する
            int iniSubdiaLength = 0;
            int iniSuperdiaLength = 0;
            {
                int rowcolLength;
                int subdiaLength;
                int superdiaLength;
                GetBandMatrixSubDiaSuperDia(matPattern, out rowcolLength, out subdiaLength, out superdiaLength);
                System.Diagnostics.Debug.WriteLine(
                    "Initial rowcolLength: {0} subdiaLength: {1} superdiaLength: {2}",
                    rowcolLength, subdiaLength, superdiaLength);

                iniSubdiaLength = subdiaLength;
                iniSuperdiaLength = superdiaLength;
            }

            // 非０要素出現順に節点番号を格納
            IList<int> newIndexs = new List<int>();
            Queue<int> check = new Queue<int>();
            int[] remain = new int[matPattern.GetLength(0)];
            for (int i = 0; i < matPattern.GetLength(0); i++)
            {
                remain[i] = i;
            }
            while (newIndexs.Count < n)
            {
                for (int iRemain = 0; iRemain < remain.Length; iRemain++)
                {
                    int i = remain[iRemain];
                    if (i == -1) continue;
                    check.Enqueue(i);
                    remain[iRemain] = -1;
                    break;
                }
                while (check.Count > 0)
                {
                    int i = check.Dequeue();
                    newIndexs.Add(i);
                    for (int iRemain = 0; iRemain < remain.Length; iRemain++)
                    {
                        int j = remain[iRemain];
                        if (j == -1) continue;
                        if (matPattern[i, j])
                        {
                            check.Enqueue(j);
                            remain[iRemain] = -1;
                        }
                    }
                }
            }
            System.Diagnostics.Debug.Assert(newIndexs.Count == n);
            ComplexSparseMatrix newA = new ComplexSparseMatrix(n, n);
            // 遅い
            //for (int row = 0; row < n; row++)
            //{
            //    for (int col = 0; col < n; col++)
            //    {
            //        newA[row, col] = A[newIndexs[row], newIndexs[col]];
            //    }
            //}
            int[] oldIndexs = new int[n];
            for (int i = 0; i < n; i++)
            {
                oldIndexs[newIndexs[i]] = i;
            }

            for (int row = 0; row < n; row++)
            {
                foreach (var pair in A.RowColIndexValues[row])
                {
                    int col = pair.Key;
                    System.Numerics.Complex value = pair.Value;
                    newA[oldIndexs[row], oldIndexs[col]] = value;
                }
            }

            // 改善できないこともあるのでチェックする
            bool improved = false;
            // 非０パターンを取得
            bool[,] newMatPattern = GetComplexMatrixNonzeroPattern(newA);
            // check
            {
                int rowcolLength;
                int subdiaLength;
                int superdiaLength;
                GetBandMatrixSubDiaSuperDia(newMatPattern, out rowcolLength, out subdiaLength, out superdiaLength);
                System.Diagnostics.Debug.WriteLine(
                    "Ordered rowcolLength: {0} subdiaLength: {1} superdiaLength: {2}",
                    rowcolLength, subdiaLength, superdiaLength);

                if (subdiaLength <= iniSubdiaLength && superdiaLength <= iniSuperdiaLength)
                {
                    improved = true;
                }
            }
            if (improved)
            {
                // 置き換え
                orderedA = newA;
                indexs = newIndexs.ToArray();

                orderedB = new System.Numerics.Complex[n];
                for (int row = 0; row < n; row++)
                {
                    orderedB[row] = B[newIndexs[row]];
                }
            }
            else
            {
                System.Diagnostics.Debug.WriteLine("band with not optimized!");
            }
            return improved;
        }

        public static void GetComplexSparseMatrixSubDiaSuperDia(
            IvyFEM.Linear.ComplexSparseMatrix A,
            out int rowcolLength,
            out int subdiaLength,
            out int superdiaLength)
        {
            // 非０要素のパターンを取得
            bool[,] matPattern = GetComplexMatrixNonzeroPattern(A);
            int n = matPattern.GetLength(0);
            // subdiagonal、superdiagonalのサイズを取得する
            GetBandMatrixSubDiaSuperDia(matPattern, out rowcolLength, out subdiaLength, out superdiaLength);
        }
    }
}
