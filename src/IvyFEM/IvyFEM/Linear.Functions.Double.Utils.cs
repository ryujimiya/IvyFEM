using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public partial class Functions
    {
        /////////////////////////////////////////////////////////////////////////////////
        // double

        // LU : Lii = 1 Lij (i=1...n -1 j=0...i-1) Uij (i=0...n-1 j=i...n-1)
        // Note: M = L * U
        public static DoubleSparseMatrix DoubleCalcILU(DoubleSparseMatrix A, int fillinLevel)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;

            DoubleSparseMatrix LU = new DoubleSparseMatrix(A);
            int[,] level = new int[n, n];
            for (int row = 0; row < n; row++)
            {
                for (int col = 0; col < n; col++)
                {
                    level[row, col] = (Math.Abs(A[row, col]) >= IvyFEM.Constants.PrecisionLowerLimit) ?
                        0 : (fillinLevel + 1);
                }
            }

            for (int i = 1; i < n; i++)
            {
                for (int k = 0; k <= (i - 1); k++)
                {
                    //if (!LU.RowColIndexValues[i].ContainsKey(k))
                    //{
                    //    continue;
                    //}
                    if (level[i, k] > fillinLevel)
                    {
                        continue;
                    }
                    LU[i, k] /= LU[k, k];
                    double LUik = LU[i, k];
                    foreach (var pair in LU.RowColIndexValues[k])
                    {
                        int j = pair.Key;
                        double LUkj = pair.Value;
                        if (j >= k + 1 && j < n)
                        {
                            //
                        }
                        else
                        {
                            continue;
                        }

                        level[i, j] = Math.Min(level[i, j], level[i, k] + level[k, j] + 1);
                        if (level[i, j] <= fillinLevel)
                        {
                            LU[i, j] -= LUik * LUkj;
                        }
                    }
                }
            }
            return LU;
        }

        public static void DoubleSolveLU(out double[] X, DoubleSparseMatrix LU, double[] B)
        {
            System.Diagnostics.Debug.Assert(LU.RowLength == LU.ColumnLength);
            int n = LU.RowLength;

            X = new double[n];
            B.CopyTo(X, 0);

            // Ly = b 
            for (int row = 1; row < n; row++)
            {
                foreach (var pair in LU.RowColIndexValues[row])
                {
                    int col = pair.Key;
                    double value = pair.Value;
                    if (col >= 0 && col < row)
                    {
                        X[row] -= value * X[col]; // LU[row, col]
                    }
                }
            }

            // Ux = y
            for (int row = n - 2; row >= 0; row--)
            {
                foreach (var pair in LU.RowColIndexValues[row])
                {
                    int col = pair.Key;
                    double value = pair.Value;
                    if (col >= row + 1 && col < n)
                    {
                        X[row] -= value * X[col]; // LU[row, col]

                    }
                }
                X[row] /= LU[row, row];
            }
        }

        public static void DoubleCalcILUWithPivoting(out DoubleSparseMatrix LU, out int[] pivot, DoubleSparseMatrix A, int fillinLevel)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;

            LU = new DoubleSparseMatrix(A);
            pivot = new int[n];
            for (int row = 0; row < n; row++)
            {
                pivot[row] = row;
            }

            int[,] level = new int[n, n];
            for (int row = 0; row < n; row++)
            {
                for (int col = 0; col < n; col++)
                {
                    level[row, col] = (Math.Abs(A[row, col]) >= IvyFEM.Constants.PrecisionLowerLimit) ?
                        0 : (fillinLevel + 1);
                }
            }

            for (int k = 0; k < (n - 1); k++)
            {
                int p = k;
                {
                    double max = Math.Abs(LU[k, k]);
                    for (int i = k + 1; i < n; i++)
                    {
                        double abs = Math.Abs(LU[i, k]);
                        if (abs > max)
                        {
                            max = abs;
                            p = i;
                        }
                    }
                }
                if (k != p)
                {
                    {
                        int tmp = pivot[k];
                        pivot[k] = pivot[p];
                        pivot[p] = tmp;
                    }
                    {
                        var tmp = LU.RowColIndexValues[k];
                        LU.RowColIndexValues[k] = LU.RowColIndexValues[p];
                        LU.RowColIndexValues[p] = tmp;
                    }
                    for (int j = 0; j < n; j++)
                    {
                        int tmp = level[k, j];
                        level[k, j] = level[p, j];
                        level[p, j] = tmp;
                    }
                }
                for (int i = k + 1; i < n; i++)
                {
                    if (level[i, k] > fillinLevel)
                    {
                        continue;
                    }
                    System.Diagnostics.Debug.Assert(LU.RowColIndexValues[k].ContainsKey(k));
                    LU[i, k] /= LU[k, k];
                    double LUik = LU[i, k];
                    foreach (var pair in LU.RowColIndexValues[k])
                    {
                        int j = pair.Key;
                        double LUkj = pair.Value;
                        if (j >= k + 1 && j < n)
                        {

                        }
                        else
                        {
                            continue;
                        }
                        level[i, j] = Math.Min(level[i, j], level[i, k] + level[k, j] + 1);
                        if (level[i, j] <= fillinLevel)
                        {
                            LU[i, j] -= LUik * LUkj; // LU[i, k] LU[k, j]
                        }
                    }
                }
            }
            for (int i = 0; i < n; i++)
            {
                System.Diagnostics.Debug.Assert(LU.RowColIndexValues[i].ContainsKey(i));
            }
        }

        // Incomplete Cholesky Factorization
        public static DoubleSparseMatrix DoubleCalcIC(DoubleSparseMatrix A)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;

            DoubleSparseMatrix LU = new DoubleSparseMatrix(n, n);
            double[] D = new double[n];

            D[0] = A[0, 0];
            LU[0, 0] = 1.0;

            // L
            for (int i = 1; i < n; i++)
            {
                for (int j = 0; j <= (i - 1); j++)
                {
                    if (!A.RowColIndexValues[i].ContainsKey(j))
                    {
                        continue;
                    }
                    double tmp = A[i, j];
                    foreach (var pair in LU.RowColIndexValues[i])
                    {
                        int k = pair.Key;
                        double LUik = pair.Value;
                        if (k >= 0 && k <= (j - 1))
                        {
                            //
                        }
                        else
                        {
                            continue;
                        }
                        tmp -= LUik * LU[j, k] * D[k];  // LU[i, k]
                    }
                    LU[i, j] = (1.0 / D[j]) * tmp;
                }

                // i == jのとき
                {
                    double tmp = A[i, i];
                    foreach (var pair in LU.RowColIndexValues[i])
                    {
                        int k = pair.Key;
                        double LUik = pair.Value;
                        if (k >= 0 && k <= (i - 1))
                        {
                            //
                        }
                        else
                        {
                            continue;
                        }
                        tmp -= LUik * LUik * D[k];  // LU[i, k]
                    }
                    D[i] = tmp;
                    LU[i, i] = 1.0;
                }
            }

            // DL^T
            for (int i = 0; i < n; i++)
            {
                foreach (var pair in LU.RowColIndexValues[i])
                {
                    int j = pair.Key;
                    double LUij = pair.Value;
                    if (j >= 0 && j <= (i - 1))
                    {
                        //
                    }
                    else
                    {
                        continue;
                    }
                    LU[j, i] = D[j] * LUij; // LU[i, j]
                }
                LU[i, i] = D[i];
            }
            return LU;
        }
    }
}
