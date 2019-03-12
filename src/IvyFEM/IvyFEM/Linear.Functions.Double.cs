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

        public static bool DoubleSolveNoPreconCG(
            out double[] X, DoubleSparseMatrix A, double[] B,
            double convRatioTolerance)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;
            System.Diagnostics.Debug.Assert(B.Length == n);
            double convRatio = convRatioTolerance;
            double tolerance = convRatio;
            const int maxIter = IvyFEM.Linear.Constants.MaxIter;
            double[] r = new double[n];
            double[] x = new double[n];
            double[] p = new double[n];
            int iter = 0;

            B.CopyTo(r, 0);
            double rr;
            double sqInvNorm0;
            {
                double sqNorm0 = IvyFEM.Lapack.Functions.ddot(r, r);
                if (sqNorm0 < IvyFEM.Constants.PrecisionLowerLimit)
                {
                    convRatio = 0;
                    X = x;
                    System.Diagnostics.Debug.WriteLine("iter = " + iter + " norm: " + convRatio);
                    return true;
                }
                rr = sqNorm0;
                sqInvNorm0 = 1.0 / sqNorm0;
            }

            r.CopyTo(p, 0);
            for (iter = 0; iter < maxIter; iter++)
            {
                double[] Ap = A * p;
                double alpha;
                {
                    double pAp = IvyFEM.Lapack.Functions.ddot(p, Ap);
                    alpha = rr / pAp;
                }
                r = IvyFEM.Lapack.Functions.daxpy(-alpha, Ap, r);
                x = IvyFEM.Lapack.Functions.daxpy(alpha, p, x);

                double rrPrev = rr;
                {
                    double sqNorm = IvyFEM.Lapack.Functions.ddot(r, r);
                    if (sqNorm * sqInvNorm0 < tolerance * tolerance)
                    {
                        convRatio = Math.Sqrt(sqNorm * sqInvNorm0);
                        X = x;
                        System.Diagnostics.Debug.WriteLine("iter = " + iter + " norm: " + convRatio);
                        return true;
                    }
                    rr = sqNorm;
                }
                double beta = rr / rrPrev;

                p = IvyFEM.Lapack.Functions.daxpy(beta, p, r);
            }

            {
                double sqNormRes = IvyFEM.Lapack.Functions.ddot(r, r);
                convRatio = Math.Sqrt(sqNormRes * sqInvNorm0);
                System.Diagnostics.Debug.WriteLine("iter = " + iter + " norm: " + convRatio);
                X = x;
            }
            System.Diagnostics.Debug.WriteLine("Not converged");
            return false;
        }

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

        public static bool DoubleSolvePreconditionedCG(
            out double[] X, DoubleSparseMatrix A, double[] B, DoubleSparseMatrix LU,
            double convRatioTolerance)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;
            System.Diagnostics.Debug.Assert(B.Length == n);
            System.Diagnostics.Debug.Assert(LU.RowLength == LU.ColumnLength);
            System.Diagnostics.Debug.Assert(LU.RowLength == n);
            double convRatio = convRatioTolerance;
            double tolerance = convRatio;
            const int maxIter = IvyFEM.Linear.Constants.MaxIter;
            double[] r = new double[n];
            double[] x = new double[n];
            double[] z = new double[n];
            double[] p = new double[n];
            int iter = 0;

            B.CopyTo(r, 0);
            double sqInvNorm0;
            {
                double sqNorm0 = IvyFEM.Lapack.Functions.ddot(r, r);
                if (sqNorm0 < IvyFEM.Constants.PrecisionLowerLimit)
                {
                    convRatio = 0;
                    X = x;
                    System.Diagnostics.Debug.WriteLine("iter = " + iter + " norm: " + convRatio);
                    return true;
                }
                sqInvNorm0 = 1.0 / sqNorm0;
            }

            // 前処理あり
            IvyFEM.Linear.Functions.DoubleSolveLU(out z, LU, r);
            // 前処理なし
            //z = r;

            z.CopyTo(p, 0);
            double rz = IvyFEM.Lapack.Functions.ddot(r, z);

            for (iter = 0; iter < maxIter; iter++)
            {
                double[] Ap = A * p;
                double alpha;
                {
                    double pAp = IvyFEM.Lapack.Functions.ddot(p, Ap);
                    alpha = rz / pAp;
                }
                r = IvyFEM.Lapack.Functions.daxpy(-alpha, Ap, r);
                x = IvyFEM.Lapack.Functions.daxpy(alpha, p, x);

                {
                    double sqNorm = IvyFEM.Lapack.Functions.ddot(r, r);
                    if (sqNorm * sqInvNorm0 < tolerance * tolerance)
                    {
                        convRatio = Math.Sqrt(sqNorm * sqInvNorm0);
                        X = x;
                        System.Diagnostics.Debug.WriteLine("iter = " + iter + " norm: " + convRatio);
                        return true;
                    }
                }

                // 前処理あり
                IvyFEM.Linear.Functions.DoubleSolveLU(out z, LU, r);
                // 前処理なし
                //z = r;

                double rzPrev = rz;
                rz = IvyFEM.Lapack.Functions.ddot(r, z);
                double beta = rz / rzPrev;

                p = IvyFEM.Lapack.Functions.daxpy(beta, p, z);
            }

            {
                double sqNormRes = IvyFEM.Lapack.Functions.ddot(r, r);
                convRatio = Math.Sqrt(sqNormRes * sqInvNorm0);
                System.Diagnostics.Debug.WriteLine("iter = " + iter + " norm: " + convRatio);
                X = x;
            }
            System.Diagnostics.Debug.WriteLine("Not converged");
            return false;
        }

        public static bool DoubleSolveCG(out double[] X, DoubleSparseMatrix A, double[] B, int fillinLevel,
            double convRatioTolerance)
        {
            int t;
            t = System.Environment.TickCount;
            DoubleSparseMatrix LU = IvyFEM.Linear.Functions.DoubleCalcILU(A, fillinLevel);
            System.Diagnostics.Debug.WriteLine("DoubleSolveCG 1: t= " + (System.Environment.TickCount - t));
            t = System.Environment.TickCount;
            bool success = IvyFEM.Linear.Functions.DoubleSolvePreconditionedCG(out X, A, B, LU,
                convRatioTolerance);
            System.Diagnostics.Debug.WriteLine("DoubleSolveCG 2: t= " + (System.Environment.TickCount - t));
            return success;
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

        public static bool DoubleSolveCGWithPivoting(out double[] X, DoubleSparseMatrix A, double[] B, 
            int fillinLevel,
            double convRatioTolerance)
        {
            int t;
            t = System.Environment.TickCount;
            DoubleSparseMatrix LU;
            int[] pivot;
            DoubleCalcILUWithPivoting(out LU, out pivot, A, fillinLevel);
            DoubleSparseMatrix pivotingA = new DoubleSparseMatrix(A.RowLength, A.ColumnLength);
            for (int i = 0; i < pivotingA.RowLength; i++)
            {
                pivotingA.RowColIndexValues[i] = A.RowColIndexValues[pivot[i]];
            }
            double[] pivotingB = new double[B.Length];
            for (int i = 0; i < B.Length; i++)
            {
                pivotingB[i] = B[pivot[i]];
            }
            A = null;
            B = null;
            System.Diagnostics.Debug.WriteLine("DoubleSolveCGWithPivoting 1: t= " + (System.Environment.TickCount - t));
            t = System.Environment.TickCount;
            bool success = IvyFEM.Linear.Functions.DoubleSolvePreconditionedCG(out X, pivotingA, pivotingB, LU,
                convRatioTolerance);
            System.Diagnostics.Debug.WriteLine("DoubleSolveCGWithPivoting 2: t= " + (System.Environment.TickCount - t));
            return success;
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

        public static bool DoubleSolveICCG(out double[] X, DoubleSparseMatrix A, double[] B,
            double convRatioTolerance)
        {
            int t;
            t = System.Environment.TickCount;
            DoubleSparseMatrix LU = IvyFEM.Linear.Functions.DoubleCalcIC(A);
            System.Diagnostics.Debug.WriteLine("DoubleSolveICCG 1: t= " + (System.Environment.TickCount - t));
            t = System.Environment.TickCount;
            bool success = IvyFEM.Linear.Functions.DoubleSolvePreconditionedCG(out X, A, B, LU,
                convRatioTolerance);
            System.Diagnostics.Debug.WriteLine("DoubleSolveICCG 2: t= " + (System.Environment.TickCount - t));
            return success;
        }

        public static bool DoubleSolveNoPreconBiCGSTAB(
            out double[] X, DoubleSparseMatrix A, double[] B,
            double convRatioTolerance)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;
            System.Diagnostics.Debug.Assert(B.Length == n);
            double convRatio = convRatioTolerance;
            double tolerance = convRatio;
            const int maxIter = IvyFEM.Linear.Constants.MaxIter;
            double[] r = new double[n];
            double[] r0 = new double[n];
            double[] x = new double[n];
            double[] p = new double[n];
            double[] s = new double[n];
            int iter = 0;

            B.CopyTo(r, 0);
            double sqInvNorm0;
            {
                double sqNorm0 = IvyFEM.Lapack.Functions.ddot(r, r);
                if (sqNorm0 < IvyFEM.Constants.PrecisionLowerLimit)
                {
                    convRatio = 0;
                    X = x;
                    System.Diagnostics.Debug.WriteLine("iter = " + iter + " norm: " + convRatio);
                    return true;
                }
                sqInvNorm0 = 1.0 / sqNorm0;
            }

            r.CopyTo(r0, 0);
            r.CopyTo(p, 0);
            double r0r = IvyFEM.Lapack.Functions.ddot(r0, r);

            for (iter = 0; iter < maxIter; iter++)
            {
                double[] Ap = A * p;
                double alpha;
                {
                    double denominator = IvyFEM.Lapack.Functions.ddot(r0, Ap);
                    alpha = r0r / denominator;
                }
                s = IvyFEM.Lapack.Functions.daxpy(-alpha, Ap, r);

                double[] As = A * s;
                double omega;
                {
                    double denominator = IvyFEM.Lapack.Functions.ddot(As, As);
                    double numerator = IvyFEM.Lapack.Functions.ddot(s, As);
                    omega = numerator / denominator;
                }

                x = IvyFEM.Lapack.Functions.daxpy(alpha, p, x);
                x = IvyFEM.Lapack.Functions.daxpy(omega, s, x);
                r = IvyFEM.Lapack.Functions.daxpy(-omega, As, s);

                {
                    double sqNorm = IvyFEM.Lapack.Functions.ddot(r, r);
                    if (sqNorm * sqInvNorm0 < tolerance * tolerance)
                    {
                        convRatio = Math.Sqrt(sqNorm * sqInvNorm0);
                        X = x;
                        System.Diagnostics.Debug.WriteLine("iter = " + iter + " norm: " + convRatio);
                        return true;
                    }
                }

                double beta;
                {
                    double r0rPrev = r0r;
                    r0r = IvyFEM.Lapack.Functions.ddot(r0, r);
                    beta = (r0r * alpha) / (r0rPrev * omega);
                }
                p = IvyFEM.Lapack.Functions.daxpy(beta, p, r);
                p = IvyFEM.Lapack.Functions.daxpy(-beta * omega, Ap, p);
            }

            {
                double sqNormRes = IvyFEM.Lapack.Functions.ddot(r, r);
                convRatio = Math.Sqrt(sqNormRes * sqInvNorm0);
                System.Diagnostics.Debug.WriteLine("iter = " + iter + " norm: " + convRatio);
                X = x;
            }
            System.Diagnostics.Debug.WriteLine("Not converged");
            return false;
        }

    }
}
