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
        // complex

        public static bool ComplexSolveNoPreconCOCG(
            out System.Numerics.Complex[] X,
            ComplexSparseMatrix A, System.Numerics.Complex[] B,
            double convRatioTolerance)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;
            System.Diagnostics.Debug.Assert(B.Length == n);
            double convRatio = convRatioTolerance;
            double convRatioTol = convRatio;
            const int maxIter = IvyFEM.Linear.Constants.MaxIter;
            System.Numerics.Complex[] r = new System.Numerics.Complex[n];
            System.Numerics.Complex[] x = new System.Numerics.Complex[n];
            System.Numerics.Complex[] p = new System.Numerics.Complex[n];
            int iter = 0;

            B.CopyTo(r, 0);
            double sqInvNorm0;
            {
                double sqNorm0 = IvyFEM.Lapack.Functions.zdotc(r, r).Real;
                if (sqNorm0 < IvyFEM.Constants.PrecisionLowerLimit)
                {
                    convRatio = 0;
                    X = x;
                    System.Diagnostics.Debug.WriteLine("iter = " + iter + " norm: " + convRatio);
                    return true;
                }
                sqInvNorm0 = 1.0 / sqNorm0;
            }

            r.CopyTo(p, 0);
            System.Numerics.Complex rr = IvyFEM.Lapack.Functions.zdotu(r, r);

            for (iter = 0; iter < maxIter; iter++)
            {
                System.Numerics.Complex[] Ap = A * p;
                System.Numerics.Complex alpha;
                {
                    System.Numerics.Complex pAp = IvyFEM.Lapack.Functions.zdotu(p, Ap);
                    alpha = rr / pAp;
                }
                r = IvyFEM.Lapack.Functions.zaxpy(-alpha, Ap, r);
                x = IvyFEM.Lapack.Functions.zaxpy(alpha, p, x);

                {
                    double sqNorm = IvyFEM.Lapack.Functions.zdotc(r, r).Real;
                    if (sqNorm * sqInvNorm0 < convRatioTol * convRatioTol)
                    {
                        convRatio = Math.Sqrt(sqNorm * sqInvNorm0);
                        X = x;
                        System.Diagnostics.Debug.WriteLine("iter = " + iter + " norm: " + convRatio);
                        return true;
                    }
                }

                System.Numerics.Complex rrPrev = rr;
                rr = IvyFEM.Lapack.Functions.zdotu(r, r);
                System.Numerics.Complex beta = rr / rrPrev;

                p = IvyFEM.Lapack.Functions.zaxpy(beta, p, r);
            }

            {
                double sqNormRes = IvyFEM.Lapack.Functions.zdotc(r, r).Real;
                convRatio = Math.Sqrt(sqNormRes * sqInvNorm0);
                System.Diagnostics.Debug.WriteLine("iter = " + iter + " norm: " + convRatio);
                X = x;
            }
            System.Diagnostics.Debug.WriteLine("Not converged");
            return false;
        }

        // LU : Lii = 1 Lij (i=1...n -1 j=0...i-1) Uij (i=0...n-1 j=i...n-1)
        // Note: M = L * U
        public static ComplexSparseMatrix ComplexCalcILU(ComplexSparseMatrix A, int fillinLevel)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;

            ComplexSparseMatrix LU = new ComplexSparseMatrix(A);
            int[,] level = new int[n, n];
            for (int row = 0; row < n; row++)
            {
                for (int col = 0; col < n; col++)
                {
                    level[row, col] = (A[row, col].Magnitude >= IvyFEM.Constants.PrecisionLowerLimit) ?
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
                    System.Numerics.Complex LUik = LU[i, k];
                    foreach (var pair in LU.RowColIndexValues[k])
                    {
                        int j = pair.Key;
                        System.Numerics.Complex LUkj = pair.Value;
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

        public static void ComplexSolveLU(
            out System.Numerics.Complex[] X, ComplexSparseMatrix LU, System.Numerics.Complex[] B)
        {
            System.Diagnostics.Debug.Assert(LU.RowLength == LU.ColumnLength);
            int n = LU.RowLength;

            X = new System.Numerics.Complex[n];
            B.CopyTo(X, 0);

            // Ly = b 
            for (int row = 1; row < n; row++)
            {
                foreach (var pair in LU.RowColIndexValues[row])
                {
                    int col = pair.Key;
                    System.Numerics.Complex value = pair.Value;
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
                    System.Numerics.Complex value = pair.Value;
                    if (col >= row + 1 && col < n)
                    {
                        X[row] -= value * X[col]; // LU[row, col]

                    }
                }
                X[row] /= LU[row, row];
            }
        }

        public static bool ComplexSolvePreconditionedCOCG(
            out System.Numerics.Complex[] X,
            ComplexSparseMatrix A, System.Numerics.Complex[] B, ComplexSparseMatrix LU,
            double convRatioTolerance)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;
            System.Diagnostics.Debug.Assert(B.Length == n);
            System.Diagnostics.Debug.Assert(LU.RowLength == LU.ColumnLength);
            System.Diagnostics.Debug.Assert(LU.RowLength == n);
            double convRatio = convRatioTolerance;
            double convRatioTol = convRatio;
            const int maxIter = IvyFEM.Linear.Constants.MaxIter;
            System.Numerics.Complex[] r = new System.Numerics.Complex[n];
            System.Numerics.Complex[] x = new System.Numerics.Complex[n];
            System.Numerics.Complex[] z = new System.Numerics.Complex[n];
            System.Numerics.Complex[] p = new System.Numerics.Complex[n];
            int iter = 0;

            B.CopyTo(r, 0);
            double sqInvNorm0;
            {
                double sqNorm0 = IvyFEM.Lapack.Functions.zdotc(r, r).Real;
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
            IvyFEM.Linear.Functions.ComplexSolveLU(out z, LU, r);
            // 前処理なし
            //z = r;

            z.CopyTo(p, 0);
            System.Numerics.Complex rz = IvyFEM.Lapack.Functions.zdotu(r, z);

            for (iter = 0; iter < maxIter; iter++)
            {
                System.Numerics.Complex[] Ap = A * p;
                System.Numerics.Complex alpha;
                {
                    System.Numerics.Complex pAp = IvyFEM.Lapack.Functions.zdotu(p, Ap);
                    alpha = rz / pAp;
                }
                r = IvyFEM.Lapack.Functions.zaxpy(-alpha, Ap, r);
                x = IvyFEM.Lapack.Functions.zaxpy(alpha, p, x);

                {
                    double sqNorm = IvyFEM.Lapack.Functions.zdotc(r, r).Real;
                    if (sqNorm * sqInvNorm0 < convRatioTol * convRatioTol)
                    {
                        convRatio = Math.Sqrt(sqNorm * sqInvNorm0);
                        X = x;
                        System.Diagnostics.Debug.WriteLine("iter = " + iter + " norm: " + convRatio);
                        return true;
                    }
                }

                // 前処理あり
                IvyFEM.Linear.Functions.ComplexSolveLU(out z, LU, r);
                // 前処理なし
                //z = r;

                System.Numerics.Complex rzPrev = rz;
                rz = IvyFEM.Lapack.Functions.zdotu(r, z);
                System.Numerics.Complex beta = rz / rzPrev;

                p = IvyFEM.Lapack.Functions.zaxpy(beta, p, z);
            }

            {
                double sqNormRes = IvyFEM.Lapack.Functions.zdotc(r, r).Real;
                convRatio = Math.Sqrt(sqNormRes * sqInvNorm0);
                System.Diagnostics.Debug.WriteLine("iter = " + iter + " norm: " + convRatio);
                X = x;
            }
            System.Diagnostics.Debug.WriteLine("Not converged");
            return false;
        }

        public static bool ComplexSolveCOCG(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B, int fillinLevel,
            double convRatioTolerance)
        {
            int t;
            t = System.Environment.TickCount;
            ComplexSparseMatrix LU = IvyFEM.Linear.Functions.ComplexCalcILU(A, fillinLevel);
            System.Diagnostics.Debug.WriteLine("ComplexSolveCOCG 1: t= " + (System.Environment.TickCount - t));
            t = System.Environment.TickCount;
            bool success = IvyFEM.Linear.Functions.ComplexSolvePreconditionedCOCG(out X, A, B, LU,
                convRatioTolerance);
            System.Diagnostics.Debug.WriteLine("ComplexSolveCOCG 2: t= " + (System.Environment.TickCount - t));
            return success;
        }

        // Incomplete Cholesky Factorization
        public static ComplexSparseMatrix ComplexCalcIC(ComplexSparseMatrix A)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;

            ComplexSparseMatrix LU = new ComplexSparseMatrix(n, n);
            System.Numerics.Complex[] D = new System.Numerics.Complex[n];

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
                    System.Numerics.Complex tmp = A[i, j];
                    foreach (var pair in LU.RowColIndexValues[i])
                    {
                        int k = pair.Key;
                        System.Numerics.Complex LUik = pair.Value;
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
                    System.Numerics.Complex tmp = A[i, i];
                    foreach (var pair in LU.RowColIndexValues[i])
                    {
                        int k = pair.Key;
                        System.Numerics.Complex LUik = pair.Value;
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
                    System.Numerics.Complex LUij = pair.Value;
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

        public static bool ComplexSolveICCOCG(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B,
            double convRatioTolerance)
        {
            int t;
            t = System.Environment.TickCount;
            ComplexSparseMatrix LU = IvyFEM.Linear.Functions.ComplexCalcIC(A);
            System.Diagnostics.Debug.WriteLine("ComplexSolveICCOCG 1: t= " + (System.Environment.TickCount - t));
            t = System.Environment.TickCount;
            bool success = IvyFEM.Linear.Functions.ComplexSolvePreconditionedCOCG(out X, A, B, LU,
                convRatioTolerance);
            System.Diagnostics.Debug.WriteLine("ComplexSolveICCG 2: t= " + (System.Environment.TickCount - t));
            return success;
        }

        public static bool ComplexSolveNoPreconBiCGSTAB(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B,
            double convRatioTolerance)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;
            System.Diagnostics.Debug.Assert(B.Length == n);
            double convRatio = convRatioTolerance;
            double tolerance = convRatio;
            const int maxIter = IvyFEM.Linear.Constants.MaxIter;
            System.Numerics.Complex[] r = new System.Numerics.Complex[n];
            System.Numerics.Complex[] r0 = new System.Numerics.Complex[n];
            System.Numerics.Complex[] x = new System.Numerics.Complex[n];
            System.Numerics.Complex[] p = new System.Numerics.Complex[n];
            System.Numerics.Complex[] s = new System.Numerics.Complex[n];
            int iter = 0;

            B.CopyTo(r, 0);
            double sqInvNorm0;
            {
                double sqNorm0 = IvyFEM.Lapack.Functions.zdotc(r, r).Real;
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
            System.Numerics.Complex r0r = IvyFEM.Lapack.Functions.zdotc(r0, r);

            for (iter = 0; iter < maxIter; iter++)
            {
                System.Numerics.Complex[] Ap = A * p;
                System.Numerics.Complex alpha;
                {
                    System.Numerics.Complex denominator = IvyFEM.Lapack.Functions.zdotc(r0, Ap);
                    alpha = r0r / denominator;
                }
                s = IvyFEM.Lapack.Functions.zaxpy(-alpha, Ap, r);

                System.Numerics.Complex[] As = A * s;
                System.Numerics.Complex omega;
                {
                    System.Numerics.Complex denominator = IvyFEM.Lapack.Functions.zdotc(As, As);
                    System.Numerics.Complex numerator = IvyFEM.Lapack.Functions.zdotc(s, As);
                    omega = numerator / denominator;
                }

                x = IvyFEM.Lapack.Functions.zaxpy(alpha, p, x);
                x = IvyFEM.Lapack.Functions.zaxpy(omega, s, x);
                r = IvyFEM.Lapack.Functions.zaxpy(-omega, As, s);

                {
                    double sqNorm = IvyFEM.Lapack.Functions.zdotc(r, r).Real;
                    if (sqNorm * sqInvNorm0 < tolerance * tolerance)
                    {
                        convRatio = Math.Sqrt(sqNorm * sqInvNorm0);
                        X = x;
                        System.Diagnostics.Debug.WriteLine("iter = " + iter + " norm: " + convRatio);
                        return true;
                    }
                }

                System.Numerics.Complex beta;
                {
                    System.Numerics.Complex r0rPrev = r0r;
                    r0r = IvyFEM.Lapack.Functions.zdotc(r0, r);
                    beta = (r0r * alpha) / (r0rPrev * omega);
                }
                p = IvyFEM.Lapack.Functions.zaxpy(beta, p, r);
                p = IvyFEM.Lapack.Functions.zaxpy(-beta * omega, Ap, p);
            }

            {
                double sqNormRes = IvyFEM.Lapack.Functions.zdotc(r, r).Real;
                convRatio = Math.Sqrt(sqNormRes * sqInvNorm0);
                System.Diagnostics.Debug.WriteLine("iter = " + iter + " norm: " + convRatio);
                X = x;
            }
            System.Diagnostics.Debug.WriteLine("Not converged");
            return false;
        }
    }
}
