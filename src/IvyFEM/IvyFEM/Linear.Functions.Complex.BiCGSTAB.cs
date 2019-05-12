using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public partial class Functions
    {
        public static bool ComplexSolvePreconditionedBiCGSTAB(
            out System.Numerics.Complex[] X,
            ComplexSparseMatrix A, System.Numerics.Complex[] B, ComplexSparseMatrix LU,
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
            System.Numerics.Complex[] Mp = new System.Numerics.Complex[n];
            System.Numerics.Complex[] Ms = new System.Numerics.Complex[n];
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

            for (iter = 0; iter < maxIter; iter++)
            {
                ComplexSolveLU(out Mp, LU, p);
                System.Numerics.Complex r0r = IvyFEM.Lapack.Functions.zdotc(r0, r);

                System.Numerics.Complex[] AMp = A * Mp;
                System.Numerics.Complex alpha;
                {
                    System.Numerics.Complex denominator = IvyFEM.Lapack.Functions.zdotc(r0, AMp);
                    alpha = r0r / denominator;
                }
                s = IvyFEM.Lapack.Functions.zaxpy(-alpha, AMp, r);

                ComplexSolveLU(out Ms, LU, s);

                System.Numerics.Complex[] AMs = A * Ms;
                System.Numerics.Complex omega;
                {
                    System.Numerics.Complex denominator = IvyFEM.Lapack.Functions.zdotc(AMs, AMs);
                    System.Numerics.Complex numerator = IvyFEM.Lapack.Functions.zdotc(s, AMs);
                    omega = numerator / denominator;
                }

                x = IvyFEM.Lapack.Functions.zaxpy(alpha, Mp, x);
                x = IvyFEM.Lapack.Functions.zaxpy(omega, Ms, x);
                r = IvyFEM.Lapack.Functions.zaxpy(-omega, AMs, s);

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
                p = IvyFEM.Lapack.Functions.zaxpy(-beta * omega, AMp, p);
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

        public static bool ComplexSolveBiCGSTAB(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B, int fillinLevel,
            double convRatioTolerance)
        {
            int t;
            t = System.Environment.TickCount;
            ComplexSparseMatrix LU = IvyFEM.Linear.Functions.ComplexCalcILU(A, fillinLevel);
            System.Diagnostics.Debug.WriteLine("ComplexSolveBiCGSTAB 1: t= " + (System.Environment.TickCount - t));
            t = System.Environment.TickCount;
            bool success = IvyFEM.Linear.Functions.ComplexSolvePreconditionedBiCGSTAB(out X, A, B, LU,
                convRatioTolerance);
            System.Diagnostics.Debug.WriteLine("ComplexSolveBiCGSTAB 2: t= " + (System.Environment.TickCount - t));
            return success;
        }

        public static bool ComplexSolveBiCGSTABWithPivoting(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B, int fillinLevel,
            double convRatioTolerance)
        {
            int t;
            t = System.Environment.TickCount;
            ComplexSparseMatrix LU;
            int[] pivot;
            ComplexCalcILUWithPivoting(out LU, out pivot, A, fillinLevel);
            ComplexSparseMatrix pivotingA = new ComplexSparseMatrix(A.RowLength, A.ColumnLength);
            for (int i = 0; i < pivotingA.RowLength; i++)
            {
                pivotingA.RowColIndexValues[i] = A.RowColIndexValues[pivot[i]];
            }
            System.Numerics.Complex[] pivotingB = new System.Numerics.Complex[B.Length];
            for (int i = 0; i < B.Length; i++)
            {
                pivotingB[i] = B[pivot[i]];
            }
            A = null;
            B = null;
            System.Diagnostics.Debug.WriteLine("ComplexSolveBiCGStab 1: t= " + (System.Environment.TickCount - t));
            t = System.Environment.TickCount;
            bool success = IvyFEM.Linear.Functions.ComplexSolvePreconditionedBiCGSTAB(out X, pivotingA, pivotingB, LU,
                convRatioTolerance);
            System.Diagnostics.Debug.WriteLine("ComplexSolveBiCGSTAB 2: t= " + (System.Environment.TickCount - t));
            return success;
        }
    }
}