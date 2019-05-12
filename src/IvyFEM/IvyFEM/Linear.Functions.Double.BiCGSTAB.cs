using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public partial class Functions
    {
        public static bool DoubleSolvePreconditionedBiCGSTAB(
            out double[] X, DoubleSparseMatrix A, double[] B, DoubleSparseMatrix LU,
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
            double[] Mp = new double[n];
            double[] Ms = new double[n];
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

            for (iter = 0; iter < maxIter; iter++)
            {
                DoubleSolveLU(out Mp, LU, p);

                double r0r = IvyFEM.Lapack.Functions.ddot(r0, r);
                double[] AMp = A * Mp;
                double alpha;
                {
                    double denominator = IvyFEM.Lapack.Functions.ddot(r0, AMp);
                    alpha = r0r / denominator;
                }
                s = IvyFEM.Lapack.Functions.daxpy(-alpha, AMp, r);

                DoubleSolveLU(out Ms, LU, s);

                double[] AMs = A * Ms;
                double omega;
                {
                    double denominator = IvyFEM.Lapack.Functions.ddot(AMs, AMs);
                    double numerator = IvyFEM.Lapack.Functions.ddot(s, AMs);
                    omega = numerator / denominator;
                }

                x = IvyFEM.Lapack.Functions.daxpy(alpha, Mp, x);
                x = IvyFEM.Lapack.Functions.daxpy(omega, Ms, x);
                r = IvyFEM.Lapack.Functions.daxpy(-omega, AMs, s);

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
                p = IvyFEM.Lapack.Functions.daxpy(-beta * omega, AMp, p);
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

        public static bool DoubleSolveBiCGSTAB(out double[] X, DoubleSparseMatrix A, double[] B, int fillinLevel,
            double convRatioTolerance)
        {
            int t;
            t = System.Environment.TickCount;
            DoubleSparseMatrix LU = IvyFEM.Linear.Functions.DoubleCalcILU(A, fillinLevel);
            System.Diagnostics.Debug.WriteLine("DoubleSolveBiCGStab 1: t= " + (System.Environment.TickCount - t));
            t = System.Environment.TickCount;
            bool success = IvyFEM.Linear.Functions.DoubleSolvePreconditionedBiCGSTAB(out X, A, B, LU,
                convRatioTolerance);
            System.Diagnostics.Debug.WriteLine("DoubleSolveBiCGSTAB 2: t= " + (System.Environment.TickCount - t));
            return success;
        }

        public static bool DoubleSolveBiCGSTABWithPivoting(out double[] X, DoubleSparseMatrix A, double[] B,
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
            System.Diagnostics.Debug.WriteLine("DoubleSolveBiCGStab 1: t= " + (System.Environment.TickCount - t));
            t = System.Environment.TickCount;
            bool success = IvyFEM.Linear.Functions.DoubleSolvePreconditionedBiCGSTAB(out X, pivotingA, pivotingB, LU,
                convRatioTolerance);
            System.Diagnostics.Debug.WriteLine("DoubleSolveBiCGSTAB 2: t= " + (System.Environment.TickCount - t));
            return success;
        }
    }
}