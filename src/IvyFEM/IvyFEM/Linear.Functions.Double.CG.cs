using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public partial class Functions
    {
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
    }
}