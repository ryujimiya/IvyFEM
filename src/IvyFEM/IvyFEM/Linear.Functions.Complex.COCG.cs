using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public partial class Functions
    {
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

    }
}