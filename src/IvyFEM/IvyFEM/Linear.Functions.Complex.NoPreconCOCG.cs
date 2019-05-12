using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public partial class Functions
    {
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
    }
}