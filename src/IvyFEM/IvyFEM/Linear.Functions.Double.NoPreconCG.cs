using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public partial class Functions
    {
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
    }
}