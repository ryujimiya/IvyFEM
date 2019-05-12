using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public partial class Functions
    {
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