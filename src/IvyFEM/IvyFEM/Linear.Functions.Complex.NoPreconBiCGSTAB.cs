using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public partial class Functions
    {
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