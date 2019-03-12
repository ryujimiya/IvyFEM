using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Lapack
{
    public class Utils
    {
        public static System.Numerics.Complex[] Conjugate(System.Numerics.Complex[] A)
        {
            return IvyFEM.Lapack.Functions.zlacgv(A);
        }

        public static double[] NormalizeDoubleVector(double[] X)
        {
            double norm = IvyFEM.Lapack.Functions.dnrm2(X);
            if (Math.Abs(norm) < IvyFEM.Constants.PrecisionLowerLimit)
            {
                // 0ベクトル
                return X;
            }
            double[] normalized = IvyFEM.Lapack.Functions.dscal(X, 1.0 / norm);
            return normalized;
        }

        public static System.Numerics.Complex[] NormalizeComplexVector(System.Numerics.Complex[] X)
        {
            System.Numerics.Complex squareNorm = IvyFEM.Lapack.Functions.zdotc(X, X);
            System.Numerics.Complex norm = System.Numerics.Complex.Sqrt(squareNorm);
            if (norm.Magnitude < IvyFEM.Constants.PrecisionLowerLimit)
            {
                // 0ベクトル
                return X;
            }
            System.Numerics.Complex[] normalized = IvyFEM.Lapack.Functions.zscal(X, 1.0 / norm);
            return normalized;
        }
    }
}
