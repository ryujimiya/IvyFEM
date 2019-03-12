using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Hyperelastic2DDerivedBaseFEM : Elastic2DBaseFEM
    {
        // Ogden
        protected void SolvePrincipalValues(
            double[,] c, out System.Numerics.Complex[] fLambdas, out System.Numerics.Complex[][] cNormals)
        {
            System.Diagnostics.Debug.Assert(c.GetLength(0) == c.GetLength(1));
            System.Diagnostics.Debug.Assert(c.GetLength(0) == 2);
            int dim3 = 3;
            int dim2 = c.GetLength(0);
            IvyFEM.Lapack.DoubleMatrix tmpC = new IvyFEM.Lapack.DoubleMatrix(dim2, dim2);
            for (int row = 0; row < dim2; row++)
            {
                for (int col = 0; col < dim2; col++)
                {
                    tmpC[row, col] = c[row, col];
                }
            }

            System.Numerics.Complex[] eVals;
            System.Numerics.Complex[][] eVecs;
            {
                int ret = IvyFEM.Lapack.Functions.dgeev(
                    tmpC.Buffer, dim2, dim2, out eVals, out eVecs);
                System.Diagnostics.Debug.Assert(ret == 0);
            }

            fLambdas = new System.Numerics.Complex[dim3];
            cNormals = new System.Numerics.Complex[dim3][];
            for (int i = 0; i < dim2; i++)
            {
                System.Numerics.Complex value = System.Numerics.Complex.Sqrt(eVals[i]);
                fLambdas[i] = value;
            }
            fLambdas[dim3 - 1] = 1.0;

            for (int i = 0; i < dim2; i++)
            {
                System.Numerics.Complex[] tmpVec = new System.Numerics.Complex[dim3];
                System.Numerics.Complex maxValue = 0;
                double maxAbs = 0;
                for (int j = 0; j < dim2; j++)
                {
                    System.Numerics.Complex value = eVecs[i][j];
                    tmpVec[j] = value;

                    double abs = value.Magnitude;
                    if (maxAbs > abs)
                    {
                        maxAbs = abs;
                        maxValue = value;
                    }
                }
                tmpVec[dim3 - 1] = 0;
                // 位相調整
                if (maxAbs >= IvyFEM.Constants.PrecisionLowerLimit)
                {
                    tmpVec = IvyFEM.Lapack.Functions.zscal(tmpVec, maxAbs / maxValue);
                }
                // 規格化
                cNormals[i] = IvyFEM.Lapack.Utils.NormalizeComplexVector(tmpVec);
            }
            cNormals[dim3 - 1] = new System.Numerics.Complex[] { 0, 0, 1 };

            {
                // 2軸は直交している
                //System.Numerics.Complex dot = IvyFEM.Lapack.Functions.zdotc(cNormals[1], cNormals[0]);
                //System.Diagnostics.Debug.Assert(dot.Magnitude < IvyFEM.Constants.PrecisionLowerLimit);
                // 右手系の固有ベクトルにする
                System.Numerics.Complex[] newVec = {
                    -System.Numerics.Complex.Conjugate(cNormals[0][1]),
                    System.Numerics.Complex.Conjugate(cNormals[0][0]), 0 };
                bool isDiff = false;
                for (int i = 0; i < dim3; i++)
                {
                    var diff = cNormals[1][i] - newVec[i];
                    if (diff.Magnitude >= IvyFEM.Constants.PrecisionLowerLimit)
                    {
                        isDiff = true;
                        break;
                    }
                }
                if (isDiff)
                {
                    //System.Diagnostics.Debug.WriteLine("change principal vec(y): (" +
                    //    cNormals[1][0] + ", " + cNormals[1][1] + ", " + cNormals[1][2] + ") --> " +
                    //    newVec[0] + ", " + newVec[1] + ", " + newVec[2] + ")");
                    cNormals[1] = newVec;
                }
            }
        }
    }
}
