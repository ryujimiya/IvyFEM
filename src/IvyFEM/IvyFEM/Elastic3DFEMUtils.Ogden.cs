using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic3DFEMUtils
    {
        // Ogden
        public static void SolvePrincipalValues(
            double[,] c, out System.Numerics.Complex[] fLambdas, out System.Numerics.Complex[][] cNormals)
        {
            System.Diagnostics.Debug.Assert(c.GetLength(0) == c.GetLength(1));
            System.Diagnostics.Debug.Assert(c.GetLength(0) == 3);
            int dim3 = 3;
            IvyFEM.Lapack.DoubleMatrix tmpC = new IvyFEM.Lapack.DoubleMatrix(dim3, dim3);
            for (int row = 0; row < dim3; row++)
            {
                for (int col = 0; col < dim3; col++)
                {
                    tmpC[row, col] = c[row, col];
                }
            }

            System.Numerics.Complex[] eVals;
            System.Numerics.Complex[][] eVecs;
            {
                int ret = IvyFEM.Lapack.Functions.dgeev(
                    tmpC.Buffer, dim3, dim3, out eVals, out eVecs);
                System.Diagnostics.Debug.Assert(ret == 0);
            }

            fLambdas = new System.Numerics.Complex[dim3];
            cNormals = new System.Numerics.Complex[dim3][];
            for (int i = 0; i < dim3; i++)
            {
                System.Numerics.Complex value = System.Numerics.Complex.Sqrt(eVals[i]);
                fLambdas[i] = value;
            }

            for (int i = 0; i < dim3; i++)
            {
                System.Numerics.Complex[] tmpVec = new System.Numerics.Complex[dim3];
                System.Numerics.Complex maxValue = 0;
                double maxAbs = 0;
                for (int j = 0; j < dim3; j++)
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
                // 位相調整
                if (maxAbs >= IvyFEM.Constants.PrecisionLowerLimit)
                {
                    tmpVec = IvyFEM.Lapack.Functions.zscal(tmpVec, maxAbs / maxValue);
                }
                // 規格化
                cNormals[i] = IvyFEM.Lapack.Utils.NormalizeComplexVector(tmpVec);
            }

            {
                /*
                // 3軸は直交している
                for (int i = 0; i < dim3; i++)
                {
                    for (int j = i + 1; j < dim3; j++)
                    {
                        System.Numerics.Complex dot;
                        dot = IvyFEM.Lapack.Functions.zdotc(cNormals[i], cNormals[j]);
                        System.Diagnostics.Debug.Assert(dot.Magnitude < IvyFEM.Constants.PrecisionLowerLimit);
                    }
                }
                */
                // 右手系の固有ベクトルにする
                OpenTK.Vector3d[] vecs = new OpenTK.Vector3d[dim3]; 
                for (int i = 0; i < dim3; i++)
                {
                    vecs[i] = new OpenTK.Vector3d(cNormals[i][0].Real, cNormals[i][1].Real, cNormals[i][2].Real);
                }
                OpenTK.Vector3d newVec3 = OpenTK.Vector3d.Cross(vecs[0], vecs[1]);
                if (OpenTK.Vector3d.Dot(newVec3, vecs[2]) < 0.0)
                {
                    cNormals[2] = cNormals[2].Select(a => -1.0 * a).ToArray();
                }
            }
        }
    }
}
