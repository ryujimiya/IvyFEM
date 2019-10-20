using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public class DoubleSparseMatrix
    {
        public int RowLength { get; protected set; } = 0;
        public int ColumnLength { get; protected set; } = 0;
        public Dictionary<int, double>[] RowColIndexValues { get; protected set; } = null;

        public DoubleSparseMatrix()
        {
            Clear();
        }

        public DoubleSparseMatrix(int rowLength, int columnLength)
        {
            Resize(rowLength, columnLength);
        }

        public DoubleSparseMatrix(DoubleSparseMatrix src)
        {
            Copy(src);
        }

        public DoubleSparseMatrix(IvyFEM.Lapack.DoubleMatrix denseM)
        {
            Resize(denseM.RowLength, denseM.ColumnLength);
            for (int row = 0; row < denseM.RowLength; row++)
            {
                for (int col = 0; col < denseM.ColumnLength; col++)
                {
                    if (Math.Abs(denseM[row, col]) >= IvyFEM.Constants.PrecisionLowerLimit)
                    {
                        this[row, col] = denseM[row, col];
                    }
                }
            }
        }

        public void Resize(int rowLength, int columnLength)
        {
            RowLength = rowLength;
            ColumnLength = columnLength;
            RowColIndexValues = new Dictionary<int, double>[RowLength];
            for (int row = 0; row < RowLength; row++)
            {
                RowColIndexValues[row] = new Dictionary<int, double>();
            }
        }

        public void Clear()
        {
            RowColIndexValues = null;
            RowLength = 0;
            ColumnLength = 0;
        }

        public double this[int row, int col]
        {
            get
            {
                if (row < 0 || RowLength <= row || col < 0 || ColumnLength <= col)
                {
                    throw new IndexOutOfRangeException();
                }
                if (RowColIndexValues == null)
                {
                    throw new InvalidOperationException();
                }
                double value = 0;
                if (RowColIndexValues[row].ContainsKey(col))
                {
                    value = RowColIndexValues[row][col];
                }
                return value;
            }
            set
            {
                if (row < 0 || RowLength <= row || col < 0 || ColumnLength <= col)
                {
                    throw new IndexOutOfRangeException();
                }
                if (RowColIndexValues == null)
                {
                    throw new InvalidOperationException();
                }
                if (RowColIndexValues[row].ContainsKey(col))
                {
                    if (Math.Abs(value) >= IvyFEM.Constants.PrecisionLowerLimit)
                    {
                        RowColIndexValues[row][col] = value;
                    }
                    else
                    {
                        RowColIndexValues[row].Remove(col);
                    }
                }
                else
                {
                    if (Math.Abs(value) >= IvyFEM.Constants.PrecisionLowerLimit)
                    {
                        RowColIndexValues[row].Add(col, value);
                    }
                }
            }
        }

        public void Copy(DoubleSparseMatrix src)
        {
            Resize(src.RowLength, src.ColumnLength);
            for (int row = 0; row < src.RowLength; row++)
            {
                var srcColIndexValues = src.RowColIndexValues[row];
                foreach (var pair in srcColIndexValues)
                {
                    int col = pair.Key;
                    double value = pair.Value;
                    RowColIndexValues[row].Add(col, value);
                }
            }
        }

        public void Zero()
        {
            for (int row = 0; row < RowLength; row++)
            {
                RowColIndexValues[row].Clear();
            }
        }

        public void Identity()
        {
            Zero();
            for (int i = 0; i < RowLength; i++)
            {
                this[i, i] = 1;
            }
        }

        /*
        public void Transpose()
        {
            throw new NotImplementedException();
        }
        */

        public static DoubleSparseMatrix operator *(DoubleSparseMatrix A, DoubleSparseMatrix B)
        {
            int aRow = A.RowLength;
            int aCol = A.ColumnLength;
            int bRow = B.RowLength;
            int bCol = B.ColumnLength;
            if (aCol != bRow)
            {
                throw new ArgumentException("Mismatched size: aCol != bRow(" + aCol + " != " + bRow + ")");
            }

            int cRow = aRow;
            int cCol = bCol;
            DoubleSparseMatrix C = new DoubleSparseMatrix(cRow, cCol);

            for (int i = 0; i < aRow; i++)
            {
                var colIndexValues1 = A.RowColIndexValues[i];
                for (int j = 0; j < bCol; j++)
                {
                    foreach (var pair1 in colIndexValues1)
                    {
                        int k = pair1.Key;
                        double value1 = pair1.Value; // Aik

                        if (B.RowColIndexValues[k].ContainsKey(j))
                        {
                            double value2 = B.RowColIndexValues[k][j]; // Bkj
                            C[i, j] += value1 * value2;
                        }
                    }

                }
            }
            return C;
        }

        public static double[] operator *(DoubleSparseMatrix A, double[] b)
        {
            double[] c = new double[A.RowLength];
            for (int row = 0; row < A.RowLength; row++)
            {
                var colIndexValues = A.RowColIndexValues[row];
                foreach (var pair in colIndexValues)
                {
                    int col = pair.Key;
                    double value = pair.Value;
                    c[row] += value * b[col];
                }
            }
            return c;
        }

        // Compressed Sparse Row
        public void GetCSR(out int[] ptrs, out int[] indexs, out double[] values)
        {
            ptrs = new int[RowLength + 1];
            IList<int> tmpIndexs = new List<int>();
            IList<double> tmpValues = new List<double>();
            for (int row = 0; row < RowLength; row++)
            {
                ptrs[row] = tmpIndexs.Count;
                foreach (var pair in RowColIndexValues[row])
                {
                    int col = pair.Key;
                    double value = pair.Value;
                    tmpIndexs.Add(col);
                    tmpValues.Add(value);
                }
            }
            ptrs[RowLength] = tmpIndexs.Count;
            indexs = tmpIndexs.ToArray();
            values = tmpValues.ToArray();
        }

        public void SetCSR(int[] ptrs, int[] indexs, double[] values)
        {
            int n = ptrs.Length - 1;
            Resize(n, n);
            for (int row = 0; row < n; row++)
            {
                int sPtr = ptrs[row];
                int ePtr = ptrs[row + 1]; // Note: APtr size is n + 1, APtr[n + 1] is always AIndexsLength
                for (int iPtr = sPtr; iPtr < ePtr; iPtr++)
                {
                    int col = indexs[iPtr];
                    double value = values[iPtr];
                    this[row, col] = value;
                }
            }
        }

        public bool IsSymmetric()
        {
            System.Diagnostics.Debug.Assert(RowLength == ColumnLength);
            bool isSymmetric = true;
            int n = RowLength;
            for (int row = 0; row < n; row++)
            {
                foreach (var pair in RowColIndexValues[row])
                {
                    int col = pair.Key;
                    if (col >= row + 1 && col < n)
                    {
                        //double diff = this[row, col] - this[col, row];
                        double diff = pair.Value - this[col, row];
                        if (Math.Abs(diff) >= IvyFEM.Constants.PrecisionLowerLimit)
                        {
                            isSymmetric = false;
                            break;
                        }
                    }
                }
                if (!isSymmetric)
                {
                    break;
                }
            }
            return isSymmetric;
        }
    }
}
