using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public class ComplexSparseMatrix
    {
        public int RowLength { get; protected set; } = 0;
        public int ColumnLength { get; protected set; } = 0;
        public Dictionary<int, System.Numerics.Complex>[] RowColIndexValues { get; protected set; } = null;

        public ComplexSparseMatrix()
        {
            Clear();
        }

        public ComplexSparseMatrix(int rowLength, int columnLength)
        {
            Resize(rowLength, columnLength);
        }

        public ComplexSparseMatrix(ComplexSparseMatrix src)
        {
            Copy(src);
        }

        public ComplexSparseMatrix(IvyFEM.Lapack.ComplexMatrix denseM)
        {
            Resize(denseM.RowLength, denseM.ColumnLength);
            for (int row = 0; row < denseM.RowLength; row++)
            {
                for (int col = 0; col < denseM.ColumnLength; col++)
                {
                    if (denseM[row, col].Magnitude >= IvyFEM.Constants.PrecisionLowerLimit)
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
            RowColIndexValues = new Dictionary<int, System.Numerics.Complex>[RowLength];
            for (int row = 0; row < RowLength; row++)
            {
                RowColIndexValues[row] = new Dictionary<int, System.Numerics.Complex>();
            }
        }

        public void Clear()
        {
            RowColIndexValues = null;
            RowLength = 0;
            ColumnLength = 0;
        }

        public System.Numerics.Complex this[int row, int col]
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
                System.Numerics.Complex value = 0;
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
                    if (value.Magnitude >= IvyFEM.Constants.PrecisionLowerLimit)
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
                    if (value.Magnitude >= IvyFEM.Constants.PrecisionLowerLimit)
                    {
                        RowColIndexValues[row].Add(col, value);
                    }
                }
            }
        }

        public void Copy(ComplexSparseMatrix src)
        {
            Resize(src.RowLength, src.ColumnLength);
            for (int row = 0; row < src.RowLength; row++)
            {
                var srcColIndexValues = src.RowColIndexValues[row];
                foreach (var pair in srcColIndexValues)
                {
                    int col = pair.Key;
                    System.Numerics.Complex value = pair.Value;
                    RowColIndexValues[row].Add(col, value);
                }
            }
        }

        public void Zero()
        {
            for (int row = 0; row < ColumnLength; row++)
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

        public static ComplexSparseMatrix operator *(ComplexSparseMatrix A, ComplexSparseMatrix B)
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
            ComplexSparseMatrix C = new ComplexSparseMatrix(cRow, cCol);

            for (int i = 0; i < aRow; i++)
            {
                var colIndexValues1 = A.RowColIndexValues[i];
                for (int j = 0; j < bCol; j++)
                {
                    foreach (var pair1 in colIndexValues1)
                    {
                        int k = pair1.Key;
                        System.Numerics.Complex value1 = pair1.Value; // Aik

                        if (B.RowColIndexValues[k].ContainsKey(j))
                        {
                            System.Numerics.Complex value2 = B.RowColIndexValues[k][j]; // Bkj
                            C[i, j] += value1 * value2;
                        }
                    }

                }
            }
            return C;
        }

        public static System.Numerics.Complex[] operator *(ComplexSparseMatrix A, System.Numerics.Complex[] b)
        {
            System.Numerics.Complex[] c = new System.Numerics.Complex[A.RowLength];
            for (int row = 0; row < A.RowLength; row++)
            {
                var colIndexValues = A.RowColIndexValues[row];
                foreach (var pair in colIndexValues)
                {
                    int col = pair.Key;
                    System.Numerics.Complex value = pair.Value;
                    c[row] += value * b[col];
                }
            }
            return c;
        }

        public static System.Numerics.Complex[] operator *(ComplexSparseMatrix A, double[] b)
        {
            System.Numerics.Complex[] c = new System.Numerics.Complex[A.RowLength];
            for (int row = 0; row < A.RowLength; row++)
            {
                var colIndexValues = A.RowColIndexValues[row];
                foreach (var pair in colIndexValues)
                {
                    int col = pair.Key;
                    System.Numerics.Complex value = pair.Value;
                    c[row] += value * b[col];
                }
            }
            return c;
        }

        // Compressed Sparse Row
        public void GetCSR(out int[] ptrs, out int[] indexs, out System.Numerics.Complex[] values)
        {
            ptrs = new int[RowLength + 1];
            IList<int> tmpIndexs = new List<int>();
            IList<System.Numerics.Complex> tmpValues = new List<System.Numerics.Complex>();
            for (int row = 0; row < RowLength; row++)
            {
                ptrs[row] = tmpIndexs.Count;
                foreach (var pair in RowColIndexValues[row])
                {
                    int col = pair.Key;
                    System.Numerics.Complex value = pair.Value;
                    tmpIndexs.Add(col);
                    tmpValues.Add(value);
                }
            }
            ptrs[RowLength] = tmpIndexs.Count;
            indexs = tmpIndexs.ToArray();
            values = tmpValues.ToArray();
        }

        public void SetCSR(int[] ptrs, int[] indexs, System.Numerics.Complex[] values)
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
                    System.Numerics.Complex value = values[iPtr];
                    this[row, col] = value;
                }
            }
        }

        public bool IsHermitian()
        {
            System.Diagnostics.Debug.Assert(RowLength == ColumnLength);
            bool isHermitian = true;
            int n = RowLength;
            for (int row = 0; row < n; row++)
            {
                foreach (var pair in RowColIndexValues[row])
                {
                    int col = pair.Key;
                    if (col >= row + 1 && col < n)
                    {
                        //System.Numerics.Complex diff =
                        //    this[row, col] - System.Numerics.Complex.Conjugate(this[col, row]);
                        System.Numerics.Complex diff =
                            pair.Value - System.Numerics.Complex.Conjugate(this[col, row]);
                        if (diff.Magnitude >= IvyFEM.Constants.PrecisionLowerLimit)
                        {
                            isHermitian = false;
                            break;
                        }
                    }
                }
                if (!isHermitian)
                {
                    break;
                }
            }
            return isHermitian;
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
                        //System.Numerics.Complex diff = this[row, col] - this[col, row];
                        System.Numerics.Complex diff = pair.Value - this[col, row];
                        if (diff.Magnitude >= IvyFEM.Constants.PrecisionLowerLimit)
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
