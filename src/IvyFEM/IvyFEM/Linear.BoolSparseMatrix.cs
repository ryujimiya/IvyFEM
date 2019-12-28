using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public class BoolSparseMatrix
    {
        public int RowLength { get; protected set; } = 0;
        public int ColumnLength { get; protected set; } = 0;
        public Dictionary<int, bool>[] RowColIndexValues { get; protected set; } = null;

        public BoolSparseMatrix()
        {
            Clear();
        }

        public BoolSparseMatrix(int rowLength, int columnLength)
        {
            Resize(rowLength, columnLength);
        }

        public BoolSparseMatrix(BoolSparseMatrix src)
        {
            Copy(src);
        }

        public void Resize(int rowLength, int columnLength)
        {
            RowLength = rowLength;
            ColumnLength = columnLength;
            RowColIndexValues = new Dictionary<int, bool>[RowLength];
            for (int row = 0; row < RowLength; row++)
            {
                RowColIndexValues[row] = new Dictionary<int, bool>();
            }
        }

        public void Clear()
        {
            RowColIndexValues = null;
            RowLength = 0;
            ColumnLength = 0;
        }

        public bool this[int row, int col]
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
                bool value = false;
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
                    if (value)
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
                    if (value)
                    {
                        RowColIndexValues[row].Add(col, value);
                    }
                }
            }
        }

        public void Copy(BoolSparseMatrix src)
        {
            Resize(src.RowLength, src.ColumnLength);
            for (int row = 0; row < src.RowLength; row++)
            {
                var srcColIndexValues = src.RowColIndexValues[row];
                foreach (var pair in srcColIndexValues)
                {
                    int col = pair.Key;
                    bool value = pair.Value;
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

        /*
        public void Identity()
        {
            throw new NotImplementedException();
        }
        */

        /*
        public void Transpose()
        {
            throw new NotImplementedException();
        }
        */

        public bool AssertSymmetric(double th)
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
                        //if (this[row, col] != this[col, row])
                        if (pair.Value != this[col, row])
                        {
                            isSymmetric = false;
                            System.Diagnostics.Debug.Assert(false);
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
