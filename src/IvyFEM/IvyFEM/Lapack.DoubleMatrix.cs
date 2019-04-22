using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Lapack
{
    /// <summary>
    /// 行列クラス(double)
    ///    1次元配列として行列データを保持します。
    ///    1次元配列は、clapackの配列数値格納順序と同じ（行データを先に格納する: Column Major Order)
    ///    既存のdouble[,]からの置き換えポイント
    ///       double[,] --> DoubleMatrix
    ///       GetLength(0) --> RowLength
    ///       GetLength(1) --> ColumnLength
    /// </summary>
    public class DoubleMatrix
    {
        public double[] Buffer { get; protected set; } = null;
        public int RowLength { get; protected set; } = 0;
        public int ColumnLength { get; protected set; } = 0;

        public DoubleMatrix()
        {
            Clear();
        }

        public DoubleMatrix(int rowLength, int columnLength)
        {
            Resize(rowLength, columnLength);
        }

        public DoubleMatrix(double[] buffer, int rowLength, int columnLength, bool alloc = true)
        {
            Copy(buffer, rowLength, columnLength, alloc);
        }

        public DoubleMatrix(DoubleMatrix src)
        {
            Copy(src);
        }

        public static explicit operator DoubleMatrix(IvyFEM.Linear.DoubleSparseMatrix sparseM)
        {
            DoubleMatrix m = new DoubleMatrix(sparseM.RowLength, sparseM.ColumnLength);
            for (int row = 0; row < sparseM.RowLength; row++)
            {
                var colIndexValues = sparseM.RowColIndexValues[row];
                foreach (var pair in colIndexValues)
                {
                    int col = pair.Key;
                    double value = pair.Value;
                    m[row, col] = value;
                }
            }
            return m;
        }

        public void Resize(int rowLength, int columnLength)
        {
            Buffer = new double[rowLength * columnLength];
            RowLength = rowLength;
            ColumnLength = columnLength;
        }

        public void Clear()
        {
            Buffer = null;
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
                return Buffer[row + col * RowLength];
            }
            set
            {
                if (row < 0 || RowLength <= row || col < 0 || ColumnLength <= col)
                {
                    throw new IndexOutOfRangeException();
                }
                Buffer[row + col * RowLength] = value;
            }
        }

        public void Copy(DoubleMatrix src)
        {
            Copy(src.Buffer, src.RowLength, src.ColumnLength, true);
        }

        public void Copy(double[] buffer, int rowLength, int columnLength, bool alloc)
        {
            System.Diagnostics.Debug.Assert(buffer.Length == rowLength * columnLength);
            if (buffer.Length != rowLength * columnLength)
            {
                return;
            }

            if (!alloc)
            {
                Buffer = buffer;
                RowLength = rowLength;
                ColumnLength = columnLength;
            }
            else
            {
                // バッファ確保
                if (RowLength == rowLength && ColumnLength == columnLength)
                {
                    // 何もしない
                }
                else if (Buffer != null && Buffer.Length == rowLength * columnLength)
                {
                    RowLength = rowLength;
                    ColumnLength = columnLength;
                }
                else
                {
                    Resize(rowLength, columnLength);
                }

                // コピー
                buffer.CopyTo(Buffer, 0);
            }
        }

        public void Zero()
        {
            int size = Buffer.Length;
            for (int i = 0; i < size; i++)
            {
                Buffer[i] = 0.0;
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

        public void Transpose()
        {
            DoubleMatrix t = new DoubleMatrix(ColumnLength, RowLength);

            for (int row = 0; row < RowLength; row++)
            {
                for (int col = 0; col < ColumnLength; col++)
                {
                    t[col, row] = this[row, col];
                }
            }

            Copy(t.Buffer, t.RowLength, t.ColumnLength, false);
        }

        public static DoubleMatrix Inverse(DoubleMatrix A)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;
            DoubleMatrix workA = new DoubleMatrix(A);
            DoubleMatrix workB = new DoubleMatrix(n, n);
            workB.Identity(); // 単位行列
            double[] a = workA.Buffer;
            double[] b = workB.Buffer;
            // [A][X] = [B]
            //  [B]の内容が書き換えられるので、matXを新たに生成せず、Bを出力に指定している
            int xRow = 0;
            int xCol = 0;
            IvyFEM.Lapack.Functions.dgesv(out b, out xRow, out xCol, a, n, n, b, n, n);

            bool alloc = false; // 指定したバッファを使用する
            DoubleMatrix X = new DoubleMatrix(b, xRow, xCol, alloc);
            return X;
        }

        public void Inverse()
        {
            DoubleMatrix ret = Inverse(this);
            Copy(ret.Buffer, ret.RowLength, ret.ColumnLength, false);
        }

        public static DoubleMatrix operator *(DoubleMatrix A, DoubleMatrix B)
        {
            double[] c;
            int cRow;
            int cCol;
            IvyFEM.Lapack.Functions.dgemmAB(out c, out cRow, out cCol,
                A.Buffer, A.RowLength, A.ColumnLength, TransposeType.Nop,
                B.Buffer, B.RowLength, B.ColumnLength, TransposeType.Nop);

            bool alloc = false;
            DoubleMatrix C = new DoubleMatrix(c, cRow, cCol, alloc);
            return C;
        }

        public static double[] operator *(DoubleMatrix A, double[] b)
        {
            double[] c;
            IvyFEM.Lapack.Functions.dgemvAX(out c, 
                A.Buffer, A.RowLength, A.ColumnLength, TransposeType.Nop,
                b);
            return c;
        }

        public bool IsSymmetric()
        {
            System.Diagnostics.Debug.Assert(RowLength == ColumnLength);
            bool isSymmetric = true;
            int n = RowLength;
            for (int row = 0; row < n; row++)
            {
                for (int col = row + 1; col < n; col++)
                {
                    double diff = this[row, col] - this[col, row];
                    if (Math.Abs(diff) >= IvyFEM.Constants.PrecisionLowerLimit)
                    {
                        isSymmetric = false;
                        break;
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
