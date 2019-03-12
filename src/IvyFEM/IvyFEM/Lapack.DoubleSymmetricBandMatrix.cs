using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Lapack
{
    public class DoubleSymmetricBandMatrix
    {
        public double[] Buffer { get; protected set; } = null;
        public int RowLength => RowColLength;
        public int ColumnLength => RowColLength;
        public int RowColLength { get; protected set; } = 0;
        public int SuperdiaLength { get; protected set; } = 0;

        public DoubleSymmetricBandMatrix()
        {
            Clear();
        }

        public DoubleSymmetricBandMatrix(int rowColLength, int superdiaLength)
        {
            Resize(rowColLength, superdiaLength);
        }

        public DoubleSymmetricBandMatrix(double[] buffer, int rowColLength, int superdiaLength, bool alloc = true)
        {
            Copy(buffer, rowColLength, superdiaLength, true);
        }

        public DoubleSymmetricBandMatrix(DoubleSymmetricBandMatrix src)
        {
            Copy(src);
        }

        public static explicit operator DoubleSymmetricBandMatrix(DoubleMatrix denseM)
        {
            DoubleSymmetricBandMatrix m = new DoubleSymmetricBandMatrix();
            System.Diagnostics.Debug.Assert(denseM.RowLength == denseM.ColumnLength);
            if (denseM.RowLength != denseM.ColumnLength)
            {
                System.Diagnostics.Debug.Assert(false);
                return m;
            }
            int rowColLength = denseM.RowLength;

            // superdia長さを取得する
            int superdiaLength = 0;
            for (int c = 0; c < rowColLength; c++)
            {
                if (c > 0)
                {
                    int cnt = 0;
                    for (int r = 0; r <= c - 1; r++)
                    {
                        // 非０要素が見つかったら抜ける
                        if (Math.Abs(denseM[r, c]) >= IvyFEM.Constants.PrecisionLowerLimit)
                        {
                            cnt = c - r;
                            break;
                        }
                    }
                    if (cnt > superdiaLength)
                    {
                        superdiaLength = cnt;
                    }
                }
            }
            System.Diagnostics.Debug.WriteLine("rowcolLength: {0} superdiaLength: {1}",
                rowColLength, superdiaLength);

            // バッファの確保
            m.Resize(rowColLength, superdiaLength);
            // 値をコピーする
            for (int c = 0; c < rowColLength; c++)
            {
                // 対角成分
                m[c, c] = denseM[c, c];

                // superdiagonal成分
                if (c > 0)
                {
                    for (int r = c - 1; r >= c - superdiaLength && r >= 0; r--)
                    {
                        m[r, c] = denseM[r, c];
                    }
                }
            }
            return m;
        }

        public static explicit operator DoubleSymmetricBandMatrix(IvyFEM.Linear.DoubleSparseMatrix sparseM)
        {
            DoubleSymmetricBandMatrix m = new DoubleSymmetricBandMatrix();
            System.Diagnostics.Debug.Assert(sparseM.RowLength == sparseM.ColumnLength);
            if (sparseM.RowLength != sparseM.ColumnLength)
            {
                System.Diagnostics.Debug.Assert(false);
                return m;
            }
            int rowColLength = sparseM.RowLength;

            // superdia長さを取得する
            int superdiaLength = 0;
            for (int c = 0; c < rowColLength; c++)
            {
                if (c > 0)
                {
                    int cnt = 0;
                    for (int r = 0; r <= c - 1; r++)
                    {
                        // 非０要素が見つかったら抜ける
                        if (Math.Abs(sparseM[r, c]) >= IvyFEM.Constants.PrecisionLowerLimit)
                        {
                            cnt = c - r;
                            break;
                        }
                    }
                    if (cnt > superdiaLength)
                    {
                        superdiaLength = cnt;
                    }
                }
            }
            System.Diagnostics.Debug.WriteLine("rowcolLength: {0} superdiaLength: {1}",
                rowColLength, superdiaLength);

            // バッファの確保
            m.Resize(rowColLength, superdiaLength);
            // 値をコピーする
            for (int c = 0; c < rowColLength; c++)
            {
                // 対角成分
                m[c, c] = sparseM[c, c];

                // superdiagonal成分
                if (c > 0)
                {
                    for (int r = c - 1; r >= c - superdiaLength && r >= 0; r--)
                    {
                        m[r, c] = sparseM[r, c];
                    }
                }
            }
            return m;
        }

        public void Resize(int rowColLength, int superdiaLength)
        {
            int bufferRowLength;
            int bufferColLength;
            GetBufferRowColLength(out bufferRowLength, out bufferColLength,
                rowColLength, superdiaLength);

            Buffer = new double[bufferRowLength * bufferColLength];
            RowColLength = rowColLength;
            SuperdiaLength = superdiaLength;
        }

        public void Clear()
        {
            Buffer = null;
            RowColLength = 0;
            SuperdiaLength = 0;
        }

        public string Dump()
        {
            string ret = "";
            string CRLF = System.Environment.NewLine;

            ret += "DoubleSymmetricBandMatrix" + CRLF;
            ret += "ColumnLength = " + ColumnLength + CRLF;
            ret += "RowLength = " + RowLength + CRLF;
            for (int col = 0; col < ColumnLength; col++)
            {
                for (int row = 0; row < RowLength; row++)
                {
                    ret += "[" + row + ", " + col + "] = " + this[row, col] + CRLF;
                }
            }
            return ret;
        }

        private static void GetBufferRowColLength(out int bufferRowLength, out int bufferColLength,
            int rowColLength, int superdiaLength)
        {
            bufferRowLength = superdiaLength + 1;
            bufferColLength = rowColLength;
        }

        public double this[int row, int col]
        {
            get
            {
                if (row < 0 || this.RowLength <= row || col < 0 || this.ColumnLength <= col)
                {
                    throw new IndexOutOfRangeException();
                }
                if (!(row >= col - SuperdiaLength && row <= col + SuperdiaLength))
                {
                    return 0;
                }
                int bufferRowLength;
                int bufferColLength;
                GetBufferRowColLength(out bufferRowLength, out bufferColLength,
                    RowColLength, SuperdiaLength);
                if (row > col)
                {
                    // 下三角形
                    return Buffer[(col - row) + SuperdiaLength + row * bufferRowLength];
                }
                // 上三角形
                return Buffer[(row - col) + SuperdiaLength + col * bufferRowLength];
            }
            set
            {
                if (row < 0 || RowLength <= row || col < 0 || ColumnLength <= col)
                {
                    throw new IndexOutOfRangeException();
                }
                if (!(row >= col - SuperdiaLength && row <= col + SuperdiaLength))
                {
                    return;
                }
                if (row > col)
                {
                    // 下三角形
                    return;
                }
                int bufferRowLength;
                int bufferColLength;
                GetBufferRowColLength(out bufferRowLength, out bufferColLength,
                    RowColLength, SuperdiaLength);

                // 上三角形
                Buffer[(row - col) + SuperdiaLength + col * bufferRowLength] = value;
            }
        }

        public void Copy(DoubleSymmetricBandMatrix src)
        {
            Copy(src.Buffer, src.RowColLength, src.SuperdiaLength, true);
        }

        public void Copy(double[] buffer,
            int rowColLength, int superdiaLength, bool alloc)
        {
            int bufferRowLength;
            int bufferColLength;
            GetBufferRowColLength(out bufferRowLength, out bufferColLength,
                RowColLength, SuperdiaLength);

            System.Diagnostics.Debug.Assert(buffer.Length == bufferRowLength * bufferColLength);
            if (buffer.Length != bufferRowLength * bufferColLength)
            {
                return;
            }

            if (!alloc)
            {
                Buffer = buffer;
                RowColLength = rowColLength;
                SuperdiaLength = superdiaLength;
            }
            else
            {
                // バッファ確保
                if (Buffer != null && Buffer.Length == bufferRowLength * bufferColLength)
                {
                    // なにもしない
                }
                else
                {
                    Resize(rowColLength, superdiaLength);
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
                Buffer[i] = 0;
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

        /*
        public static DoubleSymmetricBandMatrix operator *(DoubleSymmetricBandMatrix A, DoubleSymmetricBandMatrix B)
        {
            throw new NotImplementedException();
        }
        */

        public static double[] operator *(DoubleSymmetricBandMatrix A, double[] b)
        {
            // 最適化する必要がある
            System.Diagnostics.Debug.Assert(A.ColumnLength == b.Length);
            double[] c = new double[b.Length];
            for (int row = 0; row < A.RowLength; row++)
            {
                for (int col = 0; col < A.ColumnLength; col++)
                {
                    c[row] += A[row, col] * b[col];
                }
            }
            return c;
        }
    }
}
