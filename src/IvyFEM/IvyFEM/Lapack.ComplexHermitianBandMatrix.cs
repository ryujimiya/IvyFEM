using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Lapack
{
    public class ComplexHermitianBandMatrix
    {
        public System.Numerics.Complex[] Buffer { get; protected set; } = null;
        public int RowLength => RowColLength;
        public int ColumnLength => RowColLength;
        public int RowColLength { get; protected set; } = 0;
        public int SuperdiaLength { get; protected set; } = 0;

        public ComplexHermitianBandMatrix()
        {
            Clear();
        }

        public ComplexHermitianBandMatrix(int rowColLength, int superdiaLength)
        {
            Resize(rowColLength, superdiaLength);
        }

        public ComplexHermitianBandMatrix(System.Numerics.Complex[] buffer, int rowColLength, int superdiaLength, bool alloc = true)
        {
            Copy(buffer, rowColLength, superdiaLength, true);
        }

        public ComplexHermitianBandMatrix(ComplexHermitianBandMatrix src)
        {
            Copy(src);
        }

        public ComplexHermitianBandMatrix(ComplexMatrix denseM)
        {
            System.Diagnostics.Debug.Assert(denseM.RowLength == denseM.ColumnLength);
            if (denseM.RowLength != denseM.ColumnLength)
            {
                System.Diagnostics.Debug.Assert(false);
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
                        if (denseM[r, c].Magnitude >= IvyFEM.Constants.PrecisionLowerLimit)
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
            //System.Diagnostics.Debug.WriteLine("rowcolLength: {0} superdiaLength: {1}",
            //    rowColLength, superdiaLength);

            // バッファの確保
            Resize(rowColLength, superdiaLength);
            // 値をコピーする
            for (int c = 0; c < rowColLength; c++)
            {
                // 対角成分
                this[c, c] = denseM[c, c];

                // superdiagonal成分
                if (c > 0)
                {
                    for (int r = c - 1; r >= c - superdiaLength && r >= 0; r--)
                    {
                        this[r, c] = denseM[r, c];
                    }
                }
            }
        }

        public ComplexHermitianBandMatrix(IvyFEM.Linear.ComplexSparseMatrix sparseM)
        {
            ComplexHermitianBandMatrix m = new ComplexHermitianBandMatrix();
            System.Diagnostics.Debug.Assert(sparseM.RowLength == sparseM.ColumnLength);
            if (sparseM.RowLength != sparseM.ColumnLength)
            {
                System.Diagnostics.Debug.Assert(false);
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
                        if (sparseM[r, c].Magnitude >= IvyFEM.Constants.PrecisionLowerLimit)
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
            Resize(rowColLength, superdiaLength);
            // 値をコピーする
            for (int c = 0; c < rowColLength; c++)
            {
                // 対角成分
                this[c, c] = sparseM[c, c];

                // superdiagonal成分
                if (c > 0)
                {
                    for (int r = c - 1; r >= c - superdiaLength && r >= 0; r--)
                    {
                        this[r, c] = sparseM[r, c];
                    }
                }
            }
        }

        public void Resize(int rowColLength, int superdiaLength)
        {
            int bufferRowLength;
            int bufferColLength;
            GetBufferRowColLength(out bufferRowLength, out bufferColLength,
                rowColLength, superdiaLength);

            Buffer = new System.Numerics.Complex[bufferRowLength * bufferColLength];
            RowColLength = rowColLength;
            SuperdiaLength = superdiaLength;
        }

        public void Clear()
        {
            Buffer = null;
            RowColLength = 0;
            SuperdiaLength = 0;
        }

        private static void GetBufferRowColLength(out int bufferRowLength, out int bufferColLength,
            int rowColLength, int superdiaLength)
        {
            bufferRowLength = superdiaLength + 1;
            bufferColLength = rowColLength;
        }

        public System.Numerics.Complex this[int row, int col]
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
                    return System.Numerics.Complex.Conjugate(
                        Buffer[(col - row) + SuperdiaLength + row * bufferRowLength]);
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

        public void Copy(ComplexHermitianBandMatrix src)
        {
            Copy(src.Buffer, src.RowColLength, src.SuperdiaLength, true);
        }

        public void Copy(System.Numerics.Complex[] buffer,
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
        public static ComplexSymmetricBandMatrix operator *(
            ComplexSymmetricBandMatrix A, ComplexSymmetricBandMatrix B)
        {
            throw new NotImplementedException();
        }
        */

        public static System.Numerics.Complex[] operator *(ComplexHermitianBandMatrix A, System.Numerics.Complex[] b)
        {
            // 最適化する必要がある
            System.Diagnostics.Debug.Assert(A.ColumnLength == b.Length);
            System.Numerics.Complex[] c = new System.Numerics.Complex[b.Length];
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
