using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Lapack
{
    public class DoubleBandMatrix
    {
        public double[] Buffer { get; protected set; } = null;
        public int RowLength => RowColLength;
        public int ColumnLength => RowColLength;
        public int RowColLength { get; protected set; } = 0;
        public int SubdiaLength { get; protected set; } = 0;
        public int SuperdiaLength { get; protected set; } = 0;

        public DoubleBandMatrix()
        {
            Clear();
        }

        public DoubleBandMatrix(int rowColLength, int subdiaLength, int superdiaLength)
        {
            Resize(rowColLength, subdiaLength, superdiaLength);
        }

        public DoubleBandMatrix(double[] buffer, int rowColLength, int subdiaLength, int superdiaLength, bool alloc = true)
        {
            Copy(buffer, rowColLength, subdiaLength, superdiaLength, true);
        }

        public DoubleBandMatrix(DoubleBandMatrix src)
        {
            Copy(src);
        }

        public static explicit operator DoubleBandMatrix(DoubleMatrix denseM)
        {
            DoubleBandMatrix m = new DoubleBandMatrix();
            System.Diagnostics.Debug.Assert(denseM.RowLength == denseM.ColumnLength);
            if (denseM.RowLength != denseM.ColumnLength)
            {
                System.Diagnostics.Debug.Assert(false);
                return m;
            }
            int rowColLength = denseM.RowLength;

            // subdia長さ、superdia長さを取得する
            int subdiaLength = 0;
            int superdiaLength = 0;
            for (int c = 0; c < rowColLength; c++)
            {
                if (c < rowColLength - 1)
                {
                    int cnt = 0;
                    for (int r = rowColLength - 1; r >= c + 1; r--)
                    {
                        // 非０要素が見つかったら抜ける
                        if (Math.Abs(denseM[r, c]) >= IvyFEM.Constants.PrecisionLowerLimit)
                        {
                            cnt = r - c;
                            break;
                        }
                    }
                    if (cnt > subdiaLength)
                    {
                        subdiaLength = cnt;
                    }
                }
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
            //System.Diagnostics.Debug.WriteLine("rowcolLength: {0} subdiaLength: {1} superdiaLength: {2}",
            //    rowColLength, subdiaLength, superdiaLength);

            // バッファの確保
            m.Resize(rowColLength, subdiaLength, superdiaLength);
            // 値をコピーする
            for (int c = 0; c < rowColLength; c++)
            {
                // 対角成分
                m[c, c] = denseM[c, c];

                // subdiagonal成分
                if (c < rowColLength - 1)
                {
                    for (int r = c + 1; r <= c + subdiaLength && r < rowColLength; r++)
                    {
                        m[r, c] = denseM[r, c];
                    }
                }
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

        public static explicit operator DoubleBandMatrix(IvyFEM.Linear.DoubleSparseMatrix sparseM)
        {
            DoubleBandMatrix m = new DoubleBandMatrix();
            System.Diagnostics.Debug.Assert(sparseM.RowLength == sparseM.ColumnLength);
            if (sparseM.RowLength != sparseM.ColumnLength)
            {
                System.Diagnostics.Debug.Assert(false);
                return m;
            }
            int rowColLength = sparseM.RowLength;

            // subdia長さ、superdia長さを取得する
            int subdiaLength = 0;
            int superdiaLength = 0;
            {
                int tmpRowColLength = 0;
                IvyFEM.Linear.Utils.GetDoubleSparseMatrixSubDiaSuperDia(
                    sparseM, out tmpRowColLength, out subdiaLength, out superdiaLength);
            }

            // バッファの確保
            m.Resize(rowColLength, subdiaLength, superdiaLength);
            // 値をコピーする
            for (int c = 0; c < rowColLength; c++)
            {
                // 対角成分
                m[c, c] = sparseM[c, c];

                // subdiagonal成分
                if (c < rowColLength - 1)
                {
                    for (int r = c + 1; r <= c + subdiaLength && r < rowColLength; r++)
                    {
                        m[r, c] = sparseM[r, c];
                    }
                }
                // superdiagonal成分
                if (c > 0)
                {
                    for (int r = c - 1; r >= c - superdiaLength && r >= 0; r--)
                    {
                        m[r, c] = sparseM[r, c];
                    }
                }
            }
            System.Diagnostics.Debug.WriteLine("cast to band matrix: rowcolLength: {0} subdiaLength: {1} superdiaLength: {2}",
                rowColLength, subdiaLength, superdiaLength);
            return m;
        }

        public void Resize(int rowColLength, int subdiaLength, int superdiaLength)
        {
            int bufferRowLength;
            int bufferColLength;
            GetBufferRowColLength(out bufferRowLength, out bufferColLength,
                rowColLength, subdiaLength, superdiaLength);

            Buffer = new double[bufferRowLength * bufferColLength];
            RowColLength = rowColLength;
            SubdiaLength = subdiaLength;
            SuperdiaLength = superdiaLength;
        }

        public void Clear()
        {
            Buffer = null;
            RowColLength = 0;
            SubdiaLength = 0;
            SuperdiaLength = 0;
        }

        private static void GetBufferRowColLength(out int bufferRowLength, out int bufferColLength,
            int rowColLength, int subdiaLength, int superdiaLength)
        {
            bufferRowLength = subdiaLength * 2 + superdiaLength + 1;
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
                if (!(row >= col - SuperdiaLength && row <= col + SubdiaLength))
                {
                    return 0;
                }
                int bufferRowLength;
                int bufferColLength;
                GetBufferRowColLength(out bufferRowLength, out bufferColLength,
                    RowColLength, SubdiaLength, SuperdiaLength);

                return Buffer[(row - col) + SubdiaLength + SuperdiaLength + col * bufferRowLength];
            }
            set
            {
                if (row < 0 || RowLength <= row || col < 0 || ColumnLength <= col)
                {
                    throw new IndexOutOfRangeException();
                }
                if (!(row >= col - SuperdiaLength && row <= col + SubdiaLength))
                {
                    return;
                }
                int bufferRowLength;
                int bufferColLength;
                GetBufferRowColLength(out bufferRowLength, out bufferColLength,
                    RowColLength, SubdiaLength, SuperdiaLength);

                Buffer[(row - col) + SubdiaLength + SuperdiaLength + col * bufferRowLength] = value;
            }
        }

        public void Copy(DoubleBandMatrix src)
        {
            Copy(src.Buffer, src.RowColLength, src.SubdiaLength, src.SuperdiaLength, true);
        }

        public void Copy(double[] buffer,
            int rowColLength, int subdiaLength, int superdiaLength, bool alloc)
        {
            int bufferRowLength;
            int bufferColLength;
            GetBufferRowColLength(out bufferRowLength, out bufferColLength,
                RowColLength, SubdiaLength, SuperdiaLength);

            System.Diagnostics.Debug.Assert(buffer.Length == bufferRowLength * bufferColLength);
            if (buffer.Length != bufferRowLength * bufferColLength)
            {
                return;
            }

            if (!alloc)
            {
                Buffer = buffer;
                RowColLength = rowColLength;
                SubdiaLength = subdiaLength;
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
                    Resize(rowColLength, subdiaLength, superdiaLength);
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
        public static DoubleBandMatrix operator *(DoubleBandMatrix A, DoubleBandMatrix B)
        {
            throw new NotImplementedException();
        }
        */

        public static double[] operator *(DoubleBandMatrix A, double[] b)
        {
            double[] c;
            IvyFEM.Lapack.Functions.dgbmvAX(out c,
                A.Buffer, A.RowLength, A.ColumnLength, A.SubdiaLength, A.SubdiaLength, TransposeType.Nop,
                b);
            return c;
        }
    }
}
