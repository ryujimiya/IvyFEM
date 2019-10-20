using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Lapack
{
    public class ComplexBandMatrix
    {
        public System.Numerics.Complex[] Buffer { get; protected set; } = null;
        public int RowLength => RowColLength;
        public int ColumnLength => RowColLength;
        public int RowColLength { get; protected set; } = 0;
        public int SubdiaLength { get; protected set; } = 0;
        public int SuperdiaLength { get; protected set; } = 0;

        public ComplexBandMatrix()
        {
            Clear();
        }

        public ComplexBandMatrix(int rowColLength, int subdiaLength, int superdiaLength)
        {
            Resize(rowColLength, subdiaLength, superdiaLength);
        }

        public ComplexBandMatrix(System.Numerics.Complex[] buffer, int rowColLength, int subdiaLength, int superdiaLength, bool alloc = true)
        {
            Copy(buffer, rowColLength, subdiaLength, superdiaLength, true);
        }

        public ComplexBandMatrix(ComplexBandMatrix src)
        {
            Copy(src);
        }

        public ComplexBandMatrix(ComplexMatrix denseM)
        {
            System.Diagnostics.Debug.Assert(denseM.RowLength == denseM.ColumnLength);
            if (denseM.RowLength != denseM.ColumnLength)
            {
                System.Diagnostics.Debug.Assert(false);
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
                        // Note: Magnitue >= 0
                        if (denseM[r, c].Magnitude >= IvyFEM.Constants.PrecisionLowerLimit)
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
                        // Note: Magnitue >= 0
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
            //System.Diagnostics.Debug.WriteLine("rowcolLength: {0} subdiaLength: {1} superdiaLength: {2}",
            //    rowColLength, subdiaLength, superdiaLength);

            // バッファの確保
            Resize(rowColLength, subdiaLength, superdiaLength);
            // 値をコピーする
            for (int c = 0; c < rowColLength; c++)
            {
                // 対角成分
                this[c, c] = denseM[c, c];

                // subdiagonal成分
                if (c < rowColLength - 1)
                {
                    for (int r = c + 1; r <= c + subdiaLength && r < rowColLength; r++)
                    {
                        this[r, c] = denseM[r, c];
                    }
                }
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

        public ComplexBandMatrix(IvyFEM.Linear.ComplexSparseMatrix sparseM)
        {
            System.Diagnostics.Debug.Assert(sparseM.RowLength == sparseM.ColumnLength);
            if (sparseM.RowLength != sparseM.ColumnLength)
            {
                System.Diagnostics.Debug.Assert(false);
            }
            int rowColLength = sparseM.RowLength;

            // subdia長さ、superdia長さを取得する
            int subdiaLength = 0;
            int superdiaLength = 0;
            {
                int tmpRowColLength = 0;
                IvyFEM.Linear.Utils.GetComplexSparseMatrixSubDiaSuperDia(
                    sparseM, out tmpRowColLength, out subdiaLength, out superdiaLength);
            }

            // バッファの確保
            Resize(rowColLength, subdiaLength, superdiaLength);
            // 値をコピーする
            for (int c = 0; c < rowColLength; c++)
            {
                // 対角成分
                this[c, c] = sparseM[c, c];

                // subdiagonal成分
                if (c < rowColLength - 1)
                {
                    for (int r = c + 1; r <= c + subdiaLength && r < rowColLength; r++)
                    {
                        this[r, c] = sparseM[r, c];
                    }
                }
                // superdiagonal成分
                if (c > 0)
                {
                    for (int r = c - 1; r >= c - superdiaLength && r >= 0; r--)
                    {
                        this[r, c] = sparseM[r, c];
                    }
                }
            }
            System.Diagnostics.Debug.WriteLine("cast to band matrix: rowcolLength: {0} subdiaLength: {1} superdiaLength: {2}",
                rowColLength, subdiaLength, superdiaLength);
        }

        public void Resize(int rowColLength, int subdiaLength, int superdiaLength)
        {
            int bufferRowLength;
            int bufferColLength;
            GetBufferRowColLength(out bufferRowLength, out bufferColLength, 
                rowColLength, subdiaLength, superdiaLength);

            Buffer = new System.Numerics.Complex[bufferRowLength * bufferColLength];
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

        public System.Numerics.Complex this[int row, int col]
        {
            get
            {
                if (row < 0 || this.RowLength <= row || col < 0 || this.ColumnLength <= col)
                {
                    throw new IndexOutOfRangeException();
                }
                if (!(row >= col - SuperdiaLength && row <= col + SubdiaLength))
                {
                    return new System.Numerics.Complex();
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

        public void Copy(ComplexBandMatrix src)
        {
            Copy(src.Buffer, src.RowColLength, src.SubdiaLength, src.SuperdiaLength,true);
        }

        public void Copy(System.Numerics.Complex[] buffer,
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
                Buffer[i] = (System.Numerics.Complex)0;
            }
        }

        public void Identity()
        {
            Zero();
            for (int i = 0; i < RowLength; i++)
            {
                this[i, i] = (System.Numerics.Complex)1;
            }
        }

        /*
        public void Transpose()
        {
            throw new NotImplementedException();
        }
        */

        public static ComplexBandMatrix Conjugate(ComplexBandMatrix A)
        {
            System.Numerics.Complex[] x = IvyFEM.Lapack.Functions.zlacgv(A.Buffer);
            ComplexBandMatrix X = new ComplexBandMatrix(x,
                A.RowColLength, A.SubdiaLength, A.SuperdiaLength, false);
            return X;
        }

        public void Conjugate()
        {
            ComplexBandMatrix ret = Conjugate(this);
            Buffer = ret.Buffer;
        }

        /*
        public static ComplexBandMatrix operator *(ComplexBandMatrix A, ComplexBandMatrix B)
        {
            throw new NotImplementedException();
        }
        */

        public static System.Numerics.Complex[] operator *(ComplexBandMatrix A, System.Numerics.Complex[] b)
        {
            System.Numerics.Complex[] c;
            IvyFEM.Lapack.Functions.zgbmvAX(out c,
                A.Buffer, A.RowLength, A.ColumnLength, A.SubdiaLength, A.SubdiaLength, TransposeType.Nop, b);
            return c;
        }
    }
}
