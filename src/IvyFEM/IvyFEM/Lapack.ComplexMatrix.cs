using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Lapack
{
    public class ComplexMatrix
    {
        public System.Numerics.Complex[] Buffer { get; protected set; } = null;
        public int RowLength { get; protected set; } = 0;
        public int ColumnLength { get; protected set; } = 0;

        public ComplexMatrix()
        {
            Clear();
        }

        public ComplexMatrix(int rowLength, int columnLength)
        {
            Resize(rowLength, columnLength);
        }

        public ComplexMatrix(System.Numerics.Complex[] buffer, int rowLength, int columnLength, bool alloc = true)
        {
            Copy(buffer, rowLength, columnLength, alloc);
        }

        public ComplexMatrix(ComplexMatrix src)
        {
            Copy(src);
        }

        public static explicit operator ComplexMatrix(DoubleMatrix doubleM)
        {
            ComplexMatrix m = new ComplexMatrix(doubleM.RowLength, doubleM.ColumnLength);
            for (int i = 0; i < m.Buffer.Length; i++)
            {
                m.Buffer[i] = (System.Numerics.Complex)doubleM.Buffer[i];
            }
            return m;
        }

        public static explicit operator ComplexMatrix(IvyFEM.Linear.ComplexSparseMatrix sparseM)
        {
            ComplexMatrix m = new ComplexMatrix(sparseM.RowLength, sparseM.ColumnLength);
            for (int row = 0; row < sparseM.RowLength; row++)
            {
                var sparseColIndexValues = sparseM.RowColIndexValues[row];
                foreach (var pair in sparseColIndexValues)
                {
                    int col = pair.Key;
                    System.Numerics.Complex value = pair.Value;
                    m[row, col] = value;
                }
            }
            return m;
        }

        public void Resize(int rowLength, int columnLength)
        {
            Buffer = new System.Numerics.Complex[rowLength * columnLength];
            RowLength = rowLength;
            ColumnLength = columnLength;
        }

        public void Clear()
        {
            Buffer = null;
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

        public void Copy(ComplexMatrix src)
        {
            Copy(src.Buffer, src.RowLength, src.ColumnLength, true);
        }

        public void Copy(System.Numerics.Complex[] buffer, int rowLength, int columnLength, bool alloc)
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

        public void Transpose()
        {
            ComplexMatrix t = new ComplexMatrix(ColumnLength, RowLength);

            for (int row = 0; row < RowLength; row++)
            {
                for (int col = 0; col < ColumnLength; col++)
                {
                    t[col, row] = this[row, col];
                }
            }
            Copy(t.Buffer, t.RowLength, t.ColumnLength, false);
        }

        public static ComplexMatrix Conjugate(ComplexMatrix A)
        {
            System.Numerics.Complex[] x = IvyFEM.Lapack.Functions.zlacgv(A.Buffer);
            ComplexMatrix X = new ComplexMatrix(x, A.RowLength, A.ColumnLength, false);
            return X;
        }

        public void Conjugate()
        {
            ComplexMatrix ret = Conjugate(this);
            Buffer = ret.Buffer;
        }

        public static ComplexMatrix Inverse(ComplexMatrix A)
        {
            System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
            int n = A.RowLength;
            ComplexMatrix workA = new ComplexMatrix(A);
            ComplexMatrix workB = new ComplexMatrix(n, n);
            workB.Identity(); // 単位行列
            System.Numerics.Complex[] a = workA.Buffer;
            System.Numerics.Complex[] b = workB.Buffer;
            // [A][X] = [B]
            //  [B]の内容が書き換えられるので、matXを新たに生成せず、Bを出力に指定している
            int xRow = 0;
            int xCol = 0;
            IvyFEM.Lapack.Functions.zgesv(out b, out xRow, out xCol, a, n, n, b, n, n);

            bool alloc = false; // 指定したバッファを使用する
            ComplexMatrix X = new ComplexMatrix(b, xRow, xCol, alloc);
            return X;
        }

        public void Inverse()
        {
            ComplexMatrix ret = Inverse(this);
            Copy(ret.Buffer, ret.RowLength, ret.ColumnLength, false);
        }

        public static ComplexMatrix operator *(ComplexMatrix A, ComplexMatrix B)
        {
            System.Numerics.Complex[] c;
            int cRow;
            int cCol;
            IvyFEM.Lapack.Functions.zgemmAB(out c, out cRow, out cCol,
                A.Buffer, A.RowLength, A.ColumnLength, TransposeType.Nop,
                B.Buffer, B.RowLength, B.ColumnLength, TransposeType.Nop);

            bool alloc = false;
            ComplexMatrix C = new ComplexMatrix(c, cRow, cCol, alloc);
            return C;
        }

        public static System.Numerics.Complex[] operator *(ComplexMatrix A, System.Numerics.Complex[] b)
        {
            System.Numerics.Complex[] c;
            IvyFEM.Lapack.Functions.zgemvAX(out c,
                A.Buffer, A.RowLength, A.ColumnLength, TransposeType.Nop, b);
            return c;
        }

        public bool IsHermitian()
        {
            System.Diagnostics.Debug.Assert(RowLength == ColumnLength);
            bool isHermitian = true;
            int n = RowLength;
            for (int row = 0; row < n; row++)
            {
                for (int col = row + 1; col < n; col++)
                {
                    System.Numerics.Complex diff =
                        this[row, col] - System.Numerics.Complex.Conjugate(this[col, row]);
                    if (diff.Magnitude >= IvyFEM.Constants.PrecisionLowerLimit)
                    {
                        isHermitian = false;
                        break;
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
                for (int col = row + 1; col < n; col++)
                {
                    System.Numerics.Complex diff = this[row, col] - this[col, row];
                    if (diff.Magnitude >= IvyFEM.Constants.PrecisionLowerLimit)
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
