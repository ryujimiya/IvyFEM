using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using LisScalar = System.Numerics.Complex;

namespace IvyFEM.Lis
{
    unsafe public class LisMatrix : IDisposable
    {
        internal NativeLisMatrix* Native = null;

        public LisMatrix(int comm = IvyFEM.Lis.Constants.LisCommWorld)
        {
            int ret = IvyFEM.Lis.Functions.MatrixCreate(comm, this);
            System.Diagnostics.Debug.Assert(ret == 0);
        }

        #region IDisposable Support
        private bool disposedValue = false; // 重複する呼び出しを検出するには

        protected virtual void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                if (disposing)
                {

                }

                int ret = IvyFEM.Lis.Functions.MatrixDestroy(this);
                System.Diagnostics.Debug.Assert(ret == 0);

                disposedValue = true;
            }
        }

        ~LisMatrix()
        {
            Dispose(false);
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }
        #endregion

        public static explicit operator LisMatrix(IvyFEM.Lapack.DoubleMatrix denseM)
        {
            LisMatrix A = new LisMatrix();
            int ret;
            System.Diagnostics.Debug.Assert(denseM.RowLength == denseM.ColumnLength);
            int n = denseM.RowLength;
            ret = A.SetSize(0, n);
            System.Diagnostics.Debug.Assert(ret == 0);
            for (int row = 0; row < n; row++)
            {
                for (int col = 0; col < n; col++)
                {
                    LisScalar value = denseM[row, col];
                    if (value.Magnitude < IvyFEM.Constants.PrecisionLowerLimit)
                    {
                        continue;
                    }
                    ret = A.SetValue(SetValueFlag.LisInsValue, col, row, value);
                    System.Diagnostics.Debug.Assert(ret == 0);
                }
            }
            ret = A.SetType(IvyFEM.Lis.MatrixType.LisMatrixCSR);
            System.Diagnostics.Debug.Assert(ret == 0);
            ret = A.Assemble();
            System.Diagnostics.Debug.Assert(ret == 0);

            return A;
        }

        public static explicit operator LisMatrix(IvyFEM.Linear.DoubleSparseMatrix sparseM)
        {
            LisMatrix A = new LisMatrix();
            int ret;
            System.Diagnostics.Debug.Assert(sparseM.RowLength == sparseM.ColumnLength);
            int n = sparseM.RowLength;
            ret = A.SetSize(0, n);
            System.Diagnostics.Debug.Assert(ret == 0);
            for (int row = 0; row < n; row++)
            {
                for (int col = 0; col < n; col++)
                {
                    LisScalar value = sparseM[row, col];
                    if (value.Magnitude < IvyFEM.Constants.PrecisionLowerLimit)
                    {
                        continue;
                    }
                    ret = A.SetValue(SetValueFlag.LisInsValue, col, row, value);
                    System.Diagnostics.Debug.Assert(ret == 0);
                }
            }
            ret = A.SetType(IvyFEM.Lis.MatrixType.LisMatrixCSR);
            System.Diagnostics.Debug.Assert(ret == 0);
            ret = A.Assemble();
            System.Diagnostics.Debug.Assert(ret == 0);

            return A;
        }

        public static explicit operator LisMatrix(IvyFEM.Lapack.ComplexMatrix denseM)
        {
            LisMatrix A = new LisMatrix();
            int ret;
            System.Diagnostics.Debug.Assert(denseM.RowLength == denseM.ColumnLength);
            int n = denseM.RowLength;
            ret = A.SetSize(0, n);
            System.Diagnostics.Debug.Assert(ret == 0);
            for (int row = 0; row < n; row++)
            {
                for (int col = 0; col < n; col++)
                {
                    LisScalar value = denseM[row, col];
                    if (value.Magnitude < IvyFEM.Constants.PrecisionLowerLimit)
                    {
                        continue;
                    }
                    ret = A.SetValue(SetValueFlag.LisInsValue, col, row, value);
                    System.Diagnostics.Debug.Assert(ret == 0);
                }
            }
            ret = A.SetType(IvyFEM.Lis.MatrixType.LisMatrixCSR);
            System.Diagnostics.Debug.Assert(ret == 0);
            ret = A.Assemble();
            System.Diagnostics.Debug.Assert(ret == 0);

            return A;
        }

        public static explicit operator LisMatrix(IvyFEM.Linear.ComplexSparseMatrix sparseM)
        {
            LisMatrix A = new LisMatrix();
            int ret;
            System.Diagnostics.Debug.Assert(sparseM.RowLength == sparseM.ColumnLength);
            int n = sparseM.RowLength;
            ret = A.SetSize(0, n);
            System.Diagnostics.Debug.Assert(ret == 0);
            for (int row = 0; row < n; row++)
            {
                for (int col = 0; col < n; col++)
                {
                    LisScalar value = sparseM[row, col];
                    if (value.Magnitude < IvyFEM.Constants.PrecisionLowerLimit)
                    {
                        continue;
                    }
                    ret = A.SetValue(SetValueFlag.LisInsValue, col, row, value);
                    System.Diagnostics.Debug.Assert(ret == 0);
                }
            }
            ret = A.SetType(IvyFEM.Lis.MatrixType.LisMatrixCSR);
            System.Diagnostics.Debug.Assert(ret == 0);
            ret = A.Assemble();
            System.Diagnostics.Debug.Assert(ret == 0);

            return A;
        }

        public int Assemble()
        {
            return IvyFEM.Lis.Functions.MatrixAssemble(this);
        }

        public int SetSize(int localN, int globalN)
        {
            return IvyFEM.Lis.Functions.MatrixSetSize(this, localN, globalN);
        }

        public int GetSize(out int localN, out int globalN)
        {
            return IvyFEM.Lis.Functions.MatrixGetSize(this, out localN, out globalN);
        }

        public int GetRange(out int @is, out int ie)
        {
            return IvyFEM.Lis.Functions.MatrixGetRange(this, out @is, out ie);
        }

        public int SetType(MatrixType matrixType)
        {
            return IvyFEM.Lis.Functions.MatrixSetType(this, matrixType);
        }

        public int GetType(out MatrixType matrixType)
        {
            return IvyFEM.Lis.Functions.MatrixGetType(this, out matrixType);
        }

        public int SetValue(SetValueFlag flag, int i, int j, LisScalar value)
        {
            return IvyFEM.Lis.Functions.MatrixSetValue(flag, i, j, value, this);
        }

        public int SetValueNew(SetValueFlag flag, int i, int j, LisScalar value)
        {
            return IvyFEM.Lis.Functions.MatrixSetValueNew(flag, i, j, value, this);
        }

        public int SetValues(SetValueFlag flag, int n, LisScalar[] values)
        {
            return IvyFEM.Lis.Functions.MatrixSetValues(flag, n, values, this);
        }

        public static int Matvec(LisMatrix A, LisVector x, LisVector y)
        {
            return IvyFEM.Lis.Functions.Matvec(A, x, y);
        }
    }
}
