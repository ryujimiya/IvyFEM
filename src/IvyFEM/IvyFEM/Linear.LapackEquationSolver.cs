using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public class BandRepeatSolveInfo
    {
        public int[] IndexsForOrderingToBand { get; set; } = null;
        public int ARow { get; set; } = 0;
        public int ACol { get; set; } = 0;
        public int ASubdia { get; set; } = 0;
        public int ASuperdia { get; set; } = 0;
        public  double[] ADoubleBuffer { get; set; } = null;
        public System.Numerics.Complex[] AComplexBuffer { get; set; } = null;
        public int[] Pivot = null;
    }

    public class LapackEquationSolver : IEquationSolver
    {
        public LapackEquationSolverMethod Method { get; set; } = LapackEquationSolverMethod.Default;
        // バンド行列
        public bool IsOrderingToBandMatrix { get; set; } = false;
        public double[][] CoordsForBandMatrix { get; set; } = null;
        // 繰り返し同じAで解く
        public bool IsRepeatSolve { get; set; } = false;
        private int RepeatCount = 0;
        private BandRepeatSolveInfo BandRepeatSolveInfo = null;

        public bool DoubleSolve(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            bool success = false;
            X = null;

            switch (Method)
            {
                case LapackEquationSolverMethod.Dense:
                    if (IsRepeatSolve)
                    {
                        throw new NotImplementedException();
                    }
                    success = DoubleDenseSolve(out X, A, B);
                    break;

                case LapackEquationSolverMethod.Default:
                case LapackEquationSolverMethod.Band:
                    success = DoubleBandSolve(out X, A, B);
                    break;

                case LapackEquationSolverMethod.PositiveDefiniteBand:
                    if (IsRepeatSolve)
                    {
                        throw new NotImplementedException();
                    }
                    success = DoublePositiveDefiniteBandSolve(out X, A, B);
                    break;

                default:
                    throw new NotImplementedException();
                    //break;
            }

            return success;
        }

        public bool ComplexSolve(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            bool success = false;
            X = null;

            switch (Method)
            {
                case LapackEquationSolverMethod.Dense:
                    if (IsRepeatSolve)
                    {
                        throw new NotImplementedException();
                    }
                    success = ComplexDenseSolve(out X, A, B);
                    break;

                case LapackEquationSolverMethod.Default:
                case LapackEquationSolverMethod.Band:
                    success = ComplexBandSolve(out X, A, B);
                    break;

                case LapackEquationSolverMethod.PositiveDefiniteBand:
                    if (IsRepeatSolve)
                    {
                        throw new NotImplementedException();
                    }
                    success = ComplexPositiveDefiniteBandSolve(out X, A, B);
                    break;

                default:
                    throw new NotImplementedException();
                    //break;
            }

            return success;

        }

        ///////////////////////////////////////////////////////////////////////////////
        // dense
        private bool DoubleDenseSolve(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            bool success = false;
            IvyFEM.Lapack.DoubleMatrix denseA = new IvyFEM.Lapack.DoubleMatrix(A);
            int xRow;
            int xCol;
            int ret = IvyFEM.Lapack.Functions.dgesv(out X, out xRow, out xCol,
                denseA.Buffer, denseA.RowLength, denseA.ColumnLength,
                B, B.Length, 1);
            System.Diagnostics.Debug.Assert(ret == 0);
            success = (ret == 0);
            return success;
        }

        private bool ComplexDenseSolve(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            bool success = false;
            IvyFEM.Lapack.ComplexMatrix denseA = new IvyFEM.Lapack.ComplexMatrix(A);
            int xRow;
            int xCol;
            int ret = IvyFEM.Lapack.Functions.zgesv(out X, out xRow, out xCol,
                denseA.Buffer, denseA.RowLength, denseA.ColumnLength,
                B, B.Length, 1);
            System.Diagnostics.Debug.Assert(ret == 0);
            success = (ret == 0);

            return success;
        }


        ///////////////////////////////////////////////////////////////////////////////
        // band
        private bool DoubleBandSolve(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            bool success = false;
            if (IsRepeatSolve)
            {
                if (RepeatCount == 0)
                {
                    BandRepeatSolveInfo = new BandRepeatSolveInfo();
                    success = DoubleBandOnceSolve(out X, A, B);
                }
                else
                {
                    success = DoubleBandRepeatSolve(out X, B);
                }
                RepeatCount++;
            }
            else
            {
                BandRepeatSolveInfo = null;
                success = DoubleBandOnceSolve(out X, A, B);
            }
            return success;
        }

        private bool DoubleBandOnceSolve(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            X = null;
            DoubleSparseMatrix AToSolve;
            double[] BToSolve;
            int[] indexs;
            if (IsOrderingToBandMatrix)
            {
                bool successOrder = false;
                if (CoordsForBandMatrix != null)
                {
                    // 座標からバンド幅縮小を試みる
                    successOrder = IvyFEM.Linear.Utils.OrderToDoubleBandMatrixByCoord(
                        out AToSolve, out BToSolve, out indexs, A, B, CoordsForBandMatrix);
                }
                else
                {
                    successOrder = IvyFEM.Linear.Utils.OrderToDoubleBandMatrix(
                        out AToSolve, out BToSolve, out indexs, A, B);
                }
                System.Diagnostics.Debug.Assert(successOrder);
                if (!successOrder)
                {
                    return false;
                }
            }
            else
            {
                AToSolve = A;
                BToSolve = B;
                indexs = null;
            }
            if (BandRepeatSolveInfo != null)
            {
                BandRepeatSolveInfo.IndexsForOrderingToBand = indexs;
            }

            bool success = false;
            double[] XToSolve;
            if (IsRepeatSolve)
            {
                success = __DoubleBandSolveDirty(out XToSolve, AToSolve, BToSolve);
            }
            else
            {
                success = __DoubleBandSolve(out XToSolve, AToSolve, BToSolve);
            }
            if (IsOrderingToBandMatrix)
            {
                X = new double[XToSolve.Length];
                for (int row = 0; row < XToSolve.Length; row++)
                {
                    X[indexs[row]] = XToSolve[row];
                }
            }
            else
            {
                X = XToSolve;
            }
            return success;
        }

        private bool __DoubleBandSolve(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            bool success = false;
            IvyFEM.Lapack.DoubleBandMatrix bandA = new IvyFEM.Lapack.DoubleBandMatrix(A);
            A = null;
            int xRow;
            int xCol;
            int ret = IvyFEM.Lapack.Functions.dgbsv(out X, out xRow, out xCol,
                bandA.Buffer, bandA.RowLength, bandA.ColumnLength,
                bandA.SubdiaLength, bandA.SuperdiaLength,
                B, B.Length, 1);
            System.Diagnostics.Debug.Assert(ret == 0);
            success = (ret == 0);
            return success;
        }

        private bool __DoubleBandSolveDirty(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            bool success = false;
            IvyFEM.Lapack.DoubleBandMatrix bandA = new IvyFEM.Lapack.DoubleBandMatrix(A);
            A = null;
            int xRow;
            int xCol;
            int[] pivot;
            int ret = IvyFEM.Lapack.Functions.dgbsv_dirty(out X, out xRow, out xCol,
                bandA.Buffer, bandA.RowLength, bandA.ColumnLength,
                bandA.SubdiaLength, bandA.SuperdiaLength,
                B, B.Length, 1,
                out pivot);
            System.Diagnostics.Debug.Assert(ret == 0);
            success = (ret == 0);

            if (BandRepeatSolveInfo != null)
            {
                BandRepeatSolveInfo.ARow = bandA.RowLength;
                BandRepeatSolveInfo.ACol = bandA.ColumnLength;
                BandRepeatSolveInfo.ASubdia = bandA.SubdiaLength;
                BandRepeatSolveInfo.ASuperdia = bandA.SuperdiaLength;
                BandRepeatSolveInfo.ADoubleBuffer = bandA.Buffer;
                BandRepeatSolveInfo.Pivot = pivot;
            }
            return success;
        }

        private bool DoubleBandRepeatSolve(out double[] X, double[]B)
        {
            bool success = false;
            int xRow;
            int xCol;
            int[] pivot = BandRepeatSolveInfo.Pivot;
            double[] A = BandRepeatSolveInfo.ADoubleBuffer;
            int aRow = BandRepeatSolveInfo.ARow;
            int aCol = BandRepeatSolveInfo.ACol;
            int subdia = BandRepeatSolveInfo.ASubdia;
            int superdia = BandRepeatSolveInfo.ASuperdia;
            int bRow = B.Length;
            double[] BToSolve = null;
            if (IsOrderingToBandMatrix)
            {
                int[] indexs = BandRepeatSolveInfo.IndexsForOrderingToBand;
                BToSolve = new double[B.Length];
                for (int row = 0; row < B.Length; row++)
                {
                    BToSolve[row] = B[indexs[row]];
                }
            }
            else
            {
                BToSolve = B;
            }

            double[] XToSolve; 
            int ret = IvyFEM.Lapack.Functions.dgbtrs(out XToSolve, out xRow, out xCol,
                pivot, A, aRow, aCol, subdia, superdia, BToSolve, bRow, 1);
            System.Diagnostics.Debug.Assert(ret == 0);
            success = (ret == 0);

            if (IsOrderingToBandMatrix)
            {
                int[] indexs = BandRepeatSolveInfo.IndexsForOrderingToBand;
                X = new double[XToSolve.Length];
                for (int row = 0; row < XToSolve.Length; row++)
                {
                    X[indexs[row]] = XToSolve[row];
                }
            }
            else
            {
                X = XToSolve;
            }

            return success;
        }

        private bool ComplexBandSolve(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            bool success = false;
            if (IsRepeatSolve)
            {
                if (RepeatCount == 0)
                {
                    BandRepeatSolveInfo = new BandRepeatSolveInfo();
                    success = ComplexBandOnceSolve(out X, A, B);
                }
                else
                {
                    success = ComplexBandRepeatSolve(out X, B);
                }
                RepeatCount++;
            }
            else
            {
                BandRepeatSolveInfo = null;
                success = ComplexBandOnceSolve(out X, A, B);
            }
            return success;
        }

        private bool ComplexBandOnceSolve(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            X = null;
            ComplexSparseMatrix AToSolve;
            System.Numerics.Complex[] BToSolve;
            int[] indexs;
            if (IsOrderingToBandMatrix)
            {
                bool successOrder = false;
                if (CoordsForBandMatrix != null)
                {
                    // 座標からバンド幅縮小を試みる
                    successOrder = IvyFEM.Linear.Utils.OrderToComplexBandMatrixByCoord(
                        out AToSolve, out BToSolve, out indexs, A, B, CoordsForBandMatrix);
                }
                else
                {
                    successOrder = IvyFEM.Linear.Utils.OrderToComplexBandMatrix(
                        out AToSolve, out BToSolve, out indexs, A, B);
                }
                System.Diagnostics.Debug.Assert(successOrder);
                if (!successOrder)
                {
                    return false;
                }
            }
            else
            {
                AToSolve = A;
                BToSolve = B;
                indexs = null;
            }
            if (BandRepeatSolveInfo != null)
            {
                BandRepeatSolveInfo.IndexsForOrderingToBand = indexs;
            }

            bool success = false;
            System.Numerics.Complex[] XToSolve;
            if (IsRepeatSolve)
            {
                success = __ComplexBandSolveDirty(out XToSolve, AToSolve, BToSolve);
            }
            else
            {
                success = __ComplexBandSolve(out XToSolve, AToSolve, BToSolve);
            }
            if (IsOrderingToBandMatrix)
            {
                X = new System.Numerics.Complex[XToSolve.Length];
                for (int row = 0; row < XToSolve.Length; row++)
                {
                    X[indexs[row]] = XToSolve[row];
                }
            }
            else
            {
                X = XToSolve;
            }
            return success;
        }

        private bool __ComplexBandSolve(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            bool success = false;
            IvyFEM.Lapack.ComplexBandMatrix bandA = new IvyFEM.Lapack.ComplexBandMatrix(A);
            A = null;
            int xRow;
            int xCol;
            int ret = IvyFEM.Lapack.Functions.zgbsv(out X, out xRow, out xCol,
                bandA.Buffer, bandA.RowLength, bandA.ColumnLength,
                bandA.SubdiaLength, bandA.SuperdiaLength,
                B, B.Length, 1);
            System.Diagnostics.Debug.Assert(ret == 0);
            success = (ret == 0);
            return success;
        }

        private bool __ComplexBandSolveDirty(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            bool success = false;
            IvyFEM.Lapack.ComplexBandMatrix bandA = new IvyFEM.Lapack.ComplexBandMatrix(A);
            A = null;
            int xRow;
            int xCol;
            int[] pivot;
            int ret = IvyFEM.Lapack.Functions.zgbsv_dirty(out X, out xRow, out xCol,
                bandA.Buffer, bandA.RowLength, bandA.ColumnLength,
                bandA.SubdiaLength, bandA.SuperdiaLength,
                B, B.Length, 1,
                out pivot);
            System.Diagnostics.Debug.Assert(ret == 0);
            success = (ret == 0);

            if (BandRepeatSolveInfo != null)
            {
                BandRepeatSolveInfo.ARow = bandA.RowLength;
                BandRepeatSolveInfo.ACol = bandA.ColumnLength;
                BandRepeatSolveInfo.ASubdia = bandA.SubdiaLength;
                BandRepeatSolveInfo.ASuperdia = bandA.SuperdiaLength;
                BandRepeatSolveInfo.AComplexBuffer = bandA.Buffer;
                BandRepeatSolveInfo.Pivot = pivot;
            }
            return success;
        }

        private bool ComplexBandRepeatSolve(out System.Numerics.Complex[] X, System.Numerics.Complex[] B)
        {
            bool success = false;
            int xRow;
            int xCol;
            int[] pivot = BandRepeatSolveInfo.Pivot;
            System.Numerics.Complex[] A = BandRepeatSolveInfo.AComplexBuffer;
            int aRow = BandRepeatSolveInfo.ARow;
            int aCol = BandRepeatSolveInfo.ACol;
            int subdia = BandRepeatSolveInfo.ASubdia;
            int superdia = BandRepeatSolveInfo.ASuperdia;
            int bRow = B.Length;
            System.Numerics.Complex[] BToSolve = null;
            if (IsOrderingToBandMatrix)
            {
                int[] indexs = BandRepeatSolveInfo.IndexsForOrderingToBand;
                BToSolve = new System.Numerics.Complex[B.Length];
                for (int row = 0; row < B.Length; row++)
                {
                    BToSolve[row] = B[indexs[row]];
                }
            }
            else
            {
                BToSolve = B;
            }

            System.Numerics.Complex[] XToSolve;
            int ret = IvyFEM.Lapack.Functions.zgbtrs(out XToSolve, out xRow, out xCol,
                pivot, A, aRow, aCol, subdia, superdia, BToSolve, bRow, 1);
            System.Diagnostics.Debug.Assert(ret == 0);
            success = (ret == 0);

            if (IsOrderingToBandMatrix)
            {
                int[] indexs = BandRepeatSolveInfo.IndexsForOrderingToBand;
                X = new System.Numerics.Complex[XToSolve.Length];
                for (int row = 0; row < XToSolve.Length; row++)
                {
                    X[indexs[row]] = XToSolve[row];
                }
            }
            else
            {
                X = XToSolve;
            }

            return success;
        }

        ///////////////////////////////////////////////////////////////////////////////
        // positive definite band
        // Symmetric / Hermitian Positive Definite Band Matrix
        private bool DoublePositiveDefiniteBandSolve(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            X = null;
            double th = 1.0e-12;
            bool isSymmetric = A.AssertSymmetric(th);
            if (!isSymmetric)
            {
                return false;
            }

            bool success = false;
            IvyFEM.Lapack.DoubleSymmetricBandMatrix pbA = new IvyFEM.Lapack.DoubleSymmetricBandMatrix(A);
            A = null;
            int xRow;
            int xCol;
            int ret = IvyFEM.Lapack.Functions.dpbsv(out X, out xRow, out xCol,
                pbA.Buffer, pbA.RowLength, pbA.ColumnLength,
                pbA.SuperdiaLength,
                B, B.Length, 1);
            System.Diagnostics.Debug.Assert(ret == 0);
            success = (ret == 0);
            return success;
        }

        private bool ComplexPositiveDefiniteBandSolve(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            X = null;
            double th = 1.0e-12;
            bool isHermitian = A.AssertHermitian(th);
            if (!isHermitian)
            {
                return false;
            }

            bool success = false;
            IvyFEM.Lapack.ComplexHermitianBandMatrix pbA = new IvyFEM.Lapack.ComplexHermitianBandMatrix(A);
            A = null;
            int xRow;
            int xCol;
            int ret = IvyFEM.Lapack.Functions.zpbsv(out X, out xRow, out xCol,
                pbA.Buffer, pbA.RowLength, pbA.ColumnLength,
                pbA.SuperdiaLength,
                B, B.Length, 1);
            System.Diagnostics.Debug.Assert(ret == 0);
            success = (ret == 0);
            return success;
        }

    }
}
