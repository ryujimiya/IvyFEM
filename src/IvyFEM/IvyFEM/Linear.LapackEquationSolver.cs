using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public class LapackEquationSolver : IEquationSolver
    {
        public LapackEquationSolverMethod Method { get; set; } = LapackEquationSolverMethod.Default;
        public bool IsOrderingToBandMatrix { get; set; } = false;

        public bool DoubleSolve(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            bool success = false;
            X = null;

            switch (Method)
            {
                case LapackEquationSolverMethod.Dense:
                    success = DoubleDenseSolve(out X, A, B);
                    break;

                case LapackEquationSolverMethod.Default:
                case LapackEquationSolverMethod.Band:
                    success = DoubleBandSolve(out X, A, B);
                    break;

                case LapackEquationSolverMethod.PositiveDefiniteBand:
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
                    success = ComplexDenseSolve(out X, A, B);
                    break;

                case LapackEquationSolverMethod.Default:
                case LapackEquationSolverMethod.Band:
                    success = ComplexBandSolve(out X, A, B);
                    break;

                case LapackEquationSolverMethod.PositiveDefiniteBand:
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
            IvyFEM.Lapack.DoubleMatrix denseA = (IvyFEM.Lapack.DoubleMatrix)A;
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
            IvyFEM.Lapack.ComplexMatrix denseA = (IvyFEM.Lapack.ComplexMatrix)A;
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
            X = null;
            DoubleSparseMatrix AToSolve;
            double[] BToSolve;
            int[] indexs;
            if (IsOrderingToBandMatrix)
            {
                bool successOrder = IvyFEM.Linear.Utils.OrderToDoubleBandMatrix(
                    out AToSolve, out BToSolve, out indexs, A, B);
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

            double[] XToSolve;
            bool success = __DoubleBandSolve(out XToSolve, AToSolve, BToSolve);
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
            IvyFEM.Lapack.DoubleBandMatrix bandA = (IvyFEM.Lapack.DoubleBandMatrix)A;
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


        private bool ComplexBandSolve(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            X = null;
            ComplexSparseMatrix AToSolve;
            System.Numerics.Complex[] BToSolve;
            int[] indexs;
            if (IsOrderingToBandMatrix)
            {
                bool successOrder = IvyFEM.Linear.Utils.OrderToComplexBandMatrix(
                    out AToSolve, out BToSolve, out indexs, A, B);
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

            System.Numerics.Complex[] XToSolve;
            bool success = __ComplexBandSolve(out XToSolve, AToSolve, BToSolve);
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
            IvyFEM.Lapack.ComplexBandMatrix bandA = (IvyFEM.Lapack.ComplexBandMatrix)A;
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

        ///////////////////////////////////////////////////////////////////////////////
        // positive definite band
        // Symmetric / Hermitian Positive Definite Band Matrix
        private bool DoublePositiveDefiniteBandSolve(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            X = null;
            bool isSymmetric = A.IsSymmetric();
            System.Diagnostics.Debug.Assert(isSymmetric);
            if (!isSymmetric)
            {
                return false;
            }

            bool success = false;
            IvyFEM.Lapack.DoubleSymmetricBandMatrix pbA = (IvyFEM.Lapack.DoubleSymmetricBandMatrix)A;
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
            bool isHermitian = A.IsHermitian();
            System.Diagnostics.Debug.Assert(isHermitian);
            if (!isHermitian)
            {
                return false;
            }

            bool success = false;
            IvyFEM.Lapack.ComplexHermitianBandMatrix pbA = (IvyFEM.Lapack.ComplexHermitianBandMatrix)A;
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
