using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public class IvyFEMEquationSolver : IEquationSolver
    {
        public IvyFEMEquationSolverMethod Method { get; set; } = IvyFEMEquationSolverMethod.Default;
        public double ConvRatioTolerance { get; set; } = IvyFEM.Linear.Constants.ConvRatioTolerance;
        public int ILUFillinLevel { get; set; } = 0;

        public bool DoubleSolve(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            bool success = false;
            X = null;

            switch (Method)
            {
                case IvyFEMEquationSolverMethod.Default:
                case IvyFEMEquationSolverMethod.NoPreconCG:
                    success = DoubleSolveNoPreconCG(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.CG:
                    success = DoubleSolveCG(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.CGWithPivoting:
                    success = DoubleSolveCGWithPivoting(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.ICCG:
                    success = DoubleSolveICCG(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.NoPreconBiCGSTAB:
                    success = DoubleSolveNoPreconBiCGSTAB(out X, A, B);
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
                case IvyFEMEquationSolverMethod.CG:
                    // エルミート行列の場合はCG
                    // 必要なケースがでてくれば実装する
                    throw new NotImplementedException();
                    //break;

                case IvyFEMEquationSolverMethod.Default:
                case IvyFEMEquationSolverMethod.NoPreconCOCG:
                    success = ComplexSolveNoPreconCOCG(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.COCG:
                    // 複素対称行列の場合はCOCG
                    success = ComplexSolveCOCG(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.ICCOCG:
                    success = ComplexSolveICCOCG(out X, A, B);
                    break;

                case IvyFEMEquationSolverMethod.NoPreconBiCGSTAB:
                    success = ComplexSolveNoPreconBiCGSTAB(out X, A, B);
                    break;

                default:
                    throw new NotImplementedException();
                    //break;
            }

            return success;
        }

        //////////////////////////////////////////////////////////////////////////////////////
        // double

        private bool DoubleSolveNoPreconCG(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            int t;
            X = null;
            //t = System.Environment.TickCount;
            //bool isSymmetric = A.IsSymmetric();
            //System.Diagnostics.Debug.Assert(isSymmetric);
            //System.Diagnostics.Debug.WriteLine("  IsSymmetric t = " + (System.Environment.TickCount - t));
            //if (!isSymmetric)
            //{
            //    return false;
            //}

            //t = System.Environment.TickCount;
            //bool success = IvyFEM.Linear.Functions.DoubleSolveNoPreconCG(out X, A, B,
            //    ConvRatioTolerance);
            //System.Diagnostics.Debug.Assert(success);
            //System.Diagnostics.Debug.WriteLine("  DoubleSolveNoPreconCG t = " + (System.Environment.TickCount - t));

            // Nativeを使う
            bool success = false;
            {
                System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
                int n = A.RowLength;
                int[] APtrs;
                int[] AIndexs;
                double[] AValues;
                t = System.Environment.TickCount;
                A.GetCSR(out APtrs, out AIndexs, out AValues);
                System.Diagnostics.Debug.WriteLine("  GetCSR t = " + (System.Environment.TickCount - t));
                t = System.Environment.TickCount;
                success = IvyFEM.Native.Functions.DoubleSolveNoPreconCG(
                    out X, n, APtrs, AIndexs, AValues, B,
                    ConvRatioTolerance);
                System.Diagnostics.Debug.Assert(success);
                System.Diagnostics.Debug.WriteLine("  DoubleSolveNoPreconCG t = " + (System.Environment.TickCount - t));
            }

            return success;
        }

        private bool DoubleSolveCG(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            int t;
            X = null;
            //t = System.Environment.TickCount;
            //bool isSymmetric = A.IsSymmetric();
            //System.Diagnostics.Debug.Assert(isSymmetric);
            //System.Diagnostics.Debug.WriteLine("  IsSymmetric t = " + (System.Environment.TickCount - t));
            //if (!isSymmetric)
            //{
            //    return false;
            //}

            //t = System.Environment.TickCount;
            //bool success = IvyFEM.Linear.Functions.DoubleSolveCG(out X, A, B, ILUFillinLevel,
            //    ConvRatioTolerance);
            //System.Diagnostics.Debug.Assert(success);
            //System.Diagnostics.Debug.WriteLine("  DoubleSolveCG t = " + (System.Environment.TickCount - t));

            // Nativeを使う
            bool success = false;
            {
                System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
                int n = A.RowLength;
                int[] APtrs;
                int[] AIndexs;
                double[] AValues;
                t = System.Environment.TickCount;
                A.GetCSR(out APtrs, out AIndexs, out AValues);
                System.Diagnostics.Debug.WriteLine("  GetCSR t = " + (System.Environment.TickCount - t));
                t = System.Environment.TickCount;
                success = IvyFEM.Native.Functions.DoubleSolveCG(
                    out X, n, APtrs, AIndexs, AValues, B, ILUFillinLevel,
                    ConvRatioTolerance);
                System.Diagnostics.Debug.Assert(success);
                System.Diagnostics.Debug.WriteLine("  DoubleSolveCG t = " + (System.Environment.TickCount - t));
            }

            return success;
        }

        private bool DoubleSolveCGWithPivoting(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            int t;
            X = null;
            //t = System.Environment.TickCount;
            //bool isSymmetric = A.IsSymmetric();
            //System.Diagnostics.Debug.Assert(isSymmetric);
            //System.Diagnostics.Debug.WriteLine("  IsSymmetric t = " + (System.Environment.TickCount - t));
            //if (!isSymmetric)
            //{
            //    return false;
            //}

            //t = System.Environment.TickCount;
            //bool success = IvyFEM.Linear.Functions.DoubleSolveCGWithPivoting(out X, A, B, ILUFillinLevel,
            //    ConvRatioTolerance);
            //System.Diagnostics.Debug.Assert(success);
            //System.Diagnostics.Debug.WriteLine("  DoubleSolveCGWithPivoting t = " + (System.Environment.TickCount - t));

            // Nativeを使う
            bool success = false;
            {
                System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
                int n = A.RowLength;
                int[] APtrs;
                int[] AIndexs;
                double[] AValues;
                t = System.Environment.TickCount;
                A.GetCSR(out APtrs, out AIndexs, out AValues);
                System.Diagnostics.Debug.WriteLine("  GetCSR t = " + (System.Environment.TickCount - t));
                t = System.Environment.TickCount;
                success = IvyFEM.Native.Functions.DoubleSolveCGWithPivoting(
                    out X, n, APtrs, AIndexs, AValues, B, ILUFillinLevel,
                    ConvRatioTolerance);
                System.Diagnostics.Debug.Assert(success);
                System.Diagnostics.Debug.WriteLine("  DoubleSolveCGWithPivoting t = " + (System.Environment.TickCount - t));
            }

            return success;
        }

        private bool DoubleSolveICCG(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            int t;
            X = null;
            //t = System.Environment.TickCount;
            //bool isSymmetric = A.IsSymmetric();
            //System.Diagnostics.Debug.Assert(isSymmetric);
            //System.Diagnostics.Debug.WriteLine("  IsSymmetric t = " + (System.Environment.TickCount - t));
            //if (!isSymmetric)
            //{
            //    return false;
            //}

            //t = System.Environment.TickCount;
            //bool success = IvyFEM.Linear.Functions.DoubleSolveICCG(out X, A, B,
            //    ConvRatioTolerance);
            //System.Diagnostics.Debug.Assert(success);
            //System.Diagnostics.Debug.WriteLine("  DoubleSolveICCG t = " + (System.Environment.TickCount - t));

            // Nativeを使う
            bool success = false;
            {
                System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
                int n = A.RowLength;
                int[] APtrs;
                int[] AIndexs;
                double[] AValues;
                t = System.Environment.TickCount;
                A.GetCSR(out APtrs, out AIndexs, out AValues);
                System.Diagnostics.Debug.WriteLine("  GetCSR t = " + (System.Environment.TickCount - t));
                t = System.Environment.TickCount;
                success = IvyFEM.Native.Functions.DoubleSolveICCG(
                    out X, n, APtrs, AIndexs, AValues, B,
                    ConvRatioTolerance);
                System.Diagnostics.Debug.Assert(success);
                System.Diagnostics.Debug.WriteLine("  DoubleSolveICCG t = " + (System.Environment.TickCount - t));
            }

            return success;
        }

        private bool DoubleSolveNoPreconBiCGSTAB(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            int t;
            X = null;

            //t = System.Environment.TickCount;
            //bool success = IvyFEM.Linear.Functions.DoubleSolveNoPreconBiCGSTAB(out X, A, B,
            //    ConvRatioTolerance);
            //System.Diagnostics.Debug.Assert(success);
            //System.Diagnostics.Debug.WriteLine("  DoubleSolveNoPreconBiCGSTAB t = " + (System.Environment.TickCount - t));

            // Nativeを使う
            bool success = false;
            {
                System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
                int n = A.RowLength;
                int[] APtrs;
                int[] AIndexs;
                double[] AValues;
                t = System.Environment.TickCount;
                A.GetCSR(out APtrs, out AIndexs, out AValues);
                System.Diagnostics.Debug.WriteLine("  GetCSR t = " + (System.Environment.TickCount - t));
                t = System.Environment.TickCount;
                success = IvyFEM.Native.Functions.DoubleSolveNoPreconBiCGSTAB(
                    out X, n, APtrs, AIndexs, AValues, B,
                    ConvRatioTolerance);
                System.Diagnostics.Debug.Assert(success);
                System.Diagnostics.Debug.WriteLine("  DoubleSolveNoPreconBiCGSTAB t = " + (System.Environment.TickCount - t));
            }

            return success;
        }

        //////////////////////////////////////////////////////////////////////////////////////
        // complex

        private bool ComplexSolveNoPreconCOCG(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            int t;
            X = null;
            //t = System.Environment.TickCount;
            //bool isSymmetric = A.IsSymmetric();
            //System.Diagnostics.Debug.Assert(isSymmetric);
            //System.Diagnostics.Debug.WriteLine("  IsSymmetric t = " + (System.Environment.TickCount - t));
            //if (!isSymmetric)
            //{
            //    return false;
            //}

            //t = System.Environment.TickCount;
            //bool success = IvyFEM.Linear.Functions.ComplexSolveNoPreconCOCG(out X, A, B, 
            //    ConvRatioTolerance);
            //System.Diagnostics.Debug.Assert(success);
            //System.Diagnostics.Debug.WriteLine("  ComplexSolveNoPreconCOCG t = " + (System.Environment.TickCount - t));

            // Nativeを使う
            bool success = false;
            {
                System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
                int n = A.RowLength;
                int[] APtrs;
                int[] AIndexs;
                System.Numerics.Complex[] AValues;
                t = System.Environment.TickCount;
                A.GetCSR(out APtrs, out AIndexs, out AValues);
                System.Diagnostics.Debug.WriteLine("  GetCSR t = " + (System.Environment.TickCount - t));
                t = System.Environment.TickCount;
                success = IvyFEM.Native.Functions.ComplexSolveNoPreconCOCG(
                    out X, n, APtrs, AIndexs, AValues, B,
                    ConvRatioTolerance);
                System.Diagnostics.Debug.Assert(success);
                System.Diagnostics.Debug.WriteLine("  ComplexSolveNoPreconCOCG t = " + (System.Environment.TickCount - t));
            }

            return success;
        }

        private bool ComplexSolveCOCG(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            int t;
            X = null;
            //t = System.Environment.TickCount;
            //bool isSymmetric = A.IsSymmetric();
            //System.Diagnostics.Debug.Assert(isSymmetric);
            //System.Diagnostics.Debug.WriteLine("  IsSymmetric t = " + (System.Environment.TickCount - t));
            //if (!isSymmetric)
            //{
            //    return false;
            //}

            //t = System.Environment.TickCount;
            //bool success = IvyFEM.Linear.Functions.ComplexSolveCOCG(out X, A, B, ILUFillinLevel,
            //    ConvRatioTolerance);
            //System.Diagnostics.Debug.Assert(success);
            //System.Diagnostics.Debug.WriteLine("  ComplexSolveCOCG t = " + (System.Environment.TickCount - t));

            // Nativeを使う
            bool success = false;
            {
                System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
                int n = A.RowLength;
                int[] APtrs;
                int[] AIndexs;
                System.Numerics.Complex[] AValues;
                t = System.Environment.TickCount;
                A.GetCSR(out APtrs, out AIndexs, out AValues);
                System.Diagnostics.Debug.WriteLine("  GetCSR t = " + (System.Environment.TickCount - t));
                t = System.Environment.TickCount;
                success = IvyFEM.Native.Functions.ComplexSolveCOCG(
                    out X, n, APtrs, AIndexs, AValues, B, ILUFillinLevel,
                    ConvRatioTolerance);
                System.Diagnostics.Debug.Assert(success);
                System.Diagnostics.Debug.WriteLine("  ComplexSolveCOCG t = " + (System.Environment.TickCount - t));
            }

            return success;
        }

        private bool ComplexSolveICCOCG(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            int t;
            X = null;
            //t = System.Environment.TickCount;
            //bool isSymmetric = A.IsSymmetric();
            //System.Diagnostics.Debug.Assert(isSymmetric);
            //System.Diagnostics.Debug.WriteLine("  IsSymmetric t = " + (System.Environment.TickCount - t));
            //if (!isSymmetric)
            //{
            //    return false;
            //}

            //t = System.Environment.TickCount;
            //bool success = IvyFEM.Linear.Functions.ComplexSolveICCOCG(out X, A, B,
            //    ConvRatioTolerance);
            //System.Diagnostics.Debug.Assert(success);
            //System.Diagnostics.Debug.WriteLine("  ComplexSolveICCOCG t = " + (System.Environment.TickCount - t));

            // Nativeを使う
            bool success = false;
            {
                System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
                int n = A.RowLength;
                int[] APtrs;
                int[] AIndexs;
                System.Numerics.Complex[] AValues;
                t = System.Environment.TickCount;
                A.GetCSR(out APtrs, out AIndexs, out AValues);
                System.Diagnostics.Debug.WriteLine("  GetCSR t = " + (System.Environment.TickCount - t));
                t = System.Environment.TickCount;
                success = IvyFEM.Native.Functions.ComplexSolveICCOCG(
                    out X, n, APtrs, AIndexs, AValues, B,
                    ConvRatioTolerance);
                System.Diagnostics.Debug.Assert(success);
                System.Diagnostics.Debug.WriteLine("  ComplexSolveICCOCG t = " + (System.Environment.TickCount - t));
            }

            return success;
        }

        private bool ComplexSolveNoPreconBiCGSTAB(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            int t;
            X = null;

            //t = System.Environment.TickCount;
            //bool success = IvyFEM.Linear.Functions.ComplexSolveNoPreconBiCGSTAB(out X, A, B,
            //    ConvRatioTolerance);
            //System.Diagnostics.Debug.Assert(success);
            //System.Diagnostics.Debug.WriteLine("  ComplexSolveNoPreconBiCGSTAB t = " + (System.Environment.TickCount - t));

            // Nativeを使う
            bool success = false;
            {
                System.Diagnostics.Debug.Assert(A.RowLength == A.ColumnLength);
                int n = A.RowLength;
                int[] APtrs;
                int[] AIndexs;
                System.Numerics.Complex[] AValues;
                t = System.Environment.TickCount;
                A.GetCSR(out APtrs, out AIndexs, out AValues);
                System.Diagnostics.Debug.WriteLine("  GetCSR t = " + (System.Environment.TickCount - t));
                t = System.Environment.TickCount;
                success = IvyFEM.Native.Functions.ComplexSolveNoPreconBiCGSTAB(
                    out X, n, APtrs, AIndexs, AValues, B,
                    ConvRatioTolerance);
                System.Diagnostics.Debug.Assert(success);
                System.Diagnostics.Debug.WriteLine("  ComplexSolveNoPreconBiCGSTAB t = " + (System.Environment.TickCount - t));
            }

            return success;
        }
    }
}
