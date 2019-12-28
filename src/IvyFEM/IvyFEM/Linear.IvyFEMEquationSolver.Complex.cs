using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public partial class IvyFEMEquationSolver
    {
        //////////////////////////////////////////////////////////////////////////////////////
        // complex

        private bool ComplexSolveNoPreconCOCG(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            int t;
            X = null;
            //t = System.Environment.TickCount;
            //double th = 1.0e-12;
            //bool isSymmetric = A.AssertSymmetric(th);
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
            //double th = 1.0e-12;
            //bool isSymmetric = A.AssertSymmetric(th);
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
            //double th = 1.0e-12;
            //bool isSymmetric = A.AssertSymmetric(th);
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

        private bool ComplexSolveBiCGSTAB(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            int t;
            X = null;

            //t = System.Environment.TickCount;
            //bool success = IvyFEM.Linear.Functions.ComplexSolveBiCGSTAB(out X, A, B, ILUFillinLevel,
            //    ConvRatioTolerance);
            //System.Diagnostics.Debug.Assert(success);
            //System.Diagnostics.Debug.WriteLine("  ComplexSolveBiCGSTAB t = " + (System.Environment.TickCount - t));

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
                success = IvyFEM.Native.Functions.ComplexSolveBiCGSTAB(
                    out X, n, APtrs, AIndexs, AValues, B, ILUFillinLevel,
                    ConvRatioTolerance);
                System.Diagnostics.Debug.Assert(success);
                System.Diagnostics.Debug.WriteLine("  ComplexSolveBiCGSTAB t = " + (System.Environment.TickCount - t));
            }

            return success;
        }

        private bool ComplexSolveBiCGSTABWithPivoting(
            out System.Numerics.Complex[] X, ComplexSparseMatrix A, System.Numerics.Complex[] B)
        {
            int t;
            X = null;

            //t = System.Environment.TickCount;
            //bool success = IvyFEM.Linear.Functions.ComplexSolveBiCGSTABWithPivoting(out X, A, B, ILUFillinLevel,
            //    ConvRatioTolerance);
            //System.Diagnostics.Debug.Assert(success);
            //System.Diagnostics.Debug.WriteLine("  ComplexSolveCGWithPivoting t = " + (System.Environment.TickCount - t));

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
                success = IvyFEM.Native.Functions.ComplexSolveBiCGSTABWithPivoting(
                    out X, n, APtrs, AIndexs, AValues, B, ILUFillinLevel,
                    ConvRatioTolerance);
                System.Diagnostics.Debug.Assert(success);
                System.Diagnostics.Debug.WriteLine("  ComplexSolveCGWithPivoting t = " + (System.Environment.TickCount - t));
            }

            return success;
        }
    }
}
