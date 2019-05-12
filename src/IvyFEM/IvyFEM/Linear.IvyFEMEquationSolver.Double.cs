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

        private bool DoubleSolveBiCGSTAB(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            int t;
            X = null;

            //t = System.Environment.TickCount;
            //bool success = IvyFEM.Linear.Functions.DoubleSolveBiCGSTAB(out X, A, B, ILUFillinLevel,
            //    ConvRatioTolerance);
            //System.Diagnostics.Debug.Assert(success);
            //System.Diagnostics.Debug.WriteLine("  DoubleSolveBiCGSTAB t = " + (System.Environment.TickCount - t));

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
                success = IvyFEM.Native.Functions.DoubleSolveBiCGSTAB(
                    out X, n, APtrs, AIndexs, AValues, B, ILUFillinLevel,
                    ConvRatioTolerance);
                System.Diagnostics.Debug.Assert(success);
                System.Diagnostics.Debug.WriteLine("  DoubleSolveBiCGSTAB t = " + (System.Environment.TickCount - t));
            }

            return success;
        }

        private bool DoubleSolveBiCGSTABWithPivoting(out double[] X, DoubleSparseMatrix A, double[] B)
        {
            int t;
            X = null;

            //t = System.Environment.TickCount;
            //bool success = IvyFEM.Linear.Functions.DoubleSolveBiCGSTABWithPivoting(out X, A, B, ILUFillinLevel,
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
                success = IvyFEM.Native.Functions.DoubleSolveBiCGSTABWithPivoting(
                    out X, n, APtrs, AIndexs, AValues, B, ILUFillinLevel,
                    ConvRatioTolerance);
                System.Diagnostics.Debug.Assert(success);
                System.Diagnostics.Debug.WriteLine("  DoubleSolveCGWithPivoting t = " + (System.Environment.TickCount - t));
            }

            return success;
        }
    }
}
