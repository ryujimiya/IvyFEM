using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Native
{
    public class Functions
    {
        ///////////////////////////////////////////////////////////////////////
        // double

        public static double DoubleDot(double[] X, double[] Y)
        {
            System.Diagnostics.Debug.Assert(X.Length == Y.Length);
            int n = X.Length;
            double value = IvyFEM.Native.ImportedFunctions.DoubleDot(n, X, Y);
            return value;
        }

        // matvec
        public static double[] DoubleMV(
            double alpha, int n, int AIndexsLength, int[] APtrs, int[] AIndexs, double[] AValues,
            double[] X,
            double beta, double[] Y)
        {
            double[] Z = new double[n];
            IvyFEM.Native.ImportedFunctions.DoubleMV(Z,
                alpha, n, AIndexsLength, APtrs, AIndexs, AValues, X,
                beta, Y);
            return Z;
        }

        public static double[] DoubleAxpy(double alpha, double[] X, double[] Y)
        {
            System.Diagnostics.Debug.Assert(X.Length == Y.Length);
            int n = X.Length;
            double[] Z = new double[n];
            IvyFEM.Native.ImportedFunctions.DoubleAxpy(Z, alpha, n, X, Y);
            return Z;
        }

        public static bool DoubleSolveNoPreconCG(out double[] X,
            int n, int[] APtrs, int[] AIndexs, double[] AValues, double[] B,
            double convRatioTolerance)
        {
            X = new double[n];
            bool success = IvyFEM.Native.ImportedFunctions.DoubleSolveNoPreconCG(
                X, n, AIndexs.Length, APtrs, AIndexs, AValues, B,
                convRatioTolerance);
            return success;
        }

        public static void DoubleCalcILU(
            out int[] LUPtrs, out int[] LUIndexs, out double[] LUValues,
            int n, int[] APtrs, int[] AIndexs, double[] AValues,
            int fillinLevel)
        {
            unsafe
            {
                int workLUIndexsLength = 0;
                int* workLUPtrs = null;
                int* workLUIndexs = null;
                double* workLUValues = null;
                IvyFEM.Native.ImportedFunctions.DoubleCalcILU(
                    &workLUIndexsLength, &workLUPtrs, &workLUIndexs, &workLUValues,
                    n, AIndexs.Length, APtrs, AIndexs, AValues,                    
                    fillinLevel);

                LUPtrs = new int[n + 1];
                LUIndexs = new int[workLUIndexsLength];
                LUValues = new double[workLUIndexsLength];
                for (int i = 0; i < LUPtrs.Length; i++)
                {
                    LUPtrs[i] = workLUPtrs[i];
                }
                for (int i = 0; i < LUIndexs.Length; i++)
                {
                    LUIndexs[i] = workLUIndexs[i];
                    LUValues[i] = workLUValues[i];
                }
                IvyFEM.Native.ImportedFunctions.DoubleDeleteCSR(workLUPtrs, workLUIndexs, workLUValues);
            }
        }

        public static void DoubleSolveLU(
            out double[] X, int n, int[] LUPtrs, int[] LUIndexs, double[] LUValues, double[] B)
        {
            int LUIndexsLength = LUIndexs.Length;
            X = new double[n];
            IvyFEM.Native.ImportedFunctions.DoubleSolveLU(
                X, n, LUIndexsLength, LUPtrs, LUIndexs, LUValues, B);
        }

        public static bool DoubleSolvePreconditionedCG(out double[] X,
            int n, int[] APtrs, int[] AIndexs, double[] AValues, double[] B,
            int[] LUPtrs, int[] LUIndexs, double[] LUValues,
            double convRatioTolerance)
        {
            X = new double[n];
            bool success = IvyFEM.Native.ImportedFunctions.DoubleSolvePreconditionedCG(
                X,
                n, AIndexs.Length, APtrs, AIndexs, AValues, B,
                LUIndexs.Length, LUPtrs, LUIndexs, LUValues,
                convRatioTolerance);
            return success;
        }

        public static bool DoubleSolveCG(out double[] X,
            int n, int[] APtrs, int[] AIndexs, double[] AValues, double[] B, int fillinLevel,
            double convRatioTolerance)
        {
            X = new double[n];
            bool success = IvyFEM.Native.ImportedFunctions.DoubleSolveCG(
                X, n, AIndexs.Length, APtrs, AIndexs, AValues, B, fillinLevel,
                convRatioTolerance);
            return success;
        }

        public static void DoubleCalcILUWithPivoting(
            out int[] LUPtrs, out int[] LUIndexs, out double[] LUValues,
            out int[] pivot,
            int n, int[] APtrs, int[] AIndexs, double[] AValues,
            int fillinLevel)
        {
            pivot = new int[n];
            unsafe
            {
                int workLUIndexsLength = 0;
                int* workLUPtrs = null;
                int* workLUIndexs = null;
                double* workLUValues = null;
                IvyFEM.Native.ImportedFunctions.DoubleCalcILUWithPivoting(
                    &workLUIndexsLength, &workLUPtrs, &workLUIndexs, &workLUValues,
                    pivot,
                    n, AIndexs.Length, APtrs, AIndexs, AValues,
                    fillinLevel);

                LUPtrs = new int[n + 1];
                LUIndexs = new int[workLUIndexsLength];
                LUValues = new double[workLUIndexsLength];
                for (int i = 0; i < LUPtrs.Length; i++)
                {
                    LUPtrs[i] = workLUPtrs[i];
                }
                for (int i = 0; i < LUIndexs.Length; i++)
                {
                    LUIndexs[i] = workLUIndexs[i];
                    LUValues[i] = workLUValues[i];
                }
                IvyFEM.Native.ImportedFunctions.DoubleDeleteCSR(workLUPtrs, workLUIndexs, workLUValues);
            }
        }

        public static bool DoubleSolveCGWithPivoting(out double[] X,
            int n, int[] APtrs, int[] AIndexs, double[] AValues, double[] B, int fillinLevel,
            double convRatioTolerance)
        {
            X = new double[n];
            bool success = IvyFEM.Native.ImportedFunctions.DoubleSolveCGWithPivoting(
                X, n, AIndexs.Length, APtrs, AIndexs, AValues, B, fillinLevel,
                convRatioTolerance);
            return success;
        }

        public static void DoubleCalcIC(
            out int[] LUPtrs, out int[] LUIndexs, out double[] LUValues,
            int n, int[] APtrs, int[] AIndexs, double[] AValues)
        {
            unsafe
            {
                int workLUIndexsLength = 0;
                int* workLUPtrs = null;
                int* workLUIndexs = null;
                double* workLUValues = null;
                IvyFEM.Native.ImportedFunctions.DoubleCalcIC(
                    &workLUIndexsLength, &workLUPtrs, &workLUIndexs, &workLUValues,
                    n, AIndexs.Length, APtrs, AIndexs, AValues);

                LUPtrs = new int[n + 1];
                LUIndexs = new int[workLUIndexsLength];
                LUValues = new double[workLUIndexsLength];
                for (int i = 0; i < LUPtrs.Length; i++)
                {
                    LUPtrs[i] = workLUPtrs[i];
                }
                for (int i = 0; i < LUIndexs.Length; i++)
                {
                    LUIndexs[i] = workLUIndexs[i];
                    LUValues[i] = workLUValues[i];
                }
                IvyFEM.Native.ImportedFunctions.DoubleDeleteCSR(workLUPtrs, workLUIndexs, workLUValues);
            }
        }

        public static bool DoubleSolveICCG(out double[] X,
            int n, int[] APtrs, int[] AIndexs, double[] AValues, double[] B,
            double convRatioTolerance)
        {
            X = new double[n];
            bool success = IvyFEM.Native.ImportedFunctions.DoubleSolveICCG(
                X, n, AIndexs.Length, APtrs, AIndexs, AValues, B,
                convRatioTolerance);
            return success;
        }

        public static bool DoubleSolveNoPreconBiCGSTAB(out double[] X,
            int n, int[] APtrs, int[] AIndexs, double[] AValues, double[] B,
            double convRatioTolerance)
        {
            X = new double[n];
            bool success = IvyFEM.Native.ImportedFunctions.DoubleSolveNoPreconBiCGSTAB(
                X, n, AIndexs.Length, APtrs, AIndexs, AValues, B,
                convRatioTolerance);
            return success;
        }

        ///////////////////////////////////////////////////////////////////////
        // complex

        // X^H * Y
        public static System.Numerics.Complex ComplexDotc(System.Numerics.Complex[] X, System.Numerics.Complex[] Y)
        {
            System.Diagnostics.Debug.Assert(X.Length == Y.Length);
            int n = X.Length;
            System.Numerics.Complex value = 0;
            unsafe
            {
                fixed (System.Numerics.Complex* XP = &X[0])
                fixed (System.Numerics.Complex* YP = &Y[0])
                {
                    IvyFEM.Native.ImportedFunctions.ComplexDotc(&value, n, XP, YP);
                }
            }
            return value;
        }

        // X^T * Y
        public static System.Numerics.Complex ComplexDotu(System.Numerics.Complex[] X, System.Numerics.Complex[] Y)
        {
            System.Diagnostics.Debug.Assert(X.Length == Y.Length);
            int n = X.Length;
            System.Numerics.Complex value = 0;
            unsafe
            {
                fixed (System.Numerics.Complex* XP = &X[0])
                fixed (System.Numerics.Complex* YP = &Y[0])
                {
                    IvyFEM.Native.ImportedFunctions.ComplexDotu(&value, n, XP, YP);
                }
            }
            return value;
        }

        // matvec
        public static System.Numerics.Complex[] ComplexMV(
            System.Numerics.Complex alpha, 
            int n, int AIndexsLength, int[] APtrs, int[] AIndexs, System.Numerics.Complex[] AValues,
            System.Numerics.Complex[] X,
            System.Numerics.Complex beta, System.Numerics.Complex[] Y)
        {
            System.Numerics.Complex[] Z = new System.Numerics.Complex[n];
            unsafe
            {
                fixed (System.Numerics.Complex* ZP = &Z[0])
                fixed (System.Numerics.Complex* AValuesP = &AValues[0])
                fixed (System.Numerics.Complex* XP = &X[0])
                fixed (System.Numerics.Complex* YP = &Y[0])
                {
                    IvyFEM.Native.ImportedFunctions.ComplexMV(
                        ZP,
                        alpha, n, AIndexsLength, APtrs, AIndexs, AValuesP,
                        XP,
                        beta, YP);
                }

            }
            return Z;
        }

        public static System.Numerics.Complex[] ComplexAxpy(
            System.Numerics.Complex alpha, System.Numerics.Complex[] X, System.Numerics.Complex[] Y)
        {
            System.Diagnostics.Debug.Assert(X.Length == Y.Length);
            int n = X.Length;
            System.Numerics.Complex[] Z = new System.Numerics.Complex[n];
            unsafe
            {
                fixed (System.Numerics.Complex* ZP = &Z[0])
                fixed (System.Numerics.Complex* XP = &X[0])
                fixed (System.Numerics.Complex* YP = &Y[0])
                {
                    IvyFEM.Native.ImportedFunctions.ComplexAxpy(ZP, alpha, n, XP, YP);
                }
            }
            return Z;
        }

        public static bool ComplexSolveNoPreconCOCG(out System.Numerics.Complex[] X,
            int n, int[] APtrs, int[] AIndexs, System.Numerics.Complex[] AValues,
            System.Numerics.Complex[] B,
            double convRatioTolerance)
        {
            X = new System.Numerics.Complex[n];
            bool success = false;
            unsafe
            {
                fixed (System.Numerics.Complex* XP = &X[0])
                fixed (System.Numerics.Complex* AValuesP = &AValues[0])
                fixed (System.Numerics.Complex* BP = &B[0])
                {
                    success = IvyFEM.Native.ImportedFunctions.ComplexSolveNoPreconCOCG(
                        XP, n, AIndexs.Length, APtrs, AIndexs, AValuesP, BP,
                        convRatioTolerance);
                }
            }
            return success;
        }

        public static void ComplexCalcILU(
            out int[] LUPtrs, out int[] LUIndexs, out System.Numerics.Complex[] LUValues,
            int n, int[] APtrs, int[] AIndexs, System.Numerics.Complex[] AValues,
            int fillinLevel)
        {
            unsafe
            {
                int workLUIndexsLength = 0;
                int* workLUPtrs = null;
                int* workLUIndexs = null;
                System.Numerics.Complex* workLUValues = null;
                fixed (System.Numerics.Complex* AValuesP = &AValues[0])
                {
                    IvyFEM.Native.ImportedFunctions.ComplexCalcILU(
                        &workLUIndexsLength, &workLUPtrs, &workLUIndexs, &workLUValues,
                        n, AIndexs.Length, APtrs, AIndexs, AValuesP,
                        fillinLevel);
                }

                LUPtrs = new int[n + 1];
                LUIndexs = new int[workLUIndexsLength];
                LUValues = new System.Numerics.Complex[workLUIndexsLength];
                for (int i = 0; i < LUPtrs.Length; i++)
                {
                    LUPtrs[i] = workLUPtrs[i];
                }
                for (int i = 0; i < LUIndexs.Length; i++)
                {
                    LUIndexs[i] = workLUIndexs[i];
                    LUValues[i] = workLUValues[i];
                }
                IvyFEM.Native.ImportedFunctions.ComplexDeleteCSR(workLUPtrs, workLUIndexs, workLUValues);
            }
        }

        public static void ComplexSolveLU(
            out System.Numerics.Complex[] X,
            int n, int[] LUPtrs, int[] LUIndexs, System.Numerics.Complex[] LUValues,
            System.Numerics.Complex[] B)
        {
            int LUIndexsLength = LUIndexs.Length;
            X = new System.Numerics.Complex[n];
            unsafe
            {
                fixed (System.Numerics.Complex* XP = &X[0])
                fixed (System.Numerics.Complex* LUValuesP = &LUValues[0])
                fixed (System.Numerics.Complex* BP = &B[0])
                {
                    IvyFEM.Native.ImportedFunctions.ComplexSolveLU(
                        XP, n, LUIndexsLength, LUPtrs, LUIndexs, LUValuesP, BP);
                }
            }
        }

        public static bool ComplexSolvePreconditionedCOCG(
            out System.Numerics.Complex[] X,
            int n, int[] APtrs, int[] AIndexs, System.Numerics.Complex[] AValues,
            System.Numerics.Complex[] B,
            int[] LUPtrs, int[] LUIndexs, System.Numerics.Complex[] LUValues,
            double convRatioTolerance)
        {
            X = new System.Numerics.Complex[n];
            bool success = false;
            unsafe
            {
                fixed (System.Numerics.Complex* XP = &X[0])
                fixed (System.Numerics.Complex* AValuesP = &AValues[0])
                fixed (System.Numerics.Complex* BP = &B[0])
                fixed (System.Numerics.Complex* LUValuesP = &LUValues[0])
                {
                    success = IvyFEM.Native.ImportedFunctions.ComplexSolvePreconditionedCOCG(
                        XP,
                        n, AIndexs.Length, APtrs, AIndexs, AValuesP,
                        BP,
                        LUIndexs.Length, LUPtrs, LUIndexs, LUValuesP,
                        convRatioTolerance);
                }
            }
            return success;
        }

        public static bool ComplexSolveCOCG(out System.Numerics.Complex[] X,
            int n, int[] APtrs, int[] AIndexs, System.Numerics.Complex[] AValues,
            System.Numerics.Complex[] B, int fillinLevel,
            double convRatioTolerance)
        {
            X = new System.Numerics.Complex[n];
            bool success = false;
            unsafe
            {
                fixed (System.Numerics.Complex* XP = &X[0])
                fixed (System.Numerics.Complex* AValuesP = &AValues[0])
                fixed (System.Numerics.Complex* BP = &B[0])
                {
                    success = IvyFEM.Native.ImportedFunctions.ComplexSolveCOCG(
                        XP, n, AIndexs.Length, APtrs, AIndexs, AValuesP, BP, fillinLevel,
                        convRatioTolerance);
                }
            }
            return success;
        }

        public static void ComplexCalcIC(
            out int[] LUPtrs, out int[] LUIndexs, out System.Numerics.Complex[] LUValues,
            int n, int[] APtrs, int[] AIndexs, System.Numerics.Complex[] AValues)
        {
            unsafe
            {
                int workLUIndexsLength = 0;
                int* workLUPtrs = null;
                int* workLUIndexs = null;
                System.Numerics.Complex* workLUValues = null;
                fixed (System.Numerics.Complex* AValuesP = &AValues[0])
                {
                    IvyFEM.Native.ImportedFunctions.ComplexCalcIC(
                        &workLUIndexsLength, &workLUPtrs, &workLUIndexs, &workLUValues,
                        n, AIndexs.Length, APtrs, AIndexs, AValuesP);
                }

                LUPtrs = new int[n + 1];
                LUIndexs = new int[workLUIndexsLength];
                LUValues = new System.Numerics.Complex[workLUIndexsLength];
                for (int i = 0; i < LUPtrs.Length; i++)
                {
                    LUPtrs[i] = workLUPtrs[i];
                }
                for (int i = 0; i < LUIndexs.Length; i++)
                {
                    LUIndexs[i] = workLUIndexs[i];
                    LUValues[i] = workLUValues[i];
                }
                IvyFEM.Native.ImportedFunctions.ComplexDeleteCSR(workLUPtrs, workLUIndexs, workLUValues);
            }
        }

        public static bool ComplexSolveICCOCG(out System.Numerics.Complex[] X,
            int n, int[] APtrs, int[] AIndexs, System.Numerics.Complex[] AValues,
            System.Numerics.Complex[] B,
            double convRatioTolerance)
        {
            X = new System.Numerics.Complex[n];
            bool success = false;
            unsafe
            {
                fixed (System.Numerics.Complex* XP = &X[0])
                fixed (System.Numerics.Complex* AValuesP = &AValues[0])
                fixed (System.Numerics.Complex* BP = &B[0])
                {
                    success = IvyFEM.Native.ImportedFunctions.ComplexSolveICCOCG(
                        XP, n, AIndexs.Length, APtrs, AIndexs, AValuesP, BP,
                        convRatioTolerance);
                }
            }
            return success;
        }

        public static bool ComplexSolveNoPreconBiCGSTAB(out System.Numerics.Complex[] X,
            int n, int[] APtrs, int[] AIndexs, System.Numerics.Complex[] AValues,
            System.Numerics.Complex[] B,
            double convRatioTolerance)
        {
            X = new System.Numerics.Complex[n];
            bool success = false;
            unsafe
            {
                fixed (System.Numerics.Complex* XP = &X[0])
                fixed (System.Numerics.Complex* AValuesP = &AValues[0])
                fixed (System.Numerics.Complex* BP = &B[0])
                {
                    success = IvyFEM.Native.ImportedFunctions.ComplexSolveNoPreconBiCGSTAB(
                        XP, n, AIndexs.Length, APtrs, AIndexs, AValuesP, BP,
                        convRatioTolerance);
                }
            }
            return success;
        }
    }
}
