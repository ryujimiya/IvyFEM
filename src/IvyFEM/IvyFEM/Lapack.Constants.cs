using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Lapack
{
    public enum MatrixLayout
    {
        RowMajor = 101,
        ColMajor = 102
    }

    public class UpperLower
    {
        public static byte Upper = Convert.ToByte('U');
        public static byte Lower = Convert.ToByte('L');
    }

    public class Job
    {
        public static byte DontCompute = Convert.ToByte('N');
        public static byte Compute = Convert.ToByte('V');
    }

    public enum TransposeType
    {
        Nop,
        Transpose,
        ComplexTranspose
    }

    public class Trans
    {
        public static byte Nop = Convert.ToByte('N');
        public static byte Transpose = Convert.ToByte('T');
        public static byte ComplexTranspose = Convert.ToByte('C');

        public static byte FromTransposeType(TransposeType transposeType)
        {
            switch (transposeType)
            {
                case TransposeType.Nop:
                    return Trans.Nop;

                case TransposeType.Transpose:
                    return Trans.Transpose;

                case TransposeType.ComplexTranspose:
                    return Trans.ComplexTranspose;
            }
            return Trans.Nop;
        }
    }
}
