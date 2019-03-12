using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Lis
{
    public class Constants
    {
        public const int LisCommWorld = 0x1;
        public const int LisMatrixOptionLen = 10;

        public const int LisOptionsLen = 27;
        public const int LisParamsLen = 15;

        public const int LisMatrixDecidingSize = -((int)MatrixType.LisMatrixRCO + 1);
        public const int LisMatrixNull = -((int)MatrixType.LisMatrixRCO + 2);
        public const MatrixType LisMatrixDefault = MatrixType.LisMatrixCSR;
        public const MatrixType LisMatrixPoint = MatrixType.LisMatrixCSR;
        public const MatrixType LisMatrixBlock = MatrixType.LisMatrixBSR;
    }

    public enum SetValueFlag
    {
        LisInsValue = 0,
        LisAddValue = 1,
        LisSubValue = 2
    }

    public enum MatrixType
    {
        LisMatrixAssembling = 0,
        LisMatrixCSR = 1,
        LisMatrixCSC = 2,
        LisMatrixMSR = 3,
        LisMatrixDIA = 4,
        LisMatrixCDS = 4,
        LisMatrixELL = 5,
        LisMatrixJAD = 6,
        LisMatrixBSR = 7,
        LisMatrixBSC = 8,
        LisMatrixVBR = 9,
        LisMatrixCOO = 10,
        LisMatrixDENSE = 11,
        LisMatrixDNS = 11,
        LisMatrixRCO = 255,
        LisMatrixTJAD = 12,
        LisMatrixBJAD = 13,
        LisMatrixBCR = 14,
        LisMatrixCJAD = 15,
        LisMatrixPCSR = 16,
        LisMatrixLCSR = 17,
        LisMatrixLJAD = 18,
        LisMatrixLBSR = 19,
        LisMatrixCDIA = 20,
        LisMatrixMSC = 21
    }
}
