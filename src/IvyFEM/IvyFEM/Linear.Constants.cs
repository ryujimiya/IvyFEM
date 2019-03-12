using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Linear
{
    public class Constants
    {
        public const double ConvRatioTolerance = 1.0e-14;
        public const int MaxIter = 1000;
    }

    public enum LapackEquationSolverMethod
    {
        Default,
        Dense,
        Band,
        PositiveDefiniteBand
    }

    public enum LisEquationSolverMethod
    {
        Default,
        CG,
        BiCG,
        CGS,
        BiCGSTAB,
        BiCGSTABl,
        GPBiCG,
        TFQMR,
        Orthominm,
        GMRESm,
        Jacobi,
        GaussSeidel,
        SOR,
        BiCGSafe,
        CR,
        BiCR,
        CRS,
        BiCRSTAB,
        GPBiCR,
        BiCRSafe,
        FGMRESm,
        IDRs,
        IDR1,
        MINRES,
        COCG,
        COCR
    }

    public enum IvyFEMEquationSolverMethod
    {
        Default,
        NoPreconCG,
        CG,
        CGWithPivoting,
        ICCG,
        NoPreconCOCG,
        COCG,
        ICCOCG,
        NoPreconBiCGSTAB
    }
}
