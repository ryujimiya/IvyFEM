using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Fluid2DFEM : Fluid2DBaseFEM
    {
        public Fluid2DFEM(FEWorld world)
        {
            World = world;
        }

        protected override void CalcAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            switch (EquationType)
            {
                case FluidEquationType.Stokes:
                    CalcStokesAB(A, B);
                    break;
                case FluidEquationType.StdGNavierStokes:
                    //CalcStdGNavierStokesAB(A, B);
                    CalcStdGNavierStokesByNewtonAB(A, B);
                    break;
                case FluidEquationType.SUPGNavierStokes:
                    CalcSUPGNavierStokesAB(A, B);
                    break;
                case FluidEquationType.StdGVorticity:
                    //CalcStdGVorticityAB(A, B);
                    CalcStdGVorticityByNewtonAB(A, B);
                    break;
                case FluidEquationType.SUPGVorticity:
                    CalcSUPGVorticityAB(A, B);
                    break;
                case FluidEquationType.StdGPressurePoisson:
                    CalcStdGPressurePoissonAB(A, B);
                    break;
                default:
                    System.Diagnostics.Debug.Assert(false);
                    break;
            }
        }

        protected override void SetSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            if (EquationType == FluidEquationType.StdGVorticity ||
                EquationType == FluidEquationType.SUPGVorticity)
            {
                SetVorticitySpecialBC(A, B);
            }
            else if (EquationType == FluidEquationType.StdGPressurePoisson)
            {
                SetPressurePoissonSpecialBC(A, B, 0.0, 0.0, 0.0, 0);
            }
        }

        protected override void PostSolve()
        {
            if (EquationType == FluidEquationType.StdGVorticity ||
                EquationType == FluidEquationType.SUPGVorticity)
            {
                VorticityPostSolve();
            }
        }

        protected override bool MustUseNonlinearIter()
        {
            if (EquationType == FluidEquationType.Stokes)
            {
                return false;
            }
            return true;
        }
    }
}
