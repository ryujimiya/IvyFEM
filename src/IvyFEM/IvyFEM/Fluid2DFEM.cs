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
                    CalcStdGNavierStokesAB(A, B);
                    //CalcStdGNavierStokesByPicardAB(A, B); // SUPGの項がないときのTest Not converge
                    break;
                case FluidEquationType.SUPGNavierStokes:
                    //CalcSUPGNavierStokesAB(A, B); // 発散する
                    CalcSUPGNavierStokesByPicardAB(A, B);
                    break;
                //case FluidEquationType.GLSNavierStokes:
                //    CalcGLSNavierStokesByPicardAB(A, B);
                //    break;
                case FluidEquationType.StdGVorticity:
                    CalcStdGVorticityAB(A, B);
                    //CalcStdGVorticityByPicardAB(A, B); // Newton-Raphsonと同等の小さいμしか解けない
                    break;
                case FluidEquationType.SUPGVorticity:
                    //CalcSUPGVorticityAB(A, B); // 発散する
                    CalcSUPGVorticityByPicardAB(A, B);
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
                SetVorticityNeumannBC(A, B);
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

        protected override bool MustUseNewtonRaphson()
        {
            if (EquationType == FluidEquationType.Stokes)
            {
                return false;
            }
            return true;
        }
    }
}
