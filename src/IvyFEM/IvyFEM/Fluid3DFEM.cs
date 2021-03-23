using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Fluid3DFEM : Fluid3DBaseFEM
    {
        public Fluid3DFEM(FEWorld world)
        {
            World = world;
        }

        protected override void CalcAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            switch (EquationType)
            {
                case FluidEquationType.StdGNavierStokes:
                    //CalcStdGNavierStokesAB(A, B);
                    //CalcStdGNavierStokesByNewtonAB(A, B); // 2Dではこちらが有効になっている
                    CalcStdGNavierStokesByNewtonAB(A, B);
                    break;
                case FluidEquationType.SUPGNavierStokes:
                    CalcSUPGNavierStokesAB(A, B);
                    break;
                default:
                    System.Diagnostics.Debug.Assert(false);
                    break;
            }
        }

        protected override void SetSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {

        }

        protected override void PostSolve()
        {

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
