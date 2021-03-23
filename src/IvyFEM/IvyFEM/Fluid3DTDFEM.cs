using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Fluid3DTDFEM : Fluid2DBaseFEM
    {
        public double TimeStep { get; private set; } = 0;
        public double NewmarkBeta { get; private set; } = 1.0 / 4.0;
        public double NewmarkGamma { get; private set; } = 1.0 / 2.0;
        public uint ValueId { get; private set; } = 0;
        public uint PrevValueId { get; private set; } = 0;

        public Fluid3DTDFEM(FEWorld world,
            double timeStep,
            double newmarkBeta, double newmarkGamma,
            uint valueId, uint prevValueId)
        {
            World = world;
            TimeStep = timeStep;
            NewmarkBeta = newmarkBeta;
            NewmarkGamma = newmarkGamma;
            ValueId = valueId;
            PrevValueId = prevValueId;
        }

        public void UpdateFieldValuesTimeDomain()
        {
            UpdateFieldValuesNewmarkBetaTimeDomain(
                U, ValueId, PrevValueId,
                TimeStep,
                NewmarkBeta, NewmarkGamma);
        }

        protected override void CalcAB(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            switch (EquationType)
            {
                case FluidEquationType.StdGNavierStokes:
                    //CalcStdGNavierStokesAB(A, B);
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
