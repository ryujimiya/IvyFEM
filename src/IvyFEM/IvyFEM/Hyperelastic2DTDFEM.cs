using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Hyperelastic2DTDFEM : Hyperelastic2DDerivedBaseFEM
    {
        public double TimeStep { get; private set; } = 0;
        public double NewmarkBeta { get; private set; } = 1.0 / 4.0;
        public double NewmarkGamma { get; private set; } = 1.0 / 2.0;
        public uint UValueId { get; set; } = 0;
        public uint PrevUValueId { get; set; } = 0;
        public uint LValueId { get; set; } = 0;

        public Hyperelastic2DTDFEM(FEWorld world,
            double timeStep,
            double newmarkBeta, double newmarkGamma,
            uint uValueId, uint prevUValueId,
            uint lValueId)
        {
            World = world;
            TimeStep = timeStep;
            NewmarkBeta = newmarkBeta;
            NewmarkGamma = newmarkGamma;
            UValueId = uValueId;
            PrevUValueId = prevUValueId;
            LValueId = lValueId;
            SetupCalcABs();
        }

        protected void SetupCalcABs()
        {
            CalcElementABs.Clear();
            CalcElementABs.Add(CalcMooneyRivlinHyperelasticElementAB);
            CalcElementABs.Add(CalcOgdenRivlinHyperelasticElementAB);
            //CalcElementABs.Add(CalcOgdenOriginalRivlinIncompressibleHyperelasticElementAB);
        }

        public void UpdateFieldValuesTimeDomain()
        {
            UpdateFieldValuesTimeDomain(
                U, UValueId, PrevUValueId,
                TimeStep,
                NewmarkBeta, NewmarkGamma);
        }
    }
}
