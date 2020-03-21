using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DTDFEM : Elastic2DBaseFEM
    {
        public double TimeStep { get; private set; } = 0;
        public double NewmarkBeta { get; private set; } = 1.0 / 4.0;
        public double NewmarkGamma { get; private set; } = 1.0 / 2.0;
        public IList<uint> UValueIds { get; set; } = new List<uint> { 0 };
        public IList<uint> PrevUValueIds { get; set; } = new List<uint> { 0 };
        public uint UValueId => UValueIds[0];
        public uint PrevUValueId => PrevUValueIds[0];
        public uint LValueId { get; private set; } = 0;

        public Elastic2DTDFEM(FEWorld world,
            double timeStep,
            double newmarkBeta, double newmarkGamma,
            uint uValueId, uint prevUValueId,
            uint lValueId = 0)
        {
            World = world;
            TimeStep = timeStep;
            NewmarkBeta = newmarkBeta;
            NewmarkGamma = newmarkGamma;
            UValueIds[0] = uValueId;
            PrevUValueIds[0] = prevUValueId;
            LValueId = lValueId;
            SetupCalcABs();
        }

        public Elastic2DTDFEM(FEWorld world,
            double timeStep,
            double newmarkBeta, double newmarkGamma,
            IList<uint> uValueIds, IList<uint> prevUValueIds,
            uint lValueId = 0)
        {
            World = world;
            TimeStep = timeStep;
            NewmarkBeta = newmarkBeta;
            NewmarkGamma = newmarkGamma;
            UValueIds = uValueIds;
            PrevUValueIds = prevUValueIds;
            LValueId = lValueId;
            SetupCalcABs();
        }

        protected void SetupCalcABs()
        {
            CalcElementABs.Clear();
            CalcElementABsForLine.Clear();

            // Linear/Saint Venant
            CalcElementABs.Add(CalcLinearElasticElementAB);
            CalcElementABs.Add(CalcSaintVenantHyperelasticElementAB);

            // Hyperelastic
            CalcElementABs.Add(CalcMooneyRivlinHyperelasticElementAB);
            CalcElementABs.Add(CalcOgdenRivlinHyperelasticElementAB);
            //CalcElementABs.Add(CalcOgdenOriginalRivlinIncompressibleHyperelasticElementAB);

            // Truss/Beam
            CalcElementABsForLine.Add(CalcTrussElementABForLine);
            CalcElementABsForLine.Add(CalcBeamElementABForLine);
            CalcElementABsForLine.Add(CalcFrameElementABForLine);
        }

        public void UpdateFieldValuesTimeDomain()
        {
            System.Diagnostics.Debug.Assert(UValueIds.Count > 0);
            System.Diagnostics.Debug.Assert(UValueIds.Count == PrevUValueIds.Count);
            for (int i = 0; i < UValueIds.Count; i++)
            {
                UpdateFieldValuesNewmarkBetaTimeDomain(
                    U, UValueIds[i], PrevUValueIds[i],
                    TimeStep,
                    NewmarkBeta, NewmarkGamma);
            }
        }

        //////////////////////////////////////////////////
        public void SetStressValue(
            uint displacementValueId, uint stressValueId, uint equivStressValueId)
        {
            Elastic2DFEMUtils.SetStressValue(
                World,
                displacementValueId, stressValueId, equivStressValueId);
        }

        public void SolvePrincipalValues(
            double[,] c, out System.Numerics.Complex[] fLambdas, out System.Numerics.Complex[][] cNormals)
        {
            Elastic2DFEMUtils.SolvePrincipalValues(
                c, out fLambdas, out cNormals);
        }
    }
}
