using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IvyFEM.Linear;

namespace IvyFEM
{
    public partial class Elastic2DFEM : Elastic2DBaseFEM
    {
        public Elastic2DFEM(FEWorld world)
        {
            World = world;
            SetupCalcABs();
        }

        protected void SetupCalcABs()
        {
            CalcElementABs.Clear();
            CalcElementABsForLine.Clear();

            // Linear/Staint Venant
            CalcElementABs.Add(CalcLinearElasticElementAB);
            CalcElementABs.Add(CalcStVenantHyperelasticElementAB);

            // Hyperelastic
            CalcElementABs.Add(CalcMooneyRivlinHyperelasticElementAB);
            CalcElementABs.Add(CalcOgdenHyperelasticElementAB);
            //CalcElementABs.Add(CalcOgdenOriginalIncompressibleHyperelasticElementAB);

            // Truss/Beam
            CalcElementABsForLine.Add(CalcTrussElementABForLine);
            CalcElementABsForLine.Add(CalcBeamElementABForLine);
            CalcElementABsForLine.Add(CalcFrameElementABForLine);
            CalcElementABsForLine.Add(CalcTimoshenkoBeamElementABForLine);
            CalcElementABsForLine.Add(CalcTimoshenkoFrameElementABForLine);
            CalcElementABsForLine.Add(CalcCorotationalFrameElementABForLine);
            CalcElementABsForLine.Add(CalcTimoshenkoCorotationalFrameElementABForLine);
            CalcElementABsForLine.Add(CalcFieldConsistentTLFrameElementABForLine);
            CalcElementABsForLine.Add(CalcTimoshenkoTLFrameElementABForLine);
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
