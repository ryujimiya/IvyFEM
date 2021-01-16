using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic3DTDFEM : Elastic3DBaseFEM
    {
        public double TimeStep { get; private set; } = 0;
        public double NewmarkBeta { get; private set; } = 1.0 / 4.0;
        public double NewmarkGamma { get; private set; } = 1.0 / 2.0;
        public IList<uint> UValueIds { get; set; } = new List<uint> { 0 };
        public IList<uint> PrevUValueIds { get; set; } = new List<uint> { 0 };
        public uint UValueId => UValueIds[0];
        public uint PrevUValueId => PrevUValueIds[0];
        public uint LValueId { get; private set; } = 0;

        public Elastic3DTDFEM(FEWorld world,
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

        public Elastic3DTDFEM(FEWorld world,
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
            //--------------------------------
            InitNodeValuessForPlate.Clear();
            InitElementValuessForPlate.Clear();
            UpdateNodeValuessForPlate.Clear();
            UpdateElementValuessForPlate.Clear();

            // Init
            // Plate
            InitNodeValuessForPlate.Add(InitMITCStVenantPlateNodeValues);
            InitElementValuessForPlate.Add(InitMITCStVenantPlateElementValues);

            // Update
            // Plate
            UpdateNodeValuessForPlate.Add(UpdateMITCStVenantPlateNodeValues);
            UpdateElementValuessForPlate.Add(UpdateMITCStVenantPlateElementValues);

            //--------------------------------
            CalcElementABs.Clear();
            CalcElementABsForPlate.Clear();

            // Linear/Saint Venant
            CalcElementABs.Add(CalcLinearElasticElementAB);
            CalcElementABs.Add(CalcStVenantHyperelasticElementAB);

            // Hyperelastic
            CalcElementABs.Add(CalcMooneyRivlinHyperelasticElementAB);
            CalcElementABs.Add(CalcOgdenRivlinHyperelasticElementAB);

            // Plate
            CalcElementABsForPlate.Add(CalcDKTPlateElementAB);
            CalcElementABsForPlate.Add(CalcMindlinPlateElementAB);
            CalcElementABsForPlate.Add(CalcMITCLinearPlateElementAB);
            CalcElementABsForPlate.Add(CalcMITCStVenantPlateElementAB);
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

        /////////////////////////////////////////////////////////
        public void SolvePrincipalValues(
            double[,] c, out System.Numerics.Complex[] fLambdas, out System.Numerics.Complex[][] cNormals)
        {
            Elastic3DFEMUtils.SolvePrincipalValues(
                c, out fLambdas, out cNormals);
        }
    }
}
