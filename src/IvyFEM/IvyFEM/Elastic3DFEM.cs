using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IvyFEM.Linear;

namespace IvyFEM
{
    public partial class Elastic3DFEM : Elastic3DBaseFEM
    {
        public Elastic3DFEM(FEWorld world)
        {
            World = world;
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

            InitNodeValuessForPlate.Add(InitMITCStVenantThicknessStretchPlateNodeValues);
            InitElementValuessForPlate.Add(InitMITCStVenantThicknessStretchPlateElementValues);

            InitNodeValuessForPlate.Add(InitMITCMooneyRivlinPlateNodeValues);
            InitElementValuessForPlate.Add(InitMITCMooneyRivlinPlateElementValues);

            // Update
            // Plate
            UpdateNodeValuessForPlate.Add(UpdateMITCStVenantPlateNodeValues);
            UpdateElementValuessForPlate.Add(UpdateMITCStVenantPlateElementValues);

            UpdateNodeValuessForPlate.Add(UpdateMITCStVenantThicknessStretchPlateNodeValues);
            UpdateElementValuessForPlate.Add(UpdateMITCStVenantThicknessStretchPlateElementValues);

            UpdateNodeValuessForPlate.Add(UpdateMITCMooneyRivlinPlateNodeValues);
            UpdateElementValuessForPlate.Add(UpdateMITCMooneyRivlinPlateElementValues);

            //--------------------------------
            CalcElementABs.Clear();
            CalcElementABsForPlate.Clear();

            // Linear/Staint Venant
            CalcElementABs.Add(CalcLinearElasticElementAB);
            CalcElementABs.Add(CalcStVenantHyperelasticElementAB);

            // Hyperelastic
            CalcElementABs.Add(CalcMooneyRivlinHyperelasticElementAB);
            CalcElementABs.Add(CalcOgdenHyperelasticElementAB);

            // Plate
            CalcElementABsForPlate.Add(CalcDKTPlateElementAB);
            CalcElementABsForPlate.Add(CalcMindlinPlateElementAB);
            CalcElementABsForPlate.Add(CalcMITCLinearPlateElementAB);
            CalcElementABsForPlate.Add(CalcMITCStVenantPlateElementAB);
            CalcElementABsForPlate.Add(CalcMITCStVenantThicknessStretchPlateElementAB);
            CalcElementABsForPlate.Add(CalcMITCMooneyRivlinPlateElementAB);
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
