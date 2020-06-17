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
            InitNodeValuess.Clear();
            InitElementValuess.Clear();
            UpdateNodeValuess.Clear();
            UpdateElementValuess.Clear();

            // Init
            // Plate
            InitNodeValuess.Add(InitMITCStVenantPlateNodeValues);
            InitElementValuess.Add(InitMITCStVenantPlateElementValues);

            // Update
            // Plate
            UpdateNodeValuess.Add(UpdateMITCStVenantPlateNodeValues);
            UpdateElementValuess.Add(UpdateMITCStVenantPlateElementValues);

            //--------------------------------
            CalcElementABs.Clear();

            // Plate
            CalcElementABs.Add(CalcDKTPlateElementAB);
            CalcElementABs.Add(CalcMindlinPlateElementAB);
            CalcElementABs.Add(CalcMITCLinearPlateElementAB);
            CalcElementABs.Add(CalcMITCStVenantPlateElementAB);
        }
    }
}
