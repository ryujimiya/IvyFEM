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
            CalcElementABs.Clear();

            // Plate
            CalcElementABs.Add(CalcDKTPlateElementAB);
        }
    }
}
