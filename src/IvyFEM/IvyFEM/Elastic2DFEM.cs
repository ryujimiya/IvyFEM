using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IvyFEM.Linear;

namespace IvyFEM
{
    public partial class Elastic2DFEM : Elastic2DDerivedBaseFEM
    {
        public Elastic2DFEM(FEWorld world)
        {
            World = world;
            SetupCalcABs();
        }

        protected void SetupCalcABs()
        {
            CalcElementABs.Clear();
            CalcElementABs.Add(CalcLinearElasticElementAB);
            CalcElementABs.Add(CalcSaintVenantHyperelasticElementAB);
        }
    }
}
