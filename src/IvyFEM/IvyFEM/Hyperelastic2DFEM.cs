using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Hyperelastic2DFEM : Hyperelastic2DDerivedBaseFEM
    {
        public Hyperelastic2DFEM(FEWorld world)
        {
            World = world;
            SetupCalcABs();
        }

        protected void SetupCalcABs()
        {
            CalcElementABs.Clear();
            CalcElementABs.Add(CalcMooneyRivlinHyperelasticElementAB);
            CalcElementABs.Add(CalcOgdenHyperelasticElementAB);
            //CalcElementABs.Add(CalcOgdenOriginalIncompressibleHyperelasticElementAB);
        }
    }
}
