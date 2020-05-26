using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MITCLinearPlateMaterial : PlateBaseMaterial
    {
        public MITCLinearPlateMaterial() : base()
        {

        }

        public MITCLinearPlateMaterial(MITCLinearPlateMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }
    }
}
