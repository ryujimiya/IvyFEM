using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MITCStVenantPlateMaterial : PlateBaseMaterial
    {
        public MITCStVenantPlateMaterial() : base()
        {

        }

        public MITCStVenantPlateMaterial(MITCStVenantPlateMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }
    }
}
