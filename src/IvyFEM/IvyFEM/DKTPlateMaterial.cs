using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class DKTPlateMaterial : PlateBaseMaterial
    {
        public DKTPlateMaterial() : base()
        {

        }

        public DKTPlateMaterial(DKTPlateMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }
    }
}
