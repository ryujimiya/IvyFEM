using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MindlinPlateMaterial : PlateBaseMaterial
    {
        public MindlinPlateMaterial() : base()
        {

        }

        public MindlinPlateMaterial(MindlinPlateMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }
    }
}
