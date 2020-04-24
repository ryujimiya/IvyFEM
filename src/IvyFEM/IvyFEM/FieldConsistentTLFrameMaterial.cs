using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class FieldConsistentTLFrameMaterial : FrameBaseMaterial
    {
        public FieldConsistentTLFrameMaterial() : base()
        {

        }

        public FieldConsistentTLFrameMaterial(FieldConsistentTLFrameMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }

    }
}
