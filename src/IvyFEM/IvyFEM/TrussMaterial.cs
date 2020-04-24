using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TrussMaterial : FrameBaseMaterial
    {
        public TrussMaterial() : base()
        {

        }

        public TrussMaterial(TrussMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }

    }
}
