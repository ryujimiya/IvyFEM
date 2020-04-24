using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TimoshenkoTLFrameMaterial : TimoshenkoFrameBaseMaterial
    {
        public TimoshenkoTLFrameMaterial() : base()
        {

        }

        public TimoshenkoTLFrameMaterial(TimoshenkoTLFrameMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }

    }
}
