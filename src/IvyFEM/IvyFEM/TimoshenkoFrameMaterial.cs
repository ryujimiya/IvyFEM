using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TimoshenkoFrameMaterial : TimoshenkoFrameBaseMaterial
    {
        public TimoshenkoFrameMaterial() : base()
        {

        }

        public TimoshenkoFrameMaterial(TimoshenkoFrameMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }

    }
}
