using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TimoshenkoCorotationalFrameMaterial : TimoshenkoFrameBaseMaterial
    {
        public TimoshenkoCorotationalFrameMaterial() : base()
        {

        }

        public TimoshenkoCorotationalFrameMaterial(TimoshenkoCorotationalFrameMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }

    }
}
