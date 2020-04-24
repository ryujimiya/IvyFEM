using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class CorotationalFrameMaterial : FrameBaseMaterial
    {
        public CorotationalFrameMaterial() : base()
        {

        }

        public CorotationalFrameMaterial(CorotationalFrameMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }

    }
}
