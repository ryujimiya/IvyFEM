using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class BeamMaterial : FrameBaseMaterial
    {
        public BeamMaterial() : base()
        {

        }

        public BeamMaterial(BeamMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }

    }
}
