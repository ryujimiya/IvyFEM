using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MITCMooneyRivlinPlateMaterial : Material
    {
        public double Thickness { get => Values[0]; set => Values[0] = value; }
        public double MassDensity { get => Values[1]; set => Values[1] = value; }
        public double C1 { get => Values[2]; set => Values[2] = value; }
        public double C2 { get => Values[3]; set => Values[3] = value; }

        public MITCMooneyRivlinPlateMaterial() : base()
        {
            int len = 4;
            Values = new double[len];

            Thickness = 0.0;
            MassDensity = 1.0;
            C1 = 0.0;
            C2 = 0.0;
        }

        public MITCMooneyRivlinPlateMaterial(MITCMooneyRivlinPlateMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }
    }
}
