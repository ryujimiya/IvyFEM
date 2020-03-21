using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TrussMaterial : Material
    {
        public double Area { get => Values[0]; set => Values[0] = value; }
        public double MassDensity { get => Values[1]; set => Values[1] = value; }
        public double Young { get => Values[2]; set => Values[2] = value; }

        public TrussMaterial()
        {
            int len = 3;
            Values = new double[len];

            Area = 0.0;
            MassDensity = 0.0;
            Young = 0.0;
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
