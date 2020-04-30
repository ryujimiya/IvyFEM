using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class FrameBaseMaterial : Material
    {
        public double Area { get => Values[0]; set => Values[0] = value; }
        public double SecondMomentOfArea { get => Values[1]; set => Values[1] = value; } // Iz
        public double MassDensity { get => Values[2]; set => Values[2] = value; }
        public double Young { get => Values[3]; set => Values[3] = value; }
        public double Poisson { get => Values[4]; set => Values[4] = value; }
        public double ShearCoefficient => Young / (2.0 * (1.0 + Poisson));

        public FrameBaseMaterial()
        {
            int len = 5;
            Values = new double[len];

            Area = 0.0;
            SecondMomentOfArea = 0.0;
            MassDensity = 0.0;
            Young = 0.0;
            Poisson = 0.0;
        }

        public FrameBaseMaterial(FrameBaseMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }

    }
}
