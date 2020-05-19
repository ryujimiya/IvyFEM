using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class PlateMaterial : Material
    {
        public double Thickness { get => Values[0]; set => Values[0] = value; }
        public double MassDensity { get => Values[1]; set => Values[1] = value; }
        public double Young { get => Values[2]; set => Values[2] = value; }
        public double Poisson { get => Values[3]; set => Values[3] = value; }
        public double ShearCorrectionFactor { get => Values[4]; set => Values[4] = value; }
        public double ShearCoefficient => Young / (2.0 * (1.0 + Poisson));

        public PlateMaterial()
        {
            int len = 5;
            Values = new double[len];

            Thickness = 0.0;
            MassDensity = 0.0;
            Young = 0.0;
            Poisson = 0.0;
            ShearCorrectionFactor = 5.0 / 6.0; // 長方形断面
        }

        public PlateMaterial(PlateMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }
    }
}
