using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TimoshenkoFrameMaterial : Material
    {
        public double Area { get => Values[0]; set => Values[0] = value; }
        public double SecondMomentOfArea { get => Values[1]; set => Values[1] = value; } // Iz
        public double PolarSecondMomentOfArea { get => Values[2]; set => Values[2] = value; } // Ix
        public double MassDensity { get => Values[3]; set => Values[3] = value; }
        public double Young { get => Values[4]; set => Values[4] = value; }
        public double Poisson { get => Values[5]; set => Values[5] = value; }
        public double TimoshenkoShearCoefficient { get => Values[6]; set => Values[6] = value; }
        public double ShearCoefficient => Young / (2.0 * (1.0 + Poisson));

        public TimoshenkoFrameMaterial()
        {
            int len = 7;
            Values = new double[len];

            Area = 0.0;
            SecondMomentOfArea = 0.0;
            PolarSecondMomentOfArea = 0.0;
            MassDensity = 0.0;
            Young = 0.0;
            Poisson = 0.0;
            TimoshenkoShearCoefficient = 5.0 / 6.0; // 長方形断面
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
