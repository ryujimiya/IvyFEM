using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MooneyRivlinHyperelasticMaterial : Material
    {
        public bool IsCompressible { get => IntValues[0] == 1; set => IntValues[0] = (value ? 1 : 0); }
        public double MassDensity { get => Values[0]; set => Values[0] = value; }
        public double GravityX { get => Values[1]; set => Values[1] = value; }
        public double GravityY { get => Values[2]; set => Values[2] = value; }
        public double D1 { get => Values[3]; set => Values[3] = value; }
        public double C1 { get => Values[4]; set => Values[4] = value; }
        public double C2 { get => Values[5]; set => Values[5] = value; }

        public MooneyRivlinHyperelasticMaterial()
        {
            int len = 6;
            Values = new double[len];
            int intLen = 1;
            IntValues = new int[intLen];

            MassDensity = 1.0;
            GravityX = 0.0;
            GravityY = 0.0;
            IsCompressible = false;
            D1 = 1.0;
            C1 = 0.0;
            C2 = 0.0;
        }

        public MooneyRivlinHyperelasticMaterial(MooneyRivlinHyperelasticMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }
    }
}
