using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class HelmholtzMaterial : Material
    {
        /// <summary>
        /// 音速
        /// </summary>
        public double Velocity { get => Values[0]; set => Values[0] = value; }
        /// <summary>
        /// 励振源 F
        /// </summary>
        public System.Numerics.Complex F
        {
            get => new System.Numerics.Complex(Values[1], Values[2]);
            set { Values[1] = value.Real; Values[2] = value.Imaginary; }
        }

        public HelmholtzMaterial() : base()
        {
            int len = 3;
            Values = new double[len];
            for (int i = 0; i < len; i++)
            {
                Values[i] = 0.0;
            }
        }
    }
}
