using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class PoissonMaterial : Material
    {
        /// <summary>
        /// 係数 α
        ///   α(∇^2(φ) = -F
        /// </summary>
        public double Alpha { get => Values[0]; set => Values[0] = value; }
        /// <summary>
        /// 励振源 F
        /// </summary>
        public double F { get => Values[1]; set => Values[1] = value; }

        public PoissonMaterial() : base()
        {
            int len = 2;
            Values = new double[len];
            for (int i = 0; i < len; i++)
            {
                Values[i] = 0.0;
            }
        }
    }
}
