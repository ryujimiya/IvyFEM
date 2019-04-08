using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class DiffusionMaterial : Material
    {
        /// <summary>
        /// 密度
        /// </summary>
        public double MassDensity { get => Values[0]; set => Values[0] = value; }
        /// <summary>
        /// 比熱
        /// </summary>
        public double Capacity { get => Values[1]; set => Values[1] = value; }
        /// <summary>
        /// 熱伝導係数
        /// </summary>
        public double DiffusionCoef { get => Values[2]; set => Values[2] = value; }
        /// <summary>
        /// 熱源
        /// </summary>
        public double F { get => Values[3]; set => Values[3] = value; }

        public DiffusionMaterial() : base()
        {
            int len = 4;
            Values = new double[len];
            for (int i = 0; i < len; i++)
            {
                Values[i] = 0.0;
            }
        }
    }
}
