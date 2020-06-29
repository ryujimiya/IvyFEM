using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class EddyCurrentMaterial : Material
    {
        /// <summary>
        /// 比透磁率
        /// </summary>
        public double Mu { get => Values[0]; set => Values[0] = value; }
        /// <summary>
        /// 導電率
        /// </summary>
        public double Sigma { get => Values[1]; set => Values[1] = value; }
        /// <summary>
        /// 印加電圧の勾配
        /// </summary>
        public double GradPhi { get => Values[2]; set => Values[2] = value; }

        public EddyCurrentMaterial() : base()
        {
            int len = 3;
            Values = new double[len];
            for (int i = 0; i < len; i++)
            {
                Values[i] = 0.0;
            }
        }

        public EddyCurrentMaterial(EddyCurrentMaterial src) :  base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }
    }
}
