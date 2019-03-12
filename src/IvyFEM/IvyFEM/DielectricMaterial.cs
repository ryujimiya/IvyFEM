using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class DielectricMaterial : Material
    {
        public double Epxx { get => Values[0]; set => Values[0] = value; }
        public double Epyy { get => Values[1]; set => Values[1] = value; }
        public double Epzz { get => Values[2]; set => Values[2] = value; }
        public double Muxx { get => Values[3]; set => Values[3] = value; }
        public double Muyy { get => Values[4]; set => Values[4] = value; }
        public double Muzz { get => Values[5]; set => Values[5] = value; }

        public DielectricMaterial() : base()
        {
            int len = 6;
            Values = new double[len];
            for (int i = 0; i < len; i++)
            {
                Values[i] = 1.0;
            }
        }
    }
}
