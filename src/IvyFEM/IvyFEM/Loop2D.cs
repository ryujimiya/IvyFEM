using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Loop2D : IObject
    {
        public uint Layer { get; set; } = 0;
        public double[] Color { get; set; } = new double[3];

        public Loop2D()
        {
            Layer = 0;
            //Color = new double[3] { 0.8, 0.8, 0.8 };
            Color = new double[3] { 0.2, 0.2, 0.2 };
        }

        public Loop2D(Loop2D src)
        {
            Copy(src);
        }

        public void Copy(IObject src)
        {
            Loop2D srcLoop = src as Loop2D;

            Layer = srcLoop.Layer;
            srcLoop.Color.CopyTo(Color, 0);
        }
    }
}
