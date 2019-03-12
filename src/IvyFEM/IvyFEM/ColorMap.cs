using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    public class ColorMap : IColorMap
    {
        public bool IsFixedMinMax { get; set; } = false;
        public double MinValue { get; set; } = 0;
        public double MaxValue { get; set; } = 0;

        public ColorMap()
        {

        }

        public ColorMap(double minValue, double maxValue)
        {
            IsFixedMinMax = true;
            MinValue = minValue;
            MaxValue = maxValue;
        }

        public ColorMap(ColorMap src)
        {
            IsFixedMinMax = src.IsFixedMinMax;
            MinValue = src.MinValue;
            MaxValue = src.MaxValue;
        }

        public double[] GetColor(double value)
        {
            double[] color = new double[3];

            double r = (value - MinValue) / (MaxValue - MinValue);
            double d = 2.0 * r - 1;
            if (r > 0.75)
            {
                color[0] = 1.0;
                color[1] = 2 - 2 * d;
                color[2] = 0;
            }
            else if (r > 0.50)
            {
                color[0] = -4 * d * d + 4 * d;
                color[1] = 1.0;
                color[2] = 0;
            }
            else if (r > 0.25)
            {
                color[0] = 0.0;
                color[1] = 1.0;
                color[2] = -4 * d * d - 4 * d;
            }
            else
            {
                color[0] = 0.0;
                color[1] = 2 + 2 * d;
                color[2] = 1;
            }
            return color;
        }

    }
}
