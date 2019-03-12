using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public interface IColorMap
    {
        bool IsFixedMinMax { get; set; }
        double MinValue { get; set; }
        double MaxValue { get; set; }

        double[] GetColor(double value);
    }
}
