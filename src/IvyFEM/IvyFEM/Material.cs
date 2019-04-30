using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Material : IObject
    {
        public double[] Values { get; protected set; } = null;
        public int[] IntValues { get; protected set; } = null;

        public Material()
        {

        }

        public Material(Material src)
        {
            Copy(src);
        }

        public virtual void Copy(IObject src)
        {
            Material srcMa = src as Material;
            Values = null;
            if (srcMa.Values != null)
            {
                Values = new double[srcMa.Values.Length];
                srcMa.Values.CopyTo(Values, 0);
            }
            IntValues = null;
            if (srcMa.IntValues != null)
            {
                IntValues = new int[srcMa.IntValues.Length];
                srcMa.IntValues.CopyTo(IntValues, 0);
            }
        }

    }
}
