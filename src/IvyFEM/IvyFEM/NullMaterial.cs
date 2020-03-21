using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class NullMaterial : Material
    {
        public NullMaterial()
        {
            int len = 0;
            Values = new double[len];
        }

        public NullMaterial(TrussMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }

    }
}
