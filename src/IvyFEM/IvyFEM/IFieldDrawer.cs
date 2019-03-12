using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public interface IFieldDrawer : IDrawer 
    {
        void Update(FEWorld world);
    }
}
