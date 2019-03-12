using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class FieldDrawerArray : DrawerArray
    {
        public FieldDrawerArray()
        {

        }

        public void Update(FEWorld world)
        {
            for (int iDraw = 0; iDraw < Drawers.Count; iDraw++)
            {
                IFieldDrawer drawer = Drawers[iDraw] as IFieldDrawer;
                drawer.Update(world);
            }
        }
    }
}
