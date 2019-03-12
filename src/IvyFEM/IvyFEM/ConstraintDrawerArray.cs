using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class ConstraintDrawerArray
    {
        public IList<ConstraintDrawer> Drawers { get; } = new List<ConstraintDrawer>();

        public ConstraintDrawerArray()
        {

        }

        public void Add(ConstraintDrawer drawer)
        {
            System.Diagnostics.Debug.Assert(drawer != null);
            Drawers.Add(drawer);
        }

        public void Clear()
        {
            Drawers.Clear();
        }

        public void Draw()
        {
            for (int iDraw = 0; iDraw < Drawers.Count; iDraw++)
            {
                Drawers[iDraw].Draw();
            }
        }
    }
}
