using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class DrawerArray
    {
        public IList<IDrawer> Drawers { get; } = new List<IDrawer>();

        public DrawerArray()
        {

        }

        public void Add(IDrawer drawer)
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

        public void DrawSelection()
        {
            for (uint iDraw = 0; iDraw < Drawers.Count; iDraw++)
            {
                Drawers[(int)iDraw].DrawSelection(iDraw);
            }
        }

        public void AddSelected(int[] selectFlg)
        {
            for (uint iDraw = 0; iDraw < Drawers.Count; iDraw++)
            {
                Drawers[(int)iDraw].AddSelected(selectFlg);
            }
        }

        public void ClearSelected()
        {
            for (uint iDraw = 0; iDraw < Drawers.Count; iDraw++)
            {
                Drawers[(int)iDraw].ClearSelected();
            }
        }

        public BoundingBox3D GetBoundingBox(OpenTK.Matrix3d rot)
        {
            if (Drawers.Count == 0)
            {
                return new BoundingBox3D(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5);
            }
            BoundingBox3D bb = Drawers[0].GetBoundingBox(rot);
            for (int iDraw = 1; iDraw < Drawers.Count; iDraw++)
            {
                bb += Drawers[iDraw].GetBoundingBox(rot);
            }
            return bb;
        }

    }
}
