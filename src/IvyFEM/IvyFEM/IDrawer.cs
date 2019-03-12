using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public interface IDrawer
    {
        RotMode SutableRotMode { get; }
        bool IsAntiAliasing { get; set; }

        BoundingBox3D GetBoundingBox(OpenTK.Matrix3d rot);

        void DrawSelection(uint idraw);
        void Draw();
        void AddSelected(int[] selectFlag);
        void ClearSelected();
    }
}
