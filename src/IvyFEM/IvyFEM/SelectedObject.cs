using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class SelectedObject
    {
        public uint NameDepth { get; set; } = 0;
        public int[] Name { get; } = new int[4];
        public OpenTK.Vector3d PickedPos { get; set; } = new OpenTK.Vector3d();
    }
}
