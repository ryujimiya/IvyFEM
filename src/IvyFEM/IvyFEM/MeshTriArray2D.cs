using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshTriArray2D
    {
        public uint Id { get; set; } = 0;
        public uint LCadId { get; set; } = 0;
        public int Layer { get; set; } = 0;
        public IList<MeshTri2D> Tris { get; set; } = new List<MeshTri2D>();
    }

}
