using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshTriArray3D
    {
        public uint Id { get; set; } = 0;
        public uint LCadId { get; set; } = 0;
        public IList<MeshTri3D> Tris { get; set; } = new List<MeshTri3D>();

        public MeshTriArray3D()
        {

        }
    }
}
