using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshBarArray3D
    {
        public uint Id { get; set; }
        public uint ECadId { get; set; }
        public IList<MeshBar> Bars { get; set; } = new List<MeshBar>();
    }
}
