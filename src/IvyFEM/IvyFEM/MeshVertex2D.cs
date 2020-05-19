using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshVertex2D
    {
        public uint Id { get; set; } = 0;
        public uint VCadId { get; set; } = 0;
        public int Layer { get; set; } = 0;
        public uint V { get; set; } = 0;
    }
}
