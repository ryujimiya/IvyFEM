using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshBarArray2D
    {
        public uint Id { get; set; } = 0;
        public uint ECadId { get; set; } = 0;
        public uint[] SEId { get; } = new uint[2];
        public uint[] LRId { get; } = new uint[2];
        public int Layer { get; set; } = 0;
        public IList<MeshBar> Bars { get; set; } = new List<MeshBar>();

        public MeshBarArray2D()
        {

        }
    }
}
