using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshBar
    {
        public uint[] V { get; } = new uint[2];
        public uint[] S2 { get; } = new uint[2];
        public uint[] R2 { get; } = new uint[2];
    }

}
