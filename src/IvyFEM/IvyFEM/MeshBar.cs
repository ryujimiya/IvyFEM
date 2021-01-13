using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshBar
    {
        // 頂点Index
        public uint[] V { get; } = new uint[2];
        // 隣接要素Index
        public uint[] S2 { get; } = new uint[2];
        // 隣接関係
        public uint[] R2 { get; } = new uint[2];
    }

}
