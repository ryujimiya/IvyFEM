using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshQuad2D
    {
        public uint[] V { get; } = new uint[4];
        public int[] G2 { get; } = new int[4];
        public uint[] S2 { get; } = new uint[4];
        public uint[] R2 { get; } = new uint[4];
    };
}
