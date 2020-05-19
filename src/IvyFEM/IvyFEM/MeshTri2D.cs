using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshTri2D
    {
        public uint[] V { get; } = new uint[3];
        public int[] G2 { get; } = new int[3];
        public uint[] S2 { get; } = new uint[3];
        public uint[] R2 { get; } = new uint[3];

        public MeshTri2D()
        {

        }

        public MeshTri2D(MeshTri2D src)
        {
            src.V.CopyTo(V, 0);
            src.G2.CopyTo(G2, 0);
            src.S2.CopyTo(S2, 0);
            src.R2.CopyTo(R2, 0);
        }
    }
}
