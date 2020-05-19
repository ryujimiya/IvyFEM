using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshTri3D
    {
        public uint[] V { get; } = new uint[3];
        public int[] G2 { get; } = new int[3];
        public uint[] S2 { get; } = new uint[3];
        public uint[] R2 { get; } = new uint[3];
        public int[] Sf { get; } = new int[2];
        public int[] Gf { get; } = new int[2];
        public uint[] Df { get; } = new uint[2];

        public MeshTri3D()
        {

        }

        public MeshTri3D(MeshTri3D src)
        {
            src.V.CopyTo(V, 0);
            src.G2.CopyTo(G2, 0);
            src.S2.CopyTo(S2, 0);
            src.R2.CopyTo(R2, 0);
            src.Sf.CopyTo(Sf, 0);
            src.Gf.CopyTo(Gf, 0);
            src.Df.CopyTo(Df, 0);
        }
    }
}
