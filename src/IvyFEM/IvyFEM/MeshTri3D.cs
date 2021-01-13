using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshTri3D
    {
        // 頂点Index
        public uint[] V { get; } = new uint[3];
        // 隣接する要素配列ID(-1:隣接要素なし、-2:自分の要素配列に隣接)
        public int[] G2 { get; } = new int[3];
        // 隣接要素Index
        public int[] S2 { get; } = new int[3];
        // 隣接関係
        public int[] R2 { get; } = new int[3];
        ////////////////////////////////////////////////
        public int[] Sf { get; } = new int[2];
        public int[] Gf { get; } = new int[2];
        public int[] Df { get; } = new int[2];

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
