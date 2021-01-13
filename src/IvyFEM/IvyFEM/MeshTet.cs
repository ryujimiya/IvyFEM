using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshTet
    {
		// 頂点のIndex
		public uint[] V { get; set; } = new uint[4];
		// 隣接要素のIndex
		public uint[] S { get; set; } = new uint[4];
		
		public int[] G { get; set; } = new int[4];
		public uint[] F { get; set; } = new uint[4];
		public int OldT { get; set; } = 0;

		public MeshTet()
        {

        }

		public MeshTet(MeshTet src)
        {
			src.V.CopyTo(V, 0);
			src.S.CopyTo(S, 0);
			src.G.CopyTo(G, 0);
			src.F.CopyTo(F, 0);
			OldT = src.OldT;
        }
	}
}
