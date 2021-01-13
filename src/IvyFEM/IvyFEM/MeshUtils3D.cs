using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshUtils3D
    {
		/// <summary>
		// 四面体のある面の番号
		// 外からみて半時計周りになるように番号が並べられている．
		/// </summary>
		public static uint[][] ElNoTetFace = new uint[4][] {
			new uint[3] { 1, 2, 3 },
			new uint[3] { 0, 3, 2 },
			new uint[3] { 0, 1, 3 },
			new uint[3] { 0, 2, 1 }
		};
	}
}
