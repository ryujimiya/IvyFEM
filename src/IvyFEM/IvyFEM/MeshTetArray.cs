using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshTetArray
    {
        public uint Id { get; set; } = 0;
        public uint SCadId { get; set; } = 0;
        public IList<MeshTet> Tets { get; set; } = new List<MeshTet>();

        public MeshTetArray()
        {

        }
    }
}
