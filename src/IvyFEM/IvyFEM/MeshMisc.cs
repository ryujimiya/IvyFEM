using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshTypeLoc
    {
        public int Type { get; set; } = 0;
        public int Loc { get; set; } = -1;
        public IList<uint> IncludeRelations { get; set; } = new List<uint>();

        public MeshTypeLoc()
        {

        }

        public MeshTypeLoc(MeshTypeLoc src)
        {
            Type = src.Type;
            Loc = src.Loc;
            IncludeRelations = new List<uint>(src.IncludeRelations);
        }
    }

    public class CadIdELen
    {
        public uint CadId { get; set; } = 0;
        public double ELen { get; set; } = 0.0;

        public CadIdELen()
        {

        }

        public CadIdELen(CadIdELen src)
        {
            CadId = src.CadId;
            ELen = src.ELen;
        }
    }
}
