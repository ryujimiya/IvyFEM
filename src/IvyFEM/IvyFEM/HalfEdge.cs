using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class HalfEdge : IIdObject
    {
        public uint Id { get; set; } = 0;
        public uint UVId { get; set; } = 0;
        public uint FHEId { get; set; } = 0; // Previous
        public uint BHEId { get; set; } = 0; // Lator
        public uint OHEId { get; set; } = 0; // Opposite Side
        public uint ULId { get; set; } = 0;
        public uint EId { get; set; } = 0;
        public bool IsSameDir { get; set; } = true;

        public HalfEdge()
        {

        }

        public HalfEdge(uint id, uint uVId, uint fHEId, uint bHEId, uint oHEId, uint uLId)
        {
            Id = id;
            UVId = uVId;
            FHEId = fHEId;
            BHEId = bHEId;
            OHEId = oHEId;
            ULId = uLId;
            EId = 0;
            IsSameDir = true;
        }

        public HalfEdge(HalfEdge src)
        {
            Copy(src);
        }

        public void Copy(IIdObject src)
        {
            HalfEdge srcHE = src as HalfEdge;
            Id = srcHE.Id;
            UVId = srcHE.UVId;
            FHEId = srcHE.FHEId;
            BHEId = srcHE.BHEId;
            OHEId = srcHE.OHEId;
            ULId = srcHE.ULId;
            EId = srcHE.EId;
            IsSameDir = srcHE.IsSameDir;
        }
    }
}
