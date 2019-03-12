using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class HalfEdge : ICadObject
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

        public void Copy(ICadObject src)
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

        public string Dump()
        {
            string ret = "";
            string CRLF = System.Environment.NewLine;
            ret += "-------------------" + CRLF;
            ret += "HalfEdge" + CRLF;
            ret += "  Id = " + Id + CRLF;
            ret += "  UVId = " + UVId + CRLF;
            ret += "  FHEId = " + FHEId + CRLF;
            ret += "  BHEId = " + BHEId + CRLF;
            ret += "  OHEId = " + OHEId + CRLF;
            ret += "  ULId = " + ULId + CRLF;
            ret += "  EId = " + EId + CRLF;
            ret += "-------------------" + CRLF;
            return ret;
        }

    }

}
