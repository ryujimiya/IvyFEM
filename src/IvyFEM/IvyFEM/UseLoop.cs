using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class UseLoop : ICadObject
    {
        public uint Id { get; set; } = 0;
        public uint LId { get; set; } = 0;
        public uint HEId { get; set; } = 0;
        public uint ChildULId { get; set; } = 0;
        public uint ParentULId { get; set; } = 0;

        public UseLoop()
        {

        }

        public UseLoop(uint id, uint hEId, uint childULId, uint parentULId)
        {
            Id = id;
            HEId = hEId;
            ChildULId = childULId;
            ParentULId = parentULId;
        }

        public UseLoop(UseLoop src)
        {
            Copy(src);
        }

        public void Copy(ICadObject src)
        {
            UseLoop srcUL = src as UseLoop;
            Id = srcUL.Id;
            LId = srcUL.LId;
            HEId = srcUL.HEId;
            ChildULId = srcUL.ChildULId;
            ParentULId = srcUL.ParentULId;
        }

        public string Dump()
        {
            string ret = "";
            string CRLF = System.Environment.NewLine;
            ret += "-------------------" + CRLF;
            ret += "UseLoop" + CRLF;
            ret += "  Id = " + Id + CRLF;
            ret += "  LId = " + LId + CRLF;
            ret += "  HEId = " + HEId + CRLF;
            ret += "  ChildULId = " + ChildULId + CRLF;
            ret += "  ParentULId = " + ParentULId + CRLF;
            ret += "-------------------" + CRLF;
            return ret;
        }
    }
}
