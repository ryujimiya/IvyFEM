using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class UseVertex : IIdObject
    {
        public uint Id { get; set; } = 0;
        public uint HEId { get; set; } = 0;
        public uint VId { get; set; } = 0;

        public UseVertex()
        {

        }

        public UseVertex(uint id, uint hEId)
        {
            Id = id;
            HEId = hEId;
            VId = 0;
        }

        public UseVertex(UseVertex src)
        {
            Copy(src);
        }

        public void Copy(IIdObject src)
        {
            UseVertex srcUV = src as UseVertex;
            Id = srcUV.Id;
            HEId = srcUV.HEId;
            VId = srcUV.VId;
        }
    }
}
