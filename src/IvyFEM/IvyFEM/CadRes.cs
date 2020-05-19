using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class AddVertexRes
    {
        public uint AddVId { get; set; } = 0;
        public uint AddEId { get; set; } = 0;
    }

    public class AddPolygonRes
    {
        public uint AddLId { get; set; } = 0;
        public IList<uint> VIds { get; } = new List<uint>();
        public IList<uint> EIds { get; } = new List<uint>();

        public AddPolygonRes()
        {

        }

        public AddPolygonRes(AddPolygonRes src)
        {
            this.AddLId = src.AddLId;
            this.VIds.Clear();
            foreach (var id in src.VIds)
            {
                this.VIds.Add(id);
            }
            this.EIds.Clear();
            foreach (var id in src.EIds)
            {
                this.EIds.Add(id);
            }
        }
    }
}
