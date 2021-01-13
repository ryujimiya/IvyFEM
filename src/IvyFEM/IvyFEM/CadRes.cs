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

    public class ConnectVertexRes
    {
        public uint VId1 { get; set; } = 0;
        public uint VId2 { get; set; } = 0;
        public uint LId { get; set; } = 0;
        public uint AddEId { get; set; } = 0;
        public uint AddLId { get; set; } = 0;
        public bool IsLeftAddL { get; set; } = true;
    }

    public class AddPolygonRes
    {
        public uint AddLId { get; set; } = 0;
        public IList<uint> VIds { get; } = new List<uint>();
        public IList<uint> EIds { get; } = new List<uint>();
    }

    public class AddSurfaceRes
    {
        public IList<uint> AddVIds { get; set; } = new List<uint>();
        public IList<uint> AddEIds { get; set; } = new List<uint>();
        public IList<uint> AddLIds { get; set; } = new List<uint>();
    }
}
