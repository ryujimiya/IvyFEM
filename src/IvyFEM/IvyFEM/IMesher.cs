using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public interface IMesher
    {
        uint Dimension { get; }
        void GetCoords(out IList<double> coord);
        IList<uint> GetIds();
        bool IsId(uint id);
        IList<uint> GetIncludeElementIds(uint id);
        void GetInfo(uint id, out uint cadId, out int layer);
        void GetConnectivity(uint id, out MeshType meshType, out int[] vertexs);
        bool GetMeshInfo(uint id,
            out uint elemCount, out MeshType meshType, out int loc, out uint cadId);
        uint GetIdFromCadId(uint cadId, CadElementType cadType);
    }
}
