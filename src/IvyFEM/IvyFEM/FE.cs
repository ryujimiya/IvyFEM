using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class FE : IObject
    {
        internal FEWorld World { get; set; }
        public int QuantityId { get; set; } = -1;
        public ElementType Type { get; protected set; } = ElementType.NotSet;
        public int Order { get; set; } = -1;
        public uint VertexCount { get; protected set; } = 0;
        public uint NodeCount { get; protected set; } = 0;
        public int[] VertexCoordIds { get; protected set; } = null;
        public int[] NodeCoordIds { get; protected set; } = null;
        public double[][] Displacements { get; protected set; } = null;
        public uint MaterialId { get; set; } = 0;
        public uint MeshId { get; set; } = 0;
        public int MeshElemId { get; set; } = -1;

        public FE()
        {

        }

        public FE(FE src)
        {
            Copy(src);
        }

        public virtual void Copy(IObject src)
        {
            FE srcFE = src as FE;

            World = srcFE.World; // shallow copy
            QuantityId = srcFE.QuantityId;
            Type = srcFE.Type;
            Order = srcFE.Order;
            NodeCount = srcFE.NodeCount;
            VertexCoordIds = null;
            if (srcFE.VertexCoordIds != null)
            {
                VertexCoordIds = new int[srcFE.VertexCoordIds.Length];
                srcFE.VertexCoordIds.CopyTo(VertexCoordIds, 0);
            }
            NodeCoordIds = null;
            if (srcFE.NodeCoordIds != null)
            {
                NodeCoordIds = new int[srcFE.NodeCoordIds.Length];
                srcFE.NodeCoordIds.CopyTo(NodeCoordIds, 0);
            }
            Displacements = null;
            if (srcFE.Displacements != null)
            {
                Displacements = new double[srcFE.Displacements.Length][];
                for (int i = 0; i < srcFE.Displacements.Length; i++)
                {
                    double[] srcU = srcFE.Displacements[i];
                    double[] u = new double[srcU.Length];
                    srcU.CopyTo(u, 0);
                    Displacements[i] = u;
                }
            }
            MaterialId = srcFE.MaterialId;
            MeshId = srcFE.MeshId;
            MeshElemId = srcFE.MeshElemId;
        }

        public void SetVertexCoordIds(int[] vertexCoordIds)
        {
            System.Diagnostics.Debug.Assert(VertexCount == vertexCoordIds.Length);
            VertexCoordIds = new int[VertexCount];
            vertexCoordIds.CopyTo(VertexCoordIds, 0);
        }

        public void SetNodeCoordIds(int[] nodeCoordIds)
        {
            System.Diagnostics.Debug.Assert(NodeCount == nodeCoordIds.Length);
            NodeCoordIds = new int[NodeCount];
            nodeCoordIds.CopyTo(NodeCoordIds, 0);
        }

        public void SetDisplacements(double[][] displacements)
        {
            System.Diagnostics.Debug.Assert(NodeCount == displacements.Length);
            Displacements = new double[NodeCount][];
            for (int i = 0; i < NodeCount; i++)
            {
                double[] srcU = displacements[i];
                double[] u = new double[srcU.Length];
                srcU.CopyTo(u, 0);
                Displacements[i] = u;
            }
        }
    }
}
