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
        public int QuantityIdBaseOffset { get; set; } = 0;
        public ElementType Type { get; protected set; } = ElementType.NotSet;
        public FiniteElementType FEType { get; set; } = FiniteElementType.ScalarLagrange;
        public int Order { get; set; } = -1;
        public uint VertexCount { get; protected set; } = 0;
        public uint NodeCount { get; protected set; } = 0;
        public uint EdgeCount { get; protected set; } = 0;
        public int[] VertexCoordIds { get; protected set; } = null;
        public int[] NodeCoordIds { get; protected set; } = null;
        public int[][] EdgeCoordIdss { get; protected set; } = null;
        public double[][] Displacements { get; protected set; } = null;
        public uint MaterialId { get; set; } = 0;
        public uint MeshId { get; set; } = 0;
        public int MeshElemId { get; set; } = -1;
        protected IInterpolate Interpolate = null;

        public FE()
        {

        }

        public FE(FE src)
        {
            Copy(src);
            CreateInterpolate();
        }

        public virtual void Copy(IObject src)
        {
            FE srcFE = src as FE;

            World = srcFE.World; // shallow copy
            QuantityId = srcFE.QuantityId;
            Type = srcFE.Type;
            FEType = srcFE.FEType;
            Order = srcFE.Order;
            VertexCount = srcFE.VertexCount;
            NodeCount = srcFE.NodeCount;
            EdgeCount = srcFE.EdgeCount;
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
            EdgeCoordIdss = null;
            if (srcFE.EdgeCoordIdss != null)
            {
                int edgeCnt = srcFE.EdgeCoordIdss.Length;
                System.Diagnostics.Debug.Assert(EdgeCount == edgeCnt);
                EdgeCoordIdss = new int[edgeCnt][];
                for (int i = 0; i < edgeCnt; i++)
                {
                    EdgeCoordIdss[i] = new int[srcFE.EdgeCoordIdss[i].Length];
                    srcFE.EdgeCoordIdss[i].CopyTo(EdgeCoordIdss[i], 0);
                }
                srcFE.EdgeCoordIdss.CopyTo(EdgeCoordIdss, 0);
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

        protected void CreateInterpolate()
        {
            if (FEType == FiniteElementType.ScalarLagrange)
            {
                CreateLagrangeInterpolate();
            }
            else if (FEType == FiniteElementType.ScalarBell)
            {
                CreateBellInterpolate();
            }
            else if (FEType == FiniteElementType.Edge)
            {
                CreateEdgeInterpolate();
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }

        }

        protected void CreateLagrangeInterpolate()
        {
            System.Diagnostics.Debug.Assert(FEType == FiniteElementType.ScalarLagrange);
            if (this is LineFE)
            {
                LineFE thisLineFE = this as LineFE;
                if (Order == 1)
                {
                    Interpolate = new LineFE1stInterpolate(thisLineFE);
                }
                else if (Order == 2)
                {
                    Interpolate = new LineFE2ndInterpolate(thisLineFE);
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
            }
            else if (this is TriangleFE)
            {
                TriangleFE thisTriFE = this as TriangleFE;
                if (Order == 1)
                {
                    Interpolate = new TriangleFE1stInterpolate(thisTriFE);
                }
                else if (Order == 2)
                {
                    Interpolate = new TriangleFE2ndInterpolate(thisTriFE);
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
        }

        protected void CreateBellInterpolate()
        {
            System.Diagnostics.Debug.Assert(FEType == FiniteElementType.ScalarBell);
            if (this is LineFE)
            {
                LineFE thisLineFE = this as LineFE;
                if (Order == 1)
                {
                    // for Bell element 5次
                    // 暫定：Lagrange線要素で代用.当然ながら形状関数は正しくない
                    Interpolate = new LineFE1stInterpolate(thisLineFE);
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
            }
            else if (this is TriangleFE)
            {
                TriangleFE thisTriFE = this as TriangleFE;
                if (Order == 5)
                {
                    // for Bell element 5次
                    Interpolate = new TriangleFEBellInterpolate(thisTriFE);
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
        }

        protected void CreateEdgeInterpolate()
        {
            System.Diagnostics.Debug.Assert(FEType == FiniteElementType.Edge);
            if (this is LineFE)
            {
                // 暫定: Lagrange線要素で代用
                LineFE thisLineFE = this as LineFE;
                if (Order == 1)
                {
                    Interpolate = new LineFE1stInterpolate(thisLineFE);
                }
                else if (Order == 2)
                {
                    Interpolate = new LineFE2ndInterpolate(thisLineFE);
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
            }
            else if (this is TriangleFE)
            {
                TriangleFE thisTriFE = this as TriangleFE;
                if (Order == 1)
                {
                    Interpolate = new TriangleFEEdge1stInterpolate(thisTriFE);
                }
                else if (Order == 2)
                {
                    Interpolate = new TriangleFEEdge2ndInterpolate(thisTriFE);
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
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

        public void SetEdgeCoordIdsFromNodeCoordIds()
        {
            if (!(Interpolate is IEdgeInterpolate))
            {
                return;
            }
            IEdgeInterpolate edgeInterpolate = Interpolate as IEdgeInterpolate;
            int edgeCnt = (int)edgeInterpolate.GetEdgeCount();
            System.Diagnostics.Debug.Assert(edgeCnt == EdgeCount);
            int[][] edgePointId = edgeInterpolate.GetEdgePointIdss();
            EdgeCoordIdss = new int[edgeCnt][];
            for (int i = 0; i < edgeCnt; i++)
            {
                int[] edgeCoIds = new int[2];
                int pointId1 = edgePointId[i][0];
                int pointId2 = edgePointId[i][1];

                edgeCoIds[0] = NodeCoordIds[pointId1];
                edgeCoIds[1] = NodeCoordIds[pointId2];
                EdgeCoordIdss[i] = edgeCoIds;
            }
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
