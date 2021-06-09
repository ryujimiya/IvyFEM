using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class LoopEdgeItr
    {
        private bool IsValid = false;
        private bool IsInitial = false;
        public bool IsChildEnd { get; private set; } = false;
        private uint HEId = 0;
        private uint ULId = 0;
        private BRep2D BRep2D = new BRep2D();

        public LoopEdgeItr(BRep2D bRep2D, uint lId)
        {
            BRep2D = bRep2D;
            IsValid = false;
            ULId = 0;
            {
                System.Diagnostics.Debug.Assert(BRep2D.Loop2UseLoop.ContainsKey(lId));
                ULId = BRep2D.Loop2UseLoop[lId];
            }
            //System.Diagnostics.Debug.Assert(BRep2D.BRep.IsUseLoopId(ULId));
            if (BRep2D.BRep.IsUseLoopId(ULId))
            {
                UseLoop uL = BRep2D.BRep.GetUseLoop(ULId);
                System.Diagnostics.Debug.Assert(uL.Id == ULId);
                HEId = uL.HEId;
                IsInitial = true;
                IsValid = true;
                IsChildEnd = false;
            }
            else
            {
                // 不正な情報
                System.Diagnostics.Debug.Assert(false);
                IsInitial = true;
                IsValid = false;
                IsChildEnd = true;
            }
        }

        public LoopEdgeItr(BRep2D bRep2D, uint hEId, uint uLId)
        {
            BRep2D = bRep2D;
            IsValid = false;
            ULId = uLId;
            HEId = hEId;
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsUseLoopId(ULId));
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(HEId));
            HalfEdge hE = BRep2D.BRep.GetHalfEdge(HEId);
            System.Diagnostics.Debug.Assert(hE.ULId == ULId);
            IsInitial = false;
            IsValid = true;
            IsChildEnd = false;
        }

        public LoopEdgeItr(LoopEdgeItr src)
        {
            IsValid = src.IsValid;
            IsInitial = src.IsInitial;
            IsChildEnd = src.IsChildEnd;
            HEId = src.HEId;
            ULId = src.ULId;
            BRep2D.Copy(src.BRep2D);
        }

        public void Begin()
        {
            if (!IsValid)
            {
                return;
            }
            IsInitial = true;
            UseLoop uL = BRep2D.BRep.GetUseLoop(ULId);
            System.Diagnostics.Debug.Assert(uL.Id == ULId);
            HEId = uL.HEId;
        }

        // ++のオブジェクト生成をしないバージョン
        public void Next()
        {
            if (!IsValid)
            {
                return;
            }
            IsInitial = false;
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(HEId));
            HalfEdge hE = BRep2D.BRep.GetHalfEdge(HEId);
            HEId = hE.FHEId;
        }

        public bool IsEnd()
        {
            if (!IsValid)
            {
                return false;
            }
            if (IsInitial)
            {
                return false;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsUseLoopId(ULId));
            UseLoop uL = BRep2D.BRep.GetUseLoop(ULId);
            System.Diagnostics.Debug.Assert(uL.Id == ULId);
            if (HEId == uL.HEId)
            {
                return true;
            }
            return false;
        }

        public bool GetEdgeId(out uint eId, out bool isSameDir)
        {
            eId = 0;
            isSameDir = true;

            if (!IsValid)
            {
                return false;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(HEId));
            HalfEdge hE = BRep2D.BRep.GetHalfEdge(HEId);
            eId = hE.EId;
            isSameDir = hE.IsSameDir;
            if (eId == 0)
            {
                return false;
            }
            return true;
        }

        public uint GetVertexId()
        {
            if (!IsValid)
            {
                return 0;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(HEId));
            HalfEdge hE = BRep2D.BRep.GetHalfEdge(HEId);
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsUseVertexId(hE.UVId));
            UseVertex uV = BRep2D.BRep.GetUseVertex(hE.UVId);
            uint vId = uV.VId;
            System.Diagnostics.Debug.Assert(BRep2D.IsElementId(CadElementType.Vertex, vId));
            return vId;
        }

        public uint GetAheadVertexId()
        {
            if (!IsValid)
            {
                return 0;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(HEId));
            HalfEdge hE = BRep2D.BRep.GetHalfEdge(HEId);
            uint fHEId = hE.FHEId;
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(fHEId));
            HalfEdge fHE = BRep2D.BRep.GetHalfEdge(fHEId);
            uint fUVId = fHE.UVId;
            UseVertex fUV = BRep2D.BRep.GetUseVertex(fUVId);
            uint vId = fUV.VId;
            System.Diagnostics.Debug.Assert(BRep2D.IsElementId(CadElementType.Vertex, vId));
            return vId;
        }

        public uint GetBehindVertexId()
        {
            if (!IsValid)
            {
                return 0;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(HEId));
            HalfEdge hE = BRep2D.BRep.GetHalfEdge(HEId);
            uint bHEId = hE.BHEId;
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(bHEId));
            HalfEdge bHE = BRep2D.BRep.GetHalfEdge(bHEId);
            uint bUVId = bHE.UVId;
            UseVertex bUV = BRep2D.BRep.GetUseVertex(bUVId);
            uint vId = bUV.VId;
            System.Diagnostics.Debug.Assert(BRep2D.IsElementId(CadElementType.Vertex, vId));
            return vId;
        }

        public bool ShiftChildLoop()
        {
            if (!IsValid)
            {
                return false;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsUseLoopId(ULId));
            UseLoop parentUL = BRep2D.BRep.GetUseLoop(ULId);
            System.Diagnostics.Debug.Assert(parentUL.Id == ULId);
            if (parentUL.ChildULId == 0)
            {
                IsChildEnd = true;
                return false;
            }
            ULId = parentUL.ChildULId;
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsUseLoopId(ULId));
            UseLoop childUL = BRep2D.BRep.GetUseLoop(ULId);
            System.Diagnostics.Debug.Assert(childUL.Id == ULId);
            HEId = childUL.HEId;
            IsInitial = true;
            IsValid = true;
            return true;
        }

        public uint GetUseLoopType()
        {
            return BRep2D.GetUseLoopType(ULId);
        }

		public uint GetHalfEdgeId()
        {
            return HEId;
        }

        public uint GetUseLoopId()
        {
            return ULId;
        }

        public uint GetLoopId()
        {
            if (!IsValid)
            {
                return 0;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsUseLoopId(ULId));
            UseLoop uL = BRep2D.BRep.GetUseLoop(ULId);
            System.Diagnostics.Debug.Assert(uL.Id == ULId);
            return uL.LId;
        }

        public bool IsParent()
        {
            if (!IsValid)
            {
                return false;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsUseLoopId(ULId));
            UseLoop uL = BRep2D.BRep.GetUseLoop(ULId);
            System.Diagnostics.Debug.Assert(uL.Id == ULId);
            return (uL.ParentULId == ULId);
        }

        public bool IsSameUseLoop(LoopEdgeItr lItr)
        {
            return (lItr.ULId == ULId);
        }

        public bool IsEdgeBothSideSameLoop()
        {
            if (!IsValid)
            {
                return true;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(HEId));
            HalfEdge hE = BRep2D.BRep.GetHalfEdge(HEId);
            uint oHEId = hE.OHEId;
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(oHEId));
            HalfEdge oHE = BRep2D.BRep.GetHalfEdge(oHEId);
            return hE.ULId == oHE.ULId;
        }

        public uint GetUseLoopVertexCount()
        {
            if (!IsValid)
            {
                return 0;
            }
            uint hEId0;
            {
                System.Diagnostics.Debug.Assert(BRep2D.BRep.IsUseLoopId(ULId));
                UseLoop uL = BRep2D.BRep.GetUseLoop(ULId);
                hEId0 = uL.HEId;
            }
            uint hEId = hEId0;
            uint iCnt = 0;
            while (true)
            {
                iCnt++;
                HalfEdge hE = BRep2D.BRep.GetHalfEdge(hEId);
                uint fHEId = hE.FHEId;
                if (fHEId == hEId0)
                {
                    break;
                }
                {
                    System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(fHEId));
                    HalfEdge fHE = BRep2D.BRep.GetHalfEdge(fHEId);
                    System.Diagnostics.Debug.Assert(fHE.ULId == ULId);
                }
                hEId = fHEId;
            }
            return iCnt;
        }


    }
}
