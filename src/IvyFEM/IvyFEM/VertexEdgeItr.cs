using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class VertexEdgeItr
    {
        private bool IsValid = false;
        private bool IsInitial = false;
        private uint UVId = 0;
        private uint HEId = 0;

        //---------------------------
        private int HEPos = 0;
        private IList<uint> HEIds = null;
        //---------------------------

        private BRep2D BRep2D = new BRep2D();

        public VertexEdgeItr(BRep2D bRep2D, uint vId)
        {
            BRep2D = bRep2D;
            IsValid = false;
            UVId = vId;
            /*
            //------------------------------------------------------
            // original
            UseVertex uV = BRep2D.BRep.GetUseVertex(UVId);
            System.Diagnostics.Debug.Assert(uV.Id == UVId);
            HEId = uV.HEId;
            //------------------------------------------------------
            */
            //------------------------------------------------------
            HEIds = BRep2D.BRep.FindHalfEdgeByUseVertex(UVId);
            HEPos = 0;
            HEId = HEIds[HEPos];
            //------------------------------------------------------

            IsInitial = true;
            IsValid = true;
        }

        public VertexEdgeItr(VertexEdgeItr src)
        {
            IsValid = src.IsValid;
            IsInitial = src.IsInitial;
            UVId = src.UVId;
            HEId = src.HEId;
            BRep2D.Copy(src.BRep2D);
        }

        public void Begin()
        {
            if (!IsValid)
            {
                return;
            }
            IsInitial = true;
            /*
            //------------------------------------------------------
            // original
            UseVertex uV = BRep2D.BRep.GetUseVertex(UVId);
            System.Diagnostics.Debug.Assert(uV.Id == UVId);
            HEId = uV.HEId;
            //------------------------------------------------------
            */
            //------------------------------------------------------
            HEPos = 0;
            HEId = HEIds[HEPos];
            //------------------------------------------------------
        }

        // ++のオブジェクトを生成しないバージョン
        public void Next()
        {
            if (!IsValid)
            {
                return;
            }
            IsInitial = false;
            /*
            //------------------------------------------------------
            // original
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(HEId));
            HalfEdge hE = BRep2D.BRep.GetHalfEdge(HEId);
            uint bHEId = hE.BHEId;
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(bHEId));
            HalfEdge bHE = BRep2D.BRep.GetHalfEdge(bHEId);
            HEId = bHE.OHEId;
            //------------------------------------------------------
            */
            //------------------------------------------------------
            HEPos++;
            if (HEPos < HEIds.Count)
            {
                HEId = HEIds[HEPos];
            }
            //------------------------------------------------------
        }

        public bool IsEnd()
        {
            if (IsInitial)
            {
                return false;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsUseVertexId(UVId));
            /*
            //------------------------------------------------------
            // original
            UseVertex uV = BRep2D.BRep.GetUseVertex(UVId);
            if (HEId == uV.HEId)
            {
                return true;
            }
            return false;
            //------------------------------------------------------
            */
            //------------------------------------------------------
            if (HEIds.Count == 0 ||
                HEPos > HEIds.Count - 1)
            {
                return true;
            }
            return false;
            //------------------------------------------------------
        }

        public uint GetHalfEdgeId()
        {
            return HEId;
        }

        public uint GetUseVertexId()
        {
            return UVId;
        }

        public uint GetLoopId()
        {
            if (!IsValid)
            {
                return 0;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(HEId));
            HalfEdge hE = BRep2D.BRep.GetHalfEdge(HEId);
            uint uLId = hE.ULId;
            if (!BRep2D.BRep.IsUseLoopId(uLId))
            {
                //System.Diagnostics.Debug.Assert(false);
                return 0;
            }
            UseLoop uL = BRep2D.BRep.GetUseLoop(uLId);
            return uL.LId;
        }

        public uint GetEdgeCount()
        {
            if (!IsValid)
            {
                return 0;
            }
            /*
            //------------------------------------------------------
            // original
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsUseVertexId(UVId));
            UseVertex uV = BRep2D.BRep.GetUseVertex(UVId);
            uint hEId0 = uV.HEId;
            uint hEId = hEId0;
            uint iCnt = 0;
            while (true)
            {
                System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(hEId));
                HalfEdge hE = BRep2D.BRep.GetHalfEdge(hEId);
                if (hE.EId == 0)
                {
                    return 0;
                }
                iCnt++;
                uint bHEId = hE.BHEId;
                System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(bHEId));
                HalfEdge bHE = BRep2D.BRep.GetHalfEdge(bHEId);
                hEId = bHE.OHEId;
                if (hEId == hEId0)
                {
                    break;
                }
            }
            return iCnt;
            //------------------------------------------------------
            */
            //------------------------------------------------------
            uint iCnt = (uint)HEIds.Count;
            return iCnt;
            //------------------------------------------------------
        }

        public bool IsParent()
        {
            if (!IsValid)
            {
                return false;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(HEId));
            HalfEdge hE = BRep2D.BRep.GetHalfEdge(HEId);
            uint uLId = hE.ULId;
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsUseLoopId(uLId));
            UseLoop uL = BRep2D.BRep.GetUseLoop(uLId);
            return (uL.ParentULId == uLId);
        }

        public bool IsSameUseLoop(VertexEdgeItr lItr)
        {
            if (!IsValid)
            {
                return false;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(HEId));
            HalfEdge hE0 = BRep2D.BRep.GetHalfEdge(HEId);

            if (!lItr.IsValid)
            {
                return false;
            }
            System.Diagnostics.Debug.Assert(lItr.BRep2D.BRep.IsHalfEdgeId(lItr.HEId));
            HalfEdge hE1 = lItr.BRep2D.BRep.GetHalfEdge(lItr.HEId);
            return (hE0.ULId == hE1.ULId);
        }

        public bool GetBehindEdgeId(out uint eId, out bool isSameDir)
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
            return true;
        }

        public bool GetAheadEdgeId(out uint eId, out bool isSameDir)
        {
            eId = 0;
            isSameDir = true;

            if (!IsValid)
            {
                return false;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(HEId));
            HalfEdge hE = BRep2D.BRep.GetHalfEdge(HEId);
            uint bHEId = hE.BHEId;
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(bHEId));
            HalfEdge bHE = BRep2D.BRep.GetHalfEdge(bHEId);
            eId = bHE.EId;
            isSameDir = !bHE.IsSameDir;
            return true;
        }
    }
}
