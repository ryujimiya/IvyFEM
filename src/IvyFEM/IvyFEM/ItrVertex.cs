using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class ItrVertex
    {
        private bool IsValid = false;
        private bool IsInitial = false;
        private uint UVId = 0;
        private uint HEId = 0;
        private BRep2D BRep2D = new BRep2D();

        public ItrVertex(BRep2D bRep2D, uint vId)
        {
            BRep2D = bRep2D;
            IsValid = false;
            UVId = vId;
            UseVertex uV = BRep2D.BRep.GetUseVertex(UVId);
            System.Diagnostics.Debug.Assert(uV.Id == UVId);
            HEId = uV.HEId;
            IsInitial = true;
            IsValid = true;
        }

        public ItrVertex(ItrVertex src)
        {
            IsValid = src.IsValid;
            IsInitial = src.IsInitial;
            UVId = src.UVId;
            HEId = src.HEId;
            BRep2D.Copy(src.BRep2D);
        }

        public static ItrVertex operator ++(ItrVertex src)
        {
            ItrVertex dest = new ItrVertex(src);
            if (!src.IsValid)
            {
                return dest;
            }
            dest.IsInitial = false;
            System.Diagnostics.Debug.Assert(dest.BRep2D.BRep.IsHalfEdgeId(dest.HEId));
            HalfEdge hE = dest.BRep2D.BRep.GetHalfEdge(dest.HEId);
            uint bHEId = hE.BHEId;
            System.Diagnostics.Debug.Assert(dest.BRep2D.BRep.IsHalfEdgeId(bHEId));
            HalfEdge bHE = dest.BRep2D.BRep.GetHalfEdge(bHEId);
            dest.HEId = bHE.OHEId;
            return dest;
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
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsUseLoopId(uLId));
            UseLoop uL = BRep2D.BRep.GetUseLoop(uLId);
            return uL.LId;
        }

        public void Begin()
        {
            if (!IsValid)
            {
                return;
            }
            IsInitial = true;
            UseVertex uV = BRep2D.BRep.GetUseVertex(UVId);
            System.Diagnostics.Debug.Assert(uV.Id == UVId);
            HEId = uV.HEId;
        }

        public bool IsEnd()
        {
            if (IsInitial)
            {
                return false;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsUseVertexId(UVId));
            UseVertex uV = BRep2D.BRep.GetUseVertex(UVId);
            if (HEId == uV.HEId)
            {
                return true;
            }
            return false;
        }

        public uint CountEdge()
        {
            if (!IsValid)
            {
                return 0;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsUseVertexId(UVId));
            UseVertex uV = BRep2D.BRep.GetUseVertex(UVId);
            uint hEId0 = uV.HEId;
            uint hEId = hEId0;
            uint iCnt = 0;
            for (;;)
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

        public bool IsSameUseLoop(ItrVertex itrl)
        {
            if (!IsValid)
            {
                return false;
            }
            System.Diagnostics.Debug.Assert(BRep2D.BRep.IsHalfEdgeId(HEId));
            HalfEdge hE0 = BRep2D.BRep.GetHalfEdge(HEId);

            if (!itrl.IsValid)
            {
                return false;
            }
            System.Diagnostics.Debug.Assert(itrl.BRep2D.BRep.IsHalfEdgeId(itrl.HEId));
            HalfEdge hE1 = itrl.BRep2D.BRep.GetHalfEdge(itrl.HEId);
            return (hE0.ULId == hE1.ULId);
        }

    }
}
