using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class ConnectVertexRes
    {
        public uint VId1 { get; set; } = 0;
        public uint VId2 { get; set; } = 0;
        public uint LId { get; set; } = 0;
        public uint AddEId { get; set; } = 0;
        public uint AddLId { get; set; } = 0;
        public bool IsLeftAddL { get; set; } = true;

        public ConnectVertexRes()
        {
        }
    }

    public class BRep2D
    {
        internal BRep BRep { get; private set; } = new BRep();
        internal Dictionary<uint, uint> Loop2UseLoop { get; private set; } = new Dictionary<uint, uint>();
        internal Dictionary<uint, uint> Edge2HalfEdge { get; private set; } = new Dictionary<uint, uint>();

        public BRep2D()
        {

        }

        public void Copy(BRep2D src)
        {
            BRep.Copy(src.BRep);
            Loop2UseLoop = new Dictionary<uint, uint>(src.Loop2UseLoop);
            Edge2HalfEdge = new Dictionary<uint, uint>(src.Edge2HalfEdge);
        }

        public void Clear()
        {
            BRep.Clear();
            Edge2HalfEdge.Clear();
            Loop2UseLoop.Clear();
        }

        public bool IsElementId(CadElementType type, uint id)
        {
            if (type == CadElementType.Vertex)
            {
                return BRep.IsUseVertexId(id);
            }
            else if (type == CadElementType.Edge)
            {
                return Edge2HalfEdge.ContainsKey(id);
            }
            else if (type == CadElementType.Loop)
            {
                return Loop2UseLoop.ContainsKey(id);
            }
            return false;
        }

        public LoopEdgeItr GetLoopEdgeItr(uint lId)
        {
            return new LoopEdgeItr(this, lId);
        }

        public VertexEdgeItr GetVertexEdgeItr(uint vId)
        {
            return new VertexEdgeItr(this, vId);
        }

        public IList<uint> GetElementIds(CadElementType type)
        {
            if (type == CadElementType.Vertex)
            {
                return BRep.GetUseVertexIds();
            }
            IList<uint> res = new List<uint>();
            if (type == CadElementType.Edge)
            {
                foreach (var key in Edge2HalfEdge.Keys)
                {
                    res.Add(key);
                }
            }
            else if (type == CadElementType.Loop)
            {
                foreach (var key in Loop2UseLoop.Keys)
                {
                    res.Add(key);
                }
            }
            return res;
        }

        internal uint GetUseLoopType(uint uLId)
        {
            System.Diagnostics.Debug.Assert(BRep.IsUseLoopId(uLId));
            uint hEId0;
            {
                UseLoop uL = BRep.GetUseLoop(uLId);
                hEId0 = uL.HEId;
                System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId0));
                HalfEdge hE0 = BRep.GetHalfEdge(hEId0);
                if (hE0.FHEId == hEId0)
                {
                    System.Diagnostics.Debug.Assert(hE0.BHEId == hEId0);
                    return 0;
                }
            }
            uint hEId = hEId0;
            while (true)
            {
                System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId));
                HalfEdge hE = BRep.GetHalfEdge(hEId);
                uint fHEId = hE.FHEId;
                {
                    uint oHEId = hE.OHEId;
                    System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(oHEId));
                    HalfEdge oHE = BRep.GetHalfEdge(oHEId);
                    if (oHE.ULId != uLId)
                    {
                        return 2;
                    }
                }
                if (fHEId == hEId0) break;
                hEId = fHEId;
            }
            return 1;
        }

        public bool RemoveVertex(uint vId)
        {
            uint uVId = vId;
            if (!BRep.IsUseVertexId(uVId))
            {
                return false;
            }

            uint nEdgeAroundVtx = 0;
            {
                VertexEdgeItr vItr = new VertexEdgeItr(this, vId);
                nEdgeAroundVtx = vItr.GetEdgeCount();
            }

            if (nEdgeAroundVtx == 0)
            {
                System.Diagnostics.Debug.Assert(BRep.IsUseVertexId(uVId));
                if (!BRep.KVEL(uVId))
                {
                    System.Diagnostics.Debug.Assert(false);
                }
            }
            else if (nEdgeAroundVtx == 2)
            {
                UseVertex uV = BRep.GetUseVertex(uVId);
                uint hEId1 = uV.HEId;
                HalfEdge hE1 = BRep.GetHalfEdge(hEId1);
                uint eId1 = hE1.EId;
                uint hEId2 = hE1.OHEId;
                HalfEdge hE2 = BRep.GetHalfEdge(hEId2);
                { 
                    uint uLId1 = hE1.ULId;
                    uint uLId2 = hE2.ULId;
                    LoopEdgeItr lItr1 = new LoopEdgeItr(this, hEId1, uLId1);
                    LoopEdgeItr lItr2 = new LoopEdgeItr(this, hEId2, uLId2);
                    uint nLV1 = lItr1.GetUseLoopVertexCount();
                    uint nLV2 = lItr2.GetUseLoopVertexCount();
                    System.Diagnostics.Debug.Assert(nLV1 > 1 && nLV2 > 1);
                    if (nLV1 == 2 || nLV2 == 2)
                    {
                        return false;
                    }
                }
                uint uVId2 = hE2.UVId;
                System.Diagnostics.Debug.Assert(BRep.IsUseVertexId(uVId2));
                uint removeHEId1 = hEId1;
                if (!BRep.KVE(removeHEId1))
                {
                    System.Diagnostics.Debug.Assert(false);
                }
                System.Diagnostics.Debug.Assert(BRep.AssertValidUse() == 0);
                Edge2HalfEdge.Remove(eId1);
            }
            System.Diagnostics.Debug.Assert(AssertValid());
            return true;
        }

        public bool AssertValid()
        {
            foreach (var pair in Loop2UseLoop)
            {
                uint lId = pair.Key;
                uint parentULId = pair.Value;
                System.Diagnostics.Debug.Assert(BRep.IsUseLoopId(parentULId));
                {
                    UseLoop uL = BRep.GetUseLoop(parentULId);
                    System.Diagnostics.Debug.Assert(uL.Id == parentULId);
                    System.Diagnostics.Debug.Assert(uL.LId == lId);
                }
                System.Diagnostics.Debug.Assert(GetUseLoopType(parentULId) == 2);
                uint uLId = parentULId;
                while (true)
                {
                    System.Diagnostics.Debug.Assert(BRep.IsUseLoopId(uLId));
                    UseLoop uL = BRep.GetUseLoop(uLId);
                    System.Diagnostics.Debug.Assert(uL.Id == uLId);
                    System.Diagnostics.Debug.Assert(uL.LId == lId);      
                    System.Diagnostics.Debug.Assert(uL.ParentULId == parentULId);
                    uLId = uL.ChildULId;
                    if (uLId == 0)
                    {
                        break;
                    }
                }
            }

            foreach (var pair in Edge2HalfEdge)
            {
                uint hEId = pair.Value;
                System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId));
            }

            IList<uint> uLIds = BRep.UseLoopArray.GetObjectIds();
            for (int i = 0; i < uLIds.Count; i++)
            {
                uint uLId = uLIds[i];
                System.Diagnostics.Debug.Assert(BRep.IsUseLoopId(uLId));
                UseLoop uL = BRep.GetUseLoop(uLId);
                System.Diagnostics.Debug.Assert(uL.Id == uLId);
                uint lId = uL.LId;
                if (lId == 0)
                {
                    System.Diagnostics.Debug.Assert(uL.ParentULId == 0);
                    //2019-03-11 RemoveElement FIX
                    //System.Diagnostics.Debug.Assert(uL.ChildULId == 0);
                    continue;
                }
                System.Diagnostics.Debug.Assert(Loop2UseLoop.ContainsKey(lId));
            }

            IList<uint> hEIds = BRep.HalfEdgeArray.GetObjectIds();
            for (int i = 0; i < hEIds.Count; i++)
            {
                uint hEId = hEIds[i];
                System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId));
                HalfEdge hEdge = BRep.GetHalfEdge(hEId);
                System.Diagnostics.Debug.Assert(hEdge.Id == hEId);

                uint vId1;
                {
                    uint uVId1 = hEdge.UVId;
                    System.Diagnostics.Debug.Assert(BRep.IsUseVertexId(uVId1));
                    UseVertex uV = BRep.GetUseVertex(uVId1);
                    System.Diagnostics.Debug.Assert(uV.Id == uVId1);
                    vId1 = uV.VId;
                }

                uint vId2;
                {
                    uint fHEId = hEdge.FHEId;
                    System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(fHEId));
                    HalfEdge cwEdge = BRep.GetHalfEdge(fHEId);
                    System.Diagnostics.Debug.Assert(cwEdge.Id == fHEId);
                    System.Diagnostics.Debug.Assert(cwEdge.BHEId == hEId);
                    System.Diagnostics.Debug.Assert(cwEdge.ULId == hEdge.ULId);
                    uint uvId2 = cwEdge.UVId;
                    System.Diagnostics.Debug.Assert(BRep.IsUseVertexId(uvId2));
                    UseVertex uV = BRep.GetUseVertex(uvId2);
                    System.Diagnostics.Debug.Assert(uV.Id == uvId2);
                    vId2 = uV.VId;
                }

                bool isSameDir = hEdge.IsSameDir;
                uint eId = hEdge.EId;
                if (eId == 0)
                {
                    System.Diagnostics.Debug.Assert(hEdge.OHEId == hEId);
                    System.Diagnostics.Debug.Assert(hEdge.BHEId == hEId);
                    System.Diagnostics.Debug.Assert(hEdge.FHEId == hEId);
                    continue;
                }

                System.Diagnostics.Debug.Assert(Edge2HalfEdge.ContainsKey(eId));

                uint sVId;
                uint eVId;
                GetEdgeVertexId(eId, out sVId, out eVId);

                if (isSameDir)
                {
                    System.Diagnostics.Debug.Assert(vId1 == sVId);
                    System.Diagnostics.Debug.Assert(vId2 == eVId);
                }
                else
                {
                    System.Diagnostics.Debug.Assert(vId1 == eVId);
                    System.Diagnostics.Debug.Assert(vId2 == sVId);
                }
            }

            IList<uint> uVIds = BRep.UseVertexArray.GetObjectIds();
            for (int i = 0; i < uVIds.Count; i++)
            {
                uint uVId = uVIds[i];
                System.Diagnostics.Debug.Assert(BRep.IsUseVertexId(uVId));
                UseVertex uV = BRep.GetUseVertex(uVId);
                System.Diagnostics.Debug.Assert(uV.Id == uVId);
                uint vId = uV.VId;
                System.Diagnostics.Debug.Assert(uVId == vId);
            }
            return true;
        }

        public uint AddVertexToLoop(uint lId)
        {
            uint uLId = 0;
            {
                if (Loop2UseLoop.ContainsKey(lId))
                {
                    uLId = Loop2UseLoop[lId];
                }
            }

            uint addUVId;
            uint addHEId;
            uint addULId;
            BRep.MVEL(out addUVId, out addHEId, out addULId, uLId);
            BRep.AssertValidUse();
            BRep.SetUseLoopLoopId(addULId, lId);

            uint addVId = addUVId;
            BRep.SetUseVertexVertexId(addUVId, addVId);
            System.Diagnostics.Debug.Assert(AssertValid());
            return addVId;
        }

        public bool GetEdgeLoopId(uint eId, out uint lLId, out uint rLId)
        {
            lLId = 0;
            rLId = 0;

            uint hEId0;
            {
                if (!Edge2HalfEdge.ContainsKey(eId))
                {
                    return false;
                }
                hEId0 = Edge2HalfEdge[eId];
            }
            //!!! 2019-03-11 RemoveElement FIX
            if (!BRep.IsHalfEdgeId(hEId0))
            {
                return false;
            }
            System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId0));
            HalfEdge hE0 = BRep.GetHalfEdge(hEId0);
            uint hEId1 = hE0.OHEId;
            System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId1));
            HalfEdge hE1 = BRep.GetHalfEdge(hEId1);
            uint lULId;
            uint rULId;
            if (hE0.IsSameDir)
            {
                System.Diagnostics.Debug.Assert(!hE1.IsSameDir);
                lULId = hE0.ULId;
                rULId = hE1.ULId;
            }
            else
            {
                System.Diagnostics.Debug.Assert(hE1.IsSameDir);
                lULId = hE1.ULId;
                rULId = hE0.ULId;
            }
            UseLoop lUL = BRep.GetUseLoop(lULId);
            System.Diagnostics.Debug.Assert(BRep.IsUseLoopId(lULId));
            UseLoop rUL = BRep.GetUseLoop(rULId);
            System.Diagnostics.Debug.Assert(BRep.IsUseLoopId(rULId));
            lLId = lUL.LId;
            rLId = rUL.LId;
            return true;
        }

        public uint GetEdgeLoopId(uint eId, bool isLeft)
        {
            uint heId0;
            {
                if (Edge2HalfEdge.ContainsKey(eId))
                {
                    heId0 = Edge2HalfEdge[eId];
                }
                else
                {
                    return 0;
                }
            }
            System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(heId0));
            HalfEdge he0 = BRep.GetHalfEdge(heId0);
            uint heId1 = he0.OHEId;
            System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(heId1));
            HalfEdge he1 = BRep.GetHalfEdge(heId1);
            uint ulId = 0;
            System.Diagnostics.Debug.Assert(he0.IsSameDir != he1.IsSameDir);
            if (he0.IsSameDir == isLeft)
            {
                ulId = he0.ULId;
            }
            else
            {
                ulId = he1.ULId;
            }
            System.Diagnostics.Debug.Assert(BRep.IsUseLoopId(ulId));
            UseLoop ul = BRep.GetUseLoop(ulId);
            return ul.LId;
        }

        public LoopEdgeItr GetSideEdgeLoopEdgeItr(uint eId, bool isLeft)
        {
            System.Diagnostics.Debug.Assert(Edge2HalfEdge.ContainsKey(eId));
            uint hEId1 = Edge2HalfEdge[eId];
            System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId1));
            HalfEdge hE1 = BRep.GetHalfEdge(hEId1);
            if (isLeft)
            {
                uint uLId1 = hE1.ULId;
                return new LoopEdgeItr(this, hEId1, uLId1);
            }
            uint hEId2 = hE1.OHEId;
            System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId2));
            HalfEdge hE2 = BRep.GetHalfEdge(hEId2);
            uint uLId2 = hE2.ULId;
            return new LoopEdgeItr(this, hEId2, uLId2);
        }

        public bool GetEdgeVertexId(uint eId, out uint vId1, out uint vId2)
        {
            vId1 = 0;
            vId2 = 0;

            uint hEId;
            {
                if (!Edge2HalfEdge.ContainsKey(eId))
                {
                    return false;
                }
                hEId = Edge2HalfEdge[eId];
            }
            // 2019-03-11 RemoveElement FIX
            if (!BRep.IsHalfEdgeId(hEId))
            {
                return false;
            }
            System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId));
            HalfEdge hE = BRep.GetHalfEdge(hEId);
            vId1 = hE.UVId;
            uint fHEId = hE.FHEId;
            System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(fHEId));
            HalfEdge fHE = BRep.GetHalfEdge(fHEId);
            vId2 = fHE.UVId;
            return true;
        }

        public uint GetEdgeVertexId(uint eId, bool isRoot)
        {
            uint hEId;
            {
                hEId = Edge2HalfEdge[eId];
            }
            if (!BRep.IsHalfEdgeId(hEId))
            {
                return 0;
            }
            HalfEdge hE = BRep.GetHalfEdge(hEId);
            uint vId1 = hE.UVId;
            uint fHEId = hE.FHEId;
            if (!BRep.IsHalfEdgeId(fHEId))
            {
                return 0;
            }
            HalfEdge fHE = BRep.GetHalfEdge(fHEId);
            uint vId2 = fHE.UVId;
            if (isRoot)
            {
                return vId1;
            }
            return vId2;
        }

        private uint GetFreeKey(Dictionary<uint, uint> dictionary)
        {
            if (dictionary.Count == 0)
            {
                return 1;
            }

            IList<uint> keys = dictionary.Keys.ToList();
            uint maxId = keys.Max();

            for (uint i = 1; i < maxId; i++)
            {
                if (!keys.Contains(i))
                {
                    return i;
                }
            }
            return maxId + 1;
        }

        public uint AddVertexToEdge(uint eId)
        {
            uint hEId;
            {
                if (!Edge2HalfEdge.ContainsKey(eId))
                {
                    return 0;
                }
                hEId = Edge2HalfEdge[eId];
            }
            uint vId1;
            uint vId2;
            if (!GetEdgeVertexId(eId, out vId1, out vId2))
            {
                return 0;
            }

            uint addUVId;
            HalfEdge hE = BRep.GetHalfEdge(hEId);
            System.Diagnostics.Debug.Assert(hE.EId == eId);

            uint oHEId = hE.OHEId;
            System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(oHEId));
            uint hEId1;
            uint hEId2;
            if (hE.IsSameDir)
            {
                hEId1 = hEId;
                hEId2 = oHEId;
            }
            else
            {
                hEId1 = oHEId;
                hEId2 = hEId;
            }
            uint addHEId1;
            uint addHEId2;
            if (!BRep.MVE(out addHEId1, out addHEId2, out addUVId, hEId1))
            {
                System.Diagnostics.Debug.Assert(false);
            }
            BRep.AssertValidUse();
            uint addEId = GetFreeKey(Edge2HalfEdge);
            BRep.SetHalfEdgeEdgeId(hEId1, eId, true);
            BRep.SetHalfEdgeEdgeId(addHEId1, addEId, true);
            BRep.SetHalfEdgeEdgeId(hEId2, addEId, false);
            BRep.SetHalfEdgeEdgeId(addHEId2, eId, false);
            uint addVId = addUVId;
            BRep.SetUseVertexVertexId(addUVId, addVId);
            Edge2HalfEdge.Add(addEId, addHEId1);
            System.Diagnostics.Debug.Assert(AssertValid());
            return addVId;
        }

        public IList<KeyValuePair<uint, bool>> GetEdgesForConnectVertex(VertexEdgeItr vItr1, VertexEdgeItr vItr2)
        {
            IList<KeyValuePair<uint, bool>> eId2Dir = new List<KeyValuePair<uint, bool>>(); ;

            uint uVId1 = vItr1.GetUseVertexId();
            uint uVId2 = vItr2.GetUseVertexId();
            if (!BRep.IsUseVertexId(uVId1))
            {
                return eId2Dir;
            }
            if (!BRep.IsUseVertexId(uVId2))
            {
                return eId2Dir;
            }
            if (uVId1 == uVId2)
            {
                return eId2Dir;
            }
            if (vItr1.GetLoopId() != vItr2.GetLoopId())
            {
                return eId2Dir;
            }

            uint hEId1 = vItr1.GetHalfEdgeId();
            uint uLId1;
            {
                System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId1));
                HalfEdge hE1 = BRep.GetHalfEdge(hEId1);
                uLId1 = hE1.ULId;
                System.Diagnostics.Debug.Assert(hE1.EId != 0);
            }
            uint hEId2 = vItr2.GetHalfEdgeId();
            uint uLId2;
            {
                System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId2));
                HalfEdge hE2 = BRep.GetHalfEdge(hEId2);
                uLId2 = hE2.ULId;
                System.Diagnostics.Debug.Assert(hE2.EId != 0);
            }

            if (uLId1 != uLId2)
            {
                return eId2Dir;
            }
            uint hEId = hEId2;
            while (true)
            {
                HalfEdge hE = BRep.GetHalfEdge(hEId);
                eId2Dir.Add(new KeyValuePair<uint, bool>(hE.EId, hE.IsSameDir));
                hEId = hE.FHEId;
                if (hEId == hEId1)
                {
                    break;
                }
            }
            return eId2Dir;
        }

        public ConnectVertexRes ConnectVertex(VertexEdgeItr vItr1, VertexEdgeItr vItr2, bool isLeftAddL)
        {
            ConnectVertexRes res = new ConnectVertexRes();

            uint uVId1 = vItr1.GetUseVertexId();
            uint uVId2 = vItr2.GetUseVertexId();
            if (!BRep.IsUseVertexId(uVId1))
            {
                return res;
            }
            if (!BRep.IsUseVertexId(uVId2))
            {
                return res;
            }
            if (uVId1 == uVId2)
            {
                return res;
            }
            res.VId1 = BRep.GetUseVertex(uVId1).VId;
            res.VId2 = BRep.GetUseVertex(uVId2).VId;

            if (vItr1.GetLoopId() != vItr2.GetLoopId())
            {
                return res;
            }
            uint lId = vItr1.GetLoopId();
            res.LId = lId;

            uint hEId1 = vItr1.GetHalfEdgeId();
            uint hEId2 = vItr2.GetHalfEdgeId();

            uint uLId1;
            bool isFloat1;
            {
                System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId1));
                HalfEdge hE1 = BRep.GetHalfEdge(hEId1);
                uLId1 = hE1.ULId;
                System.Diagnostics.Debug.Assert(BRep.IsUseLoopId(uLId1));
                isFloat1 = (hE1.EId == 0);
            }
            uint uLId2;
            bool isFloat2;
            {
                System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId2));
                HalfEdge hE2 = BRep.GetHalfEdge(hEId2);
                uLId2 = hE2.ULId;
                System.Diagnostics.Debug.Assert(BRep.IsUseLoopId(uLId2));
                isFloat2 = (hE2.EId == 0);
            }

            if (uLId1 != uLId2)
            {
                if (isFloat1 && !isFloat2)
                {
                    uint addHEId1;
                    if (!BRep.MEKLOneFloatingVertex(out addHEId1, hEId2, hEId1))
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    uint addEId = GetFreeKey(Edge2HalfEdge);
                    BRep.SetHalfEdgeEdgeId(addHEId1, addEId, false);
                    BRep.SetHalfEdgeEdgeId(hEId1, addEId, true);
                    Edge2HalfEdge.Add(addEId, hEId1);
                    System.Diagnostics.Debug.Assert(AssertValid());
                    res.AddEId = addEId;
                    return res;
                }
                else if (isFloat2 && !isFloat1)
                {
                    uint addHEId1;
                    if (!BRep.MEKLOneFloatingVertex(out addHEId1, hEId1, hEId2))
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    uint addEId = GetFreeKey(Edge2HalfEdge);
                    BRep.SetHalfEdgeEdgeId(addHEId1, addEId, true);
                    BRep.SetHalfEdgeEdgeId(hEId2, addEId, false);
                    Edge2HalfEdge.Add(addEId, addHEId1);
                    System.Diagnostics.Debug.Assert(AssertValid());
                    res.AddEId = addEId;
                    return res;
                }
                else if (isFloat1 && isFloat2)
                {
                    if (!BRep.MEKLTwoFloatingVertex(hEId1, hEId2))
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    uint addEId = GetFreeKey(Edge2HalfEdge);
                    BRep.SetHalfEdgeEdgeId(hEId1, addEId, true);
                    BRep.SetHalfEdgeEdgeId(hEId2, addEId, false);
                    Edge2HalfEdge.Add(addEId, hEId1);
                    System.Diagnostics.Debug.Assert(AssertValid());
                    res.AddEId = addEId;
                    return res;
                }
                else
                {
                    uint addHEId1;
                    uint addHEId2;
                    if (!BRep.MEKL(out addHEId1, out addHEId2, hEId1, hEId2))
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    uint addEId = GetFreeKey(Edge2HalfEdge);
                    BRep.SetHalfEdgeEdgeId(addHEId1, addEId, true);
                    BRep.SetHalfEdgeEdgeId(addHEId2, addEId, false);
                    Edge2HalfEdge.Add(addEId, addHEId1);
                    if (lId != 0)
                    {
                        uint parentULId1 = BRep.GetUseLoop(uLId1).ParentULId;
                        System.Diagnostics.Debug.Assert(parentULId1 != 0);
                        System.Diagnostics.Debug.Assert(Loop2UseLoop.ContainsKey(lId));
                        Loop2UseLoop[lId] = parentULId1;
                    }
                    else
                    {

                    }
                    System.Diagnostics.Debug.Assert(AssertValid());
                    res.AddEId = addEId;
                    return res;
                }
            }
            else if (uLId1 == uLId2)
            {
                uint addHEId1;
                uint addHEId2;
                uint addULId;
                if (!BRep.MEL(out addHEId1, out addHEId2, out addULId, hEId1, hEId2))
                {
                    System.Diagnostics.Debug.Assert(false);
                }
                System.Diagnostics.Debug.Assert(BRep.AssertValidUse() == 0);
                uint addEId = GetFreeKey(Edge2HalfEdge);
                BRep.SetHalfEdgeEdgeId(addHEId1, addEId, true);
                BRep.SetHalfEdgeEdgeId(addHEId2, addEId, false);

                Edge2HalfEdge.Add(addEId, addHEId1);
                bool iflag = true;
                if (lId == 0)
                { 
                    iflag = !isLeftAddL;
                }
                else
                {
                    uint parentULId =Loop2UseLoop[lId];
                    if (uLId1 != parentULId)
                    {
                        if (isLeftAddL)
                        {
                            BRep.SwapUseLoop(addULId, uLId1);
                            BRep.SetUseLoopLoopId(addULId, lId);
                            BRep.AssertValidUse();
                            iflag = false;
                        }
                        else
                        {

                        }
                    }
                    else
                    {

                    }
                }

                uint addLId = GetFreeKey(Loop2UseLoop);
                if (iflag)
                {
                    Loop2UseLoop.Add(addLId, addULId);
                    BRep.SetUseLoopLoopId(addULId, addLId);
                }
                else
                {
                    Loop2UseLoop.Add(addLId, uLId1);
                    BRep.SetUseLoopLoopId(uLId1, addLId);
                }
                System.Diagnostics.Debug.Assert(AssertValid());
                res.AddEId = addEId;
                res.AddLId = addLId;
                res.IsLeftAddL = !iflag;
                return res;
            }
            System.Diagnostics.Debug.Assert(false);
            return res;
        }

        public IList<KeyValuePair<uint, bool>> GetEdgesForRemoveEdge(uint eId)
        {
            IList<KeyValuePair<uint, bool>> id2Dir = new List<KeyValuePair<uint, bool>>();
            uint hEId1;
            {
                if (!Edge2HalfEdge.ContainsKey(eId))
                {
                    return id2Dir;
                }
                hEId1 = Edge2HalfEdge[eId];
            }
            System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId1));
            HalfEdge hE1 = BRep.GetHalfEdge(hEId1);
            uint hEId2 = hE1.OHEId;
            if (hEId2 == hEId1)
            {
                return id2Dir;
            }
            System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId2));
            HalfEdge hE2 = BRep.GetHalfEdge(hEId2);
            if (hE1.ULId != hE2.ULId)
            {
                return id2Dir;
            }

            uint hEId = hE2.FHEId;
            while (true)
            {
                System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId));
                HalfEdge hE = BRep.GetHalfEdge(hEId);
                id2Dir.Add(new KeyValuePair<uint, bool>(hE.EId, hE.IsSameDir));
                uint fHEId = hE.FHEId;
                System.Diagnostics.Debug.Assert(fHEId != hEId2);
                if (fHEId == hEId1) break;
                hEId = fHEId;
            }
            return id2Dir;
        }

        public bool RemoveEdge(uint eId, bool isDelCP)
        {
            uint hEId1 = 0;
            {
                if (!Edge2HalfEdge.ContainsKey(eId))
                {
                    return false;
                }
                hEId1 = Edge2HalfEdge[eId];
            }
            uint hEId2;
            {
                System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId1));
                HalfEdge hE1 = BRep.GetHalfEdge(hEId1);
                hEId2 = hE1.OHEId;
            }
            uint uLId1;
            uint lId1;
            bool isEdgeVertex2 = false;
            {
                HalfEdge hE1 = BRep.GetHalfEdge(hEId1);
                System.Diagnostics.Debug.Assert(hE1.OHEId == hEId2);
                uLId1 = hE1.ULId;
                UseLoop uL1 = BRep.GetUseLoop(uLId1);
                lId1 = uL1.LId;
                if (hE1.FHEId == hEId2)
                {
                    isEdgeVertex2 = true;
                }
            }

            uint uLId2;
            uint lId2;
            bool isEdgeVertex1 = false;
            {
                HalfEdge hE2 = BRep.GetHalfEdge(hEId2);
                System.Diagnostics.Debug.Assert(hE2.OHEId == hEId1);
                uLId2 = hE2.ULId;
                UseLoop uL2 = BRep.GetUseLoop(uLId2);
                lId2 = uL2.LId;
                if (hE2.FHEId == hEId1)
                {
                    isEdgeVertex1 = true;
                }
            }
            if (lId1 != lId2)
            {
                System.Diagnostics.Debug.Assert(uLId1 != uLId2);
                System.Diagnostics.Debug.Assert(!isEdgeVertex1 && !isEdgeVertex2);
                uint newLId;
                if (lId1 != 0)
                {
                    if (lId2 == 0)
                    {
                        newLId = 0;
                    }
                    else
                    {
                        newLId = lId1;
                    }
                }
                else
                {
                    newLId = 0;
                }
                {
                    uint uLId = uLId1;
                    UseLoop uL = BRep.GetUseLoop(uLId);
                    uLId = uL.ParentULId;
                    while (true)
                    {
                        if (uLId == 0)
                        {
                            break;
                        }
                        BRep.SetUseLoopLoopId(uLId, newLId);
                        System.Diagnostics.Debug.Assert(BRep.IsUseLoopId(uLId));
                        UseLoop tmpUL = BRep.GetUseLoop(uLId);
                        uLId = tmpUL.ChildULId;
                    }
                }
                {
                    uint uLId = uLId2;
                    UseLoop uL = BRep.GetUseLoop(uLId);
                    uLId = uL.ParentULId;
                    while (true)
                    {
                        if (uLId == 0)
                        {
                            break;
                        }
                        BRep.SetUseLoopLoopId(uLId, newLId);
                        System.Diagnostics.Debug.Assert(BRep.IsUseLoopId(uLId));
                        UseLoop tmpUL = BRep.GetUseLoop(uLId);
                        uLId = tmpUL.ChildULId;
                    }
                }
                if (lId1 != 0)
                {
                    if (!BRep.KEL(hEId1))
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    if (lId2 == 0)
                    {
                        Loop2UseLoop.Remove(lId1);
                    }
                    else
                    {
                        Loop2UseLoop.Remove(lId1);
                        Loop2UseLoop.Remove(lId2);
                        System.Diagnostics.Debug.Assert(BRep.IsUseLoopId(uLId1));
                        UseLoop uL = BRep.GetUseLoop(uLId1);
                        System.Diagnostics.Debug.Assert(BRep.IsUseLoopId(uL.ParentULId));
                        Loop2UseLoop.Add(lId1, uL.ParentULId);
                    }
                }
                else
                {
                    System.Diagnostics.Debug.Assert(lId2 != 0);
                    if (!BRep.KEL(hEId2))
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    Loop2UseLoop.Remove(lId2);
                }
                Edge2HalfEdge.Remove(eId);
                System.Diagnostics.Debug.Assert(AssertValid());
                return true;
            }
            if (uLId1 == uLId2)
            {
                System.Diagnostics.Debug.Assert(lId1 == lId2);
                if (!isEdgeVertex1 && !isEdgeVertex2)
                {
                    if (lId1 != 0)
                    {
                        uint addULId;
                        if (!BRep.KEML(out addULId, hEId1))
                        {
                            System.Diagnostics.Debug.Assert(false);
                        }
                        if (isDelCP)
                        {
                            System.Diagnostics.Debug.Assert(GetUseLoopType(addULId) == 2);
                            if (!BRep.SwapUseLoop(addULId, uLId1))
                            {
                                System.Diagnostics.Debug.Assert(false);
                            }
                            Loop2UseLoop[lId1] = addULId;
                        }
                        BRep.SetUseLoopLoopId(addULId, lId1);
                    }
                    else
                    {
                        uint addULId;
                        if (!BRep.KEML(out addULId, hEId1))
                        {
                            System.Diagnostics.Debug.Assert(false);
                        }
                        BRep.SetUseLoopLoopId(addULId, 0);
                    }
                }
                else if (isEdgeVertex1 && isEdgeVertex2)
                {
                    uint addULId;
                    if (!BRep.KEMLTwoFloatingVertex(out addULId, hEId1))
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    BRep.SetHalfEdgeEdgeId(hEId1, 0, true);
                    BRep.SetHalfEdgeEdgeId(hEId2, 0, true);
                    BRep.SetUseLoopLoopId(addULId, lId1);
                }
                else if (isEdgeVertex1 && !isEdgeVertex2)
                {
                    uint id_ul_add;
                    if (!BRep.KEMLOneFloatingVertex(out id_ul_add, hEId1))
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    BRep.SetUseLoopLoopId(id_ul_add, lId1);
                }
                else if (isEdgeVertex2 && !isEdgeVertex1)
                {
                    uint id_ul_add;
                    if (!BRep.KEMLOneFloatingVertex(out id_ul_add, hEId2))
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    BRep.SetUseLoopLoopId(id_ul_add, lId1);
                }
                Edge2HalfEdge.Remove(eId);
                System.Diagnostics.Debug.Assert(AssertValid());
                return true;
            }
            return false;
        }

        public bool SwapLoopEdgeItr(LoopEdgeItr lItr, uint toLId)
        {
            uint toParentULId;
            {
                toParentULId = Loop2UseLoop[toLId];
                UseLoop uL = BRep.GetUseLoop(toParentULId);
                System.Diagnostics.Debug.Assert(uL.ParentULId == toParentULId);
            }
            uint fromULId = lItr.GetUseLoopId();
            {
                System.Diagnostics.Debug.Assert(BRep.IsUseLoopId(fromULId));
                UseLoop uL = BRep.GetUseLoop(fromULId);
                System.Diagnostics.Debug.Assert(uL.ParentULId != fromULId);
            }
            BRep.MoveUseLoop(fromULId, toParentULId);
            BRep.SetUseLoopLoopId(fromULId, toLId);
            return true;
        }

        public bool MakeHoleFromLoop(uint lId)
        {
            uint uLId1 = 0;
            {
                if (!Loop2UseLoop.ContainsKey(lId))
                {
                    return false;
                }
                uLId1 = Loop2UseLoop[lId];
            }
            System.Diagnostics.Debug.Assert(BRep.IsUseLoopId(uLId1));
            BRep.SetUseLoopLoopId(uLId1, 0);
            Loop2UseLoop.Remove(lId);
            return true;
        }

        public bool Serialize(Serializer arch)
        {
            if (arch.IsLoading)
            {   // 読み込み時の処理
                Clear();

                string className;
                string[] values;
                
                className = arch.ReadDepthClassName();
                if (className != "BRep2D")
                {
                    return true;
                }
                {
                    int ne;
                    values = arch.GetValues();
                    ne = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(ne >= 0);
                    for (uint ie = 0; ie < ne; ie++)
                    {
                        int itmp;
                        int eId;
                        int hEId;
                        values = arch.GetValues();
                        itmp = int.Parse(values[0]);
                        eId = int.Parse(values[1]);
                        hEId = int.Parse(values[2]);
                        Edge2HalfEdge.Add((uint)eId, (uint)hEId);
                    }
                }
                {
                    int nl;
                    values = arch.GetValues();
                    nl = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(nl >= 0);
                    for (uint il = 0; il < nl; il++)
                    {
                        int itmp;
                        int lId;
                        int uLId;
                        values = arch.GetValues();
                        itmp = int.Parse(values[0]);
                        lId = int.Parse(values[1]);
                        uLId = int.Parse(values[2]);
                        Loop2UseLoop.Add((uint)lId, (uint)uLId);
                    }
                }
                ////////////////
                int nuv;
                int nhe;
                int nul;
                {
                    values = arch.GetValues();
                    nuv = int.Parse(values[0]);
                    nhe = int.Parse(values[1]);
                    nul = int.Parse(values[2]);
                    System.Diagnostics.Debug.Assert(nuv > 0);
                    System.Diagnostics.Debug.Assert(nhe > 0);
                    System.Diagnostics.Debug.Assert(nul > 0);
                    /*
                    m_BRep.m_UseVertexSet.Reserve(nuv * 2);
                    m_BRep.m_HalfEdgeSet.Reserve(nhe * 2);
                    m_BRep.m_UseLoopSet.Reserve(nul * 2);
                    */
                }
                ////////////////////////////////////////////////
                arch.ShiftDepth(true);
                for (int iuv = 0; iuv < nuv; iuv++)
                {
                    className = arch.ReadDepthClassName();
                    System.Diagnostics.Debug.Assert(className == "CUseVertex");

                    int id;
                    values = arch.GetValues();
                    id = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(id > 0);

                    int vId;
                    values = arch.GetValues();
                    vId = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(vId > 0);

                    int hEId;
                    values = arch.GetValues();
                    hEId = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(hEId > 0);

                    UseVertex uv = new UseVertex((uint)id, (uint)hEId);
                    uint tmpId = BRep.UseVertexArray.AddObject(uv);
                    System.Diagnostics.Debug.Assert(tmpId == id);
                    BRep.SetUseVertexVertexId((uint)id, (uint)vId);
                }
                for (int ihe = 0; ihe < nhe; ihe++)
                {
                    className = arch.ReadDepthClassName();
                    System.Diagnostics.Debug.Assert(className == "CHalfEdge");

                    int id;
                    values = arch.GetValues();
                    id = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(id > 0);

                    int eId;
                    values = arch.GetValues();
                    eId = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(eId >= 0);

                    int iIsSameDir;
                    values = arch.GetValues();
                    iIsSameDir = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(iIsSameDir >= 0);

                    int uvId;
                    values = arch.GetValues();
                    uvId = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(uvId > 0);

                    int fHEId;
                    int ccwHEId;
                    int oHEId;
                    values = arch.GetValues();
                    fHEId = int.Parse(values[0]);
                    ccwHEId = int.Parse(values[1]);
                    oHEId = int.Parse(values[2]);
                    System.Diagnostics.Debug.Assert(
                        fHEId > 0 && ccwHEId > 0 && oHEId > 0);

                    int id_ul;
                    values = arch.GetValues();
                    id_ul = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(id_ul > 0);

                    bool is_same_dir = (iIsSameDir != 0);
                    HalfEdge he = new HalfEdge(
                        (uint)id, (uint)uvId, (uint)fHEId, (uint)ccwHEId, (uint)oHEId, (uint)id_ul);
                    uint tmp_id = BRep.HalfEdgeArray.AddObject(he);
                    System.Diagnostics.Debug.Assert(tmp_id == id);
                    BRep.SetHalfEdgeEdgeId((uint)id, (uint)eId, is_same_dir);
                }
                for (int iul = 0; iul < nul; iul++)
                {
                    className = arch.ReadDepthClassName();
                    System.Diagnostics.Debug.Assert(className == "CUseLoop");

                    int id;
                    values = arch.GetValues();
                    id = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(id > 0);

                    int lId;
                    values = arch.GetValues();
                    lId = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(lId >= 0);

                    int hEId;
                    values = arch.GetValues();
                    hEId = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(hEId > 0);

                    int cULId;
                    values = arch.GetValues();
                    cULId = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(cULId >= 0);

                    int pULId;
                    values = arch.GetValues();
                    pULId = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(pULId >= 0);

                    UseLoop ul = new UseLoop((uint)id, (uint)hEId, (uint)cULId, (uint)pULId);
                    uint tmpId = BRep.UseLoopArray.AddObject(ul);
                    System.Diagnostics.Debug.Assert(tmpId == id);
                    BRep.SetUseLoopLoopId((uint)id, (uint)lId);
                }
                arch.ShiftDepth(false);
                AssertValid();
                return true;
            }
            else
            {
                // 書き込み時の処理
                string line;

                // クラスの名前の指定，サイズの指定
                arch.WriteDepthClassName("BRep2D");
                {
                    line = string.Format("{0}", Edge2HalfEdge.Count);
                    arch.WriteLine(line);

                    uint icnt = 0;
                    foreach (var pair in Edge2HalfEdge)
                    {
                        line = string.Format("{0} {1} {2}", icnt, pair.Key, pair.Value);
                        arch.WriteLine(line);
                        icnt++;
                    }
                }
                {
                    line = string.Format("{0}", Loop2UseLoop.Count);
                    arch.WriteLine(line);

                    uint icnt = 0;
                    foreach (var pair in Loop2UseLoop)
                    {
                        line = string.Format("{0} {1} {2}", icnt, pair.Key, pair.Value);
                        arch.WriteLine(line);
                        icnt++;
                    }
                }

                line = string.Format("{0} {1} {2}", 
                    BRep.UseVertexArray.GetObjectIds().Count,
                    BRep.HalfEdgeArray.GetObjectIds().Count,
                    BRep.UseLoopArray.GetObjectIds().Count);
                arch.WriteLine(line);

                arch.ShiftDepth(true);
                ////////////////
                // UseVertexの出力
                {
                    IList<uint> ids = BRep.UseVertexArray.GetObjectIds();
                    foreach (uint uvId in ids)
                    {
                        System.Diagnostics.Debug.Assert(BRep.IsUseVertexId(uvId));
                        UseVertex uv = BRep.GetUseVertex(uvId);
                        System.Diagnostics.Debug.Assert(uv.Id == uvId);

                        arch.WriteDepthClassName("CUseVertex");

                        line = string.Format("{0}", uvId);
                        arch.WriteLine(line);

                        line = string.Format("{0}", uv.VId);
                        arch.WriteLine(line);

                        line = string.Format("{0}", uv.HEId);
                        arch.WriteLine(line);
                    }
                }
                // HalfEdgeの出力
                {
                    IList<uint> ids = BRep.HalfEdgeArray.GetObjectIds();
                    foreach (uint hEId in ids)
                    {
                        System.Diagnostics.Debug.Assert(BRep.IsHalfEdgeId(hEId));
                        HalfEdge he = BRep.GetHalfEdge(hEId);
                        System.Diagnostics.Debug.Assert(he.Id == hEId);

                        arch.WriteDepthClassName("CHalfEdge");

                        line = string.Format("{0}", hEId);
                        arch.WriteLine(line);

                        line = string.Format("{0}", he.EId);
                        arch.WriteLine(line);

                        line = string.Format("{0}", (he.IsSameDir ? 1 : 0));
                        arch.WriteLine(line);

                        line = string.Format("{0}", he.UVId);
                        arch.WriteLine(line);

                        line = string.Format("{0} {1} {2}", he.FHEId, he.BHEId, he.OHEId);
                        arch.WriteLine(line);

                        line = string.Format("{0}", he.ULId);
                        arch.WriteLine(line);
                    }
                }
                // UseLoopの出力
                {
                    IList<uint> ids = BRep.UseLoopArray.GetObjectIds();
                    foreach (uint uLId in ids)
                    {
                        System.Diagnostics.Debug.Assert(BRep.IsUseLoopId(uLId));
                        UseLoop ul = BRep.GetUseLoop(uLId);
                        System.Diagnostics.Debug.Assert(ul.Id == uLId);

                        arch.WriteDepthClassName("CUseLoop");

                        line = string.Format("{0}", uLId);
                        arch.WriteLine(line);

                        line = string.Format("{0}", ul.LId);
                        arch.WriteLine(line);

                        line = string.Format("{0}", ul.HEId);
                        arch.WriteLine(line);

                        line = string.Format("{0}", ul.ChildULId);
                        arch.WriteLine(line);

                        line = string.Format("{0}", ul.ParentULId);
                        arch.WriteLine(line);
                    }
                }
                arch.ShiftDepth(false);
            }
            return true;
        }

    }
}
