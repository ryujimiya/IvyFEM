using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class BRep
    {
        public IdObjectArray<UseLoop> UseLoopArray { get; private set; } = new IdObjectArray<UseLoop>();
        public IdObjectArray<HalfEdge> HalfEdgeArray { get; private set; } = new IdObjectArray<HalfEdge>();
        public IdObjectArray<UseVertex> UseVertexArray { get; private set; } = new IdObjectArray<UseVertex>();

        public BRep()
        {

        }

        public void Copy(BRep src)
        {
            var ulArray = UseLoopArray;
            var srcUlArray = src.UseLoopArray;
            IList<uint> srcUlIds = srcUlArray.GetObjectIds();
            ulArray.Clear();
            foreach (uint id in srcUlIds)
            {
                var obj = srcUlArray.GetObject(id);
                var tmpObj = new UseLoop(obj);
                ulArray.AddObject(obj);
            }

            var heArray = HalfEdgeArray;
            var srcHeArray = src.HalfEdgeArray;
            IList<uint> srcHeIds = srcHeArray.GetObjectIds();
            heArray.Clear();
            foreach (uint id in srcHeIds)
            {
                var obj = srcHeArray.GetObject(id);
                var tmpObj = new HalfEdge(obj);
                heArray.AddObject(obj);
            }

            var uvArray = UseVertexArray;
            var srcUvArray = src.UseVertexArray;
            IList<uint> srcUvIds = srcUvArray.GetObjectIds();
            uvArray.Clear();
            foreach (uint id in srcUvIds)
            {
                var obj = srcUvArray.GetObject(id);
                var tmpObj = new UseVertex(obj);
                uvArray.AddObject(obj);
            }
        }

        public void Clear()
        {
            UseLoopArray.Clear();
            HalfEdgeArray.Clear();
            UseVertexArray.Clear();
        }

        public bool IsUseLoopId(uint uLId)
        {
            return UseLoopArray.IsObjectId(uLId);
        }

        public bool IsHalfEdgeId(uint hEId)
        {
            return HalfEdgeArray.IsObjectId(hEId);
        }

        public bool IsUseVertexId(uint uVId)
        {
            return UseVertexArray.IsObjectId(uVId);
        }

        public IList<uint> GetUseVertexIds()
        {
            return UseVertexArray.GetObjectIds();
        }

        public UseLoop GetUseLoop(uint uLId)
        {
            System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId));
            if (!UseLoopArray.IsObjectId(uLId))
            {
                throw new InvalidOperationException("object id is invalid");
            }
            return UseLoopArray.GetObject(uLId);
        }

        public UseVertex GetUseVertex(uint uVId)
        {
            System.Diagnostics.Debug.Assert(UseVertexArray.IsObjectId(uVId));
            if (!UseVertexArray.IsObjectId(uVId))
            {
                throw new InvalidOperationException("object id is invalid");
            }
            return UseVertexArray.GetObject(uVId);
        }

        public HalfEdge GetHalfEdge(uint hEId)
        {
            System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId));
            if (!HalfEdgeArray.IsObjectId(hEId))
            {
                throw new InvalidOperationException("object id is invalid");
            }
            return HalfEdgeArray.GetObject(hEId);
        }

        public bool SetUseLoopLoopId(uint uLId, uint lId)
        {
            if (!UseLoopArray.IsObjectId(uLId))
            {
                return false;
            }
            UseLoop uL = UseLoopArray.GetObject(uLId);
            uL.LId = lId;
            if (uL.ParentULId == 0 || uL.ParentULId == uLId)
            {
                if (lId == 0)
                {
                    uL.ParentULId = 0;
                }
                else
                {
                    uL.ParentULId = uLId;
                }
            }
            return true;
        }

        public bool SetUseVertexVertexId(uint uVId, uint vId)
        {
            if (!UseVertexArray.IsObjectId(uVId))
            {
                return false;
            }
            UseVertex uV = UseVertexArray.GetObject(uVId);
            uV.VId = vId;
            return true;
        }

        public bool SetHalfEdgeEdgeId(uint hEId, uint eId, bool isSameDir)
        {
            if (!HalfEdgeArray.IsObjectId(hEId))
            {
                return false;
            }
            HalfEdge hE = HalfEdgeArray.GetObject(hEId);
            hE.EId = eId;
            hE.IsSameDir = isSameDir;
            return true;
        }

        public IList<uint> FindHalfEdgeByEdge(uint eId)
        {
            IList<uint> res = new List<uint>();
            IList<uint> hEIds = HalfEdgeArray.GetObjectIds();
            for (int i = 0; i < hEIds.Count; i++)
            {
                uint hEId = hEIds[i];
                HalfEdge hE = HalfEdgeArray.GetObject(hEId);
                if (hE.EId != eId)
                {
                    continue;
                }
                res.Add(hEId);
            }
            return res;
        }

        public IList<uint> FindHalfEdgeByVertex(uint vId)
        {
            IList<uint> res = new List<uint>();
            IList<uint> uVIds = UseVertexArray.GetObjectIds();
            for (int i = 0; i < uVIds.Count; i++)
            {
                uint uVId = uVIds[i];
                UseVertex uV = UseVertexArray.GetObject(uVId);
                if (uV.VId != vId)
                {
                    continue;
                }
                uint hEId = uV.HEId;
                uint hEId0 = hEId;
                res.Add(hEId);
                while (true)
                {
                    HalfEdge hE = HalfEdgeArray.GetObject(hEId);
                    uint oHEId = hE.OHEId;
                    HalfEdge oHE = HalfEdgeArray.GetObject(oHEId);
                    uint nextHEId = oHE.FHEId;
                    if (nextHEId == hEId0)
                    {
                        break;
                    }
                    hEId = nextHEId;
                    res.Add(hEId);
                }
            }
            return res;
        }

        //-----------------------------------------------------
        public IList<uint> FindHalfEdgeByUseVertex(uint uVId)
        {
            IList<uint> res = new List<uint>();
            IList<uint> hEIds = HalfEdgeArray.GetObjectIds();
            foreach (uint hEId in hEIds)
            {
                HalfEdge hE = HalfEdgeArray.GetObject(hEId);
                if (hE.UVId == uVId)
                {
                    res.Add(hEId);
                }
            }
            return res;
        }
        //-----------------------------------------------------

        /// <summary>
        /// heidの起点を消去して２つの辺を１つにする
        /// </summary>
        /// <param name="hEId1"></param>
        /// <returns></returns>
        public bool KVE(uint hEId1)
        {
            uint uVid1;
            uint hEId2;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId1));
                if (!HalfEdgeArray.IsObjectId(hEId1))
                {
                    return false;
                }
                HalfEdge hE1 = HalfEdgeArray.GetObject(hEId1);
                System.Diagnostics.Debug.Assert(hE1.Id == hEId1);
                hEId2 = hE1.OHEId;
                uVid1 = hE1.UVId;
            }
            uint uVId2;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId2));
                HalfEdge hE2 = HalfEdgeArray.GetObject(hEId2);
                System.Diagnostics.Debug.Assert(hE2.Id == hEId2);
                System.Diagnostics.Debug.Assert(hE2.OHEId == hEId1);
                uVId2 = hE2.UVId;
            }

            uint bHEId1;
            uint fHEId1;
            uint uLId1;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId1));
                HalfEdge hE1 = HalfEdgeArray.GetObject(hEId1);
                System.Diagnostics.Debug.Assert(hE1.Id == hEId1);
                System.Diagnostics.Debug.Assert(hE1.UVId == uVid1);
                uLId1 = hE1.ULId;
                bHEId1 = hE1.BHEId;
                fHEId1 = hE1.FHEId;
                {
                    System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId1));
                    HalfEdge bHE1 = HalfEdgeArray.GetObject(bHEId1);
                    System.Diagnostics.Debug.Assert(bHE1.Id == bHEId1);
                    System.Diagnostics.Debug.Assert(bHE1.FHEId == hEId1);
                    System.Diagnostics.Debug.Assert(bHE1.ULId == uLId1);
                }
                {
                    System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(fHEId1));
                    HalfEdge fHE1 = HalfEdgeArray.GetObject(fHEId1);
                    System.Diagnostics.Debug.Assert(fHE1.Id == fHEId1);
                    System.Diagnostics.Debug.Assert(fHE1.BHEId == hEId1);
                    System.Diagnostics.Debug.Assert(fHE1.ULId == uLId1);
                    System.Diagnostics.Debug.Assert(fHE1.UVId == uVId2);
                }
                System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId1));
            }

            uint bHEId2;
            uint fHEId2;
            uint uLId2;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId2));
                HalfEdge hE2 = HalfEdgeArray.GetObject(hEId2);
                System.Diagnostics.Debug.Assert(hE2.Id == hEId2);
                System.Diagnostics.Debug.Assert(hE2.UVId == uVId2);
                uLId2 = hE2.ULId;
                bHEId2 = hE2.BHEId;
                fHEId2 = hE2.FHEId;
                {
                    System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId2));
                    HalfEdge bHE2 = HalfEdgeArray.GetObject(bHEId2);
                    System.Diagnostics.Debug.Assert(bHE2.Id == bHEId2);
                    System.Diagnostics.Debug.Assert(bHE2.FHEId == hEId2);
                    System.Diagnostics.Debug.Assert(bHE2.ULId == uLId2);
                }
                {
                    System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(fHEId2));
                    HalfEdge fHE2 = HalfEdgeArray.GetObject(fHEId2);
                    System.Diagnostics.Debug.Assert(fHE2.Id == fHEId2);
                    System.Diagnostics.Debug.Assert(fHE2.BHEId == hEId2);
                    System.Diagnostics.Debug.Assert(fHE2.ULId == uLId2);
                    System.Diagnostics.Debug.Assert(fHE2.UVId == uVid1);
                }
                System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId2));
            }

            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId1));
                HalfEdge bHE1 = HalfEdgeArray.GetObject(bHEId1);
                if (bHE1.OHEId != fHEId2)
                {
                    return false;
                }
                {
                    System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(fHEId2));
                    HalfEdge fHE2 = HalfEdgeArray.GetObject(fHEId2);
                    System.Diagnostics.Debug.Assert(fHE2.OHEId == bHEId1);
                }
            }

            System.Diagnostics.Debug.Assert(bHEId1 != fHEId2);

            if (hEId1 != bHEId2)
            {
                {
                    HalfEdge bHE1 = HalfEdgeArray.GetObject(bHEId1);
                    bHE1.FHEId = fHEId1;
                }
                {
                    HalfEdge fHE2 = HalfEdgeArray.GetObject(fHEId2);
                    fHE2.UVId = uVId2;
                    fHE2.BHEId = bHEId2;
                }
                {
                    HalfEdge fHE1 = HalfEdgeArray.GetObject(fHEId1);
                    fHE1.BHEId = bHEId1;
                }
                {
                    HalfEdge bHE2 = HalfEdgeArray.GetObject(bHEId2);
                    bHE2.FHEId = fHEId2;
                }
            }
            else
            {
                System.Diagnostics.Debug.Assert(hEId2 == fHEId1);
                {
                    HalfEdge bHE1 = HalfEdgeArray.GetObject(bHEId1);
                    bHE1.FHEId = fHEId2;
                }
                {
                    HalfEdge fHE2 = HalfEdgeArray.GetObject(fHEId2);
                    fHE2.UVId = uVId2;
                    fHE2.BHEId = bHEId1;
                }
            }

            HalfEdgeArray.DeleteObject(hEId1);
            HalfEdgeArray.DeleteObject(hEId2);

            UseVertexArray.DeleteObject(uVid1);
            {
                UseVertex uV2 = UseVertexArray.GetObject(uVId2);
                uV2.HEId = fHEId2;
            }

            {
                UseLoop uL1 = UseLoopArray.GetObject(uLId1);
                uL1.HEId = bHEId1;
            }
            {
                UseLoop uL2 = UseLoopArray.GetObject(uLId2);
                uL2.HEId = fHEId2;
            }

            return true;
        }

        /// <summary>
        /// 頂点とエッジを追加する
        ///   辺を２つに分ける
        ///   入力ではhe2はhe1に向かいあう半辺
        ///   出力ではhe1の前にaddhe1,he1の向かいにaddhe2,addhe2の後ろにhe2
        /// </summary>
        /// <param name="addHEId1"></param>
        /// <param name="addHEId2"></param>
        /// <param name="addUVId"></param>
        /// <param name="hEId1"></param>
        /// <returns></returns>
        public bool MVE(out uint addHEId1, out uint addHEId2, out uint addUVId, uint hEId1)
        {
            HalfEdge hE1 = HalfEdgeArray.GetObject(hEId1);
            System.Diagnostics.Debug.Assert(hE1.Id == hEId1);
            uint hEId2 = hE1.OHEId;
            HalfEdge hE2 = HalfEdgeArray.GetObject(hEId2);
            System.Diagnostics.Debug.Assert(hE2.Id == hEId2);

            uint uLId1 = hE1.ULId;
            uint uLId2 = hE2.ULId;

            uint cwHEId1 = hE1.FHEId;
            uint cwHEId2 = hE2.FHEId;

            addUVId = UseVertexArray.GetFreeObjectId();
            IList<uint> hEFreeIds = HalfEdgeArray.GetFreeObjectIds(2);
            System.Diagnostics.Debug.Assert(hEFreeIds.Count == 2);
            addHEId1 = hEFreeIds[0];
            addHEId2 = hEFreeIds[1];
            {
                UseVertex addUV = new UseVertex(addUVId, addHEId1);
                uint tmpId = UseVertexArray.AddObject(addUV);
                System.Diagnostics.Debug.Assert(tmpId == addUVId);
            }
            {
                HalfEdge addHE1 = new HalfEdge(addHEId1, addUVId, cwHEId1, hEId1, hEId2, uLId1);
                uint tmpId = HalfEdgeArray.AddObject(addHE1);
                System.Diagnostics.Debug.Assert(tmpId == addHEId1);
            }
            {
                HalfEdge addHE2 = new HalfEdge(addHEId2, addUVId, cwHEId2, hEId2, hEId1, uLId2);
                uint tmpId = HalfEdgeArray.AddObject(addHE2);
                System.Diagnostics.Debug.Assert(tmpId == addHEId2);
            }
            {
                HalfEdge hE1_ = HalfEdgeArray.GetObject(hEId1);
                hE1_.FHEId = addHEId1;
                hE1_.OHEId = addHEId2;
            }
            {
                HalfEdge hE2_ = HalfEdgeArray.GetObject(hEId2);
                hE2_.FHEId = addHEId2;
                hE2_.OHEId = addHEId1;
                hE2_.EId = 0;
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(cwHEId1));
                HalfEdge cwHE1 = HalfEdgeArray.GetObject(cwHEId1);
                System.Diagnostics.Debug.Assert(cwHE1.Id == cwHEId1);
                cwHE1.BHEId = addHEId1;
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(cwHEId2));
                HalfEdge cwHE2 = HalfEdgeArray.GetObject(cwHEId2);
                System.Diagnostics.Debug.Assert(cwHE2.Id == cwHEId2);
                cwHE2.BHEId = addHEId2;
            }
            return true;
        }

        /// <summary>
        /// 頂点の削除
        /// </summary>
        /// <param name="delUVId"></param>
        /// <returns></returns>
        public bool KVEL(uint delUVId)
        {
            uint delHEId;
            uint delULId;
            {
                UseVertex uV = UseVertexArray.GetObject(delUVId);
                delHEId = uV.HEId;
                HalfEdge hE = HalfEdgeArray.GetObject(delHEId);
                delULId = hE.ULId;
                System.Diagnostics.Debug.Assert(hE.BHEId == delHEId);
                System.Diagnostics.Debug.Assert(hE.FHEId == delHEId);
                System.Diagnostics.Debug.Assert(hE.OHEId == delHEId);
            }

            uint parentULId;
            {
                UseLoop delUL = UseLoopArray.GetObject(delULId);
                parentULId = delUL.ParentULId;
                if (parentULId != delULId && parentULId != 0)
                {
                    System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(parentULId));
                    UseLoop parentUL = UseLoopArray.GetObject(parentULId);
                    System.Diagnostics.Debug.Assert(parentUL.ParentULId == parentULId);
                    System.Diagnostics.Debug.Assert(parentUL.ChildULId != 0);
                }
            }

            if (parentULId != delULId && parentULId != 0)
            {
                uint uLId = parentULId;
                while (true)
                {
                    UseLoop uL = UseLoopArray.GetObject(uLId);
                    uint childULId = uL.ChildULId;
                    System.Diagnostics.Debug.Assert(childULId != 0);
                    if (childULId == delULId)
                    {
                        UseLoop childUL = UseLoopArray.GetObject(childULId);
                        uint childChildULId = childUL.ChildULId;
                        uL.ChildULId = childChildULId;
                        break;
                    }
                    uLId = childULId;
                }
            }

            UseLoopArray.DeleteObject(delULId);
            UseVertexArray.DeleteObject(delUVId);
            HalfEdgeArray.DeleteObject(delHEId);

            return true;
        }

        /// <summary>
        /// 頂点を作る
        /// 位相ループに位相頂点を付け加える
        /// </summary>
        /// <param name="addUVId"></param>
        /// <param name="addHEId"></param>
        /// <param name="addULId"></param>
        /// <param name="uLId1"></param>
        /// <returns></returns>
        public bool MVEL(out uint addUVId, out uint addHEId, out uint addULId, uint uLId1)
        {
            addUVId = UseVertexArray.GetFreeObjectId();
            addHEId = HalfEdgeArray.GetFreeObjectId();
            addULId = UseLoopArray.GetFreeObjectId();

            {
                UseVertex addUV = new UseVertex(addUVId, addHEId);
                uint tmpId = UseVertexArray.AddObject(addUV);
                System.Diagnostics.Debug.Assert(tmpId == addUVId);
            }
            {
                HalfEdge addHE = new HalfEdge(addHEId, addUVId, addHEId, addHEId, addHEId, addULId);
                uint tmpId = HalfEdgeArray.AddObject(addHE);
                System.Diagnostics.Debug.Assert(tmpId == addHEId);
            }
            {
                UseLoop addUL = new UseLoop(addULId, addHEId, 0, uLId1);
                uint tmpId = UseLoopArray.AddObject(addUL);
                System.Diagnostics.Debug.Assert(tmpId == addULId);
            }
            if (uLId1 != 0)
            {
                uint uLId = uLId1;
                while (true)
                {
                    System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId));
                    UseLoop uL = UseLoopArray.GetObject(uLId);
                    System.Diagnostics.Debug.Assert(uL.Id == uLId);
                    if (uL.ChildULId == 0)
                    {
                        uL.ChildULId = addULId;
                        break;
                    }
                    uLId = uL.ChildULId;
                }
            }
            return true;
        }

        /// <summary>
        /// エッジを削除してループを作る
        /// </summary>
        /// <param name="addULId"></param>
        /// <param name="hEId1"></param>
        /// <returns></returns>
        public bool KEML(out uint addULId, uint hEId1)
        {
            addULId = 0;

            uint fHEId1;
            uint bHEId1;
            uint uVId1;
            uint uLId1;
            uint hEId2;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId1));
                if (!HalfEdgeArray.IsObjectId(hEId1))
                {
                    return false;
                }
                HalfEdge hE = HalfEdgeArray.GetObject(hEId1);
                System.Diagnostics.Debug.Assert(hE.Id == hEId1);
                bHEId1 = hE.BHEId;
                fHEId1 = hE.FHEId;
                uVId1 = hE.UVId;
                uLId1 = hE.ULId;
                hEId2 = hE.OHEId;
            }

            uint fHEId2;
            uint bHEId2;
            uint uVId2;
            uint uLId2;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId2));
                if (!HalfEdgeArray.IsObjectId(hEId2))
                {
                    return false;
                }
                HalfEdge hE = HalfEdgeArray.GetObject(hEId2);
                System.Diagnostics.Debug.Assert(hE.Id == hEId2);
                bHEId2 = hE.BHEId;
                fHEId2 = hE.FHEId;
                uVId2 = hE.UVId;
                uLId2 = hE.ULId;
            }
            System.Diagnostics.Debug.Assert(uLId1 == uLId2);
            {
                System.Diagnostics.Debug.Assert(fHEId1 != hEId1);
                System.Diagnostics.Debug.Assert(bHEId1 != hEId1);
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(fHEId1));
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId1));
                HalfEdge bHE = HalfEdgeArray.GetObject(bHEId1);
                System.Diagnostics.Debug.Assert(bHE.Id == bHEId1);
                System.Diagnostics.Debug.Assert(bHE.FHEId == hEId1);
                System.Diagnostics.Debug.Assert(bHE.ULId == uLId1);
            }
            {
                System.Diagnostics.Debug.Assert(fHEId2 != hEId2);
                System.Diagnostics.Debug.Assert(bHEId2 != hEId2);
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(fHEId2));
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId2));
                HalfEdge bHE = HalfEdgeArray.GetObject(bHEId2);
                System.Diagnostics.Debug.Assert(bHE.Id == bHEId2);
                System.Diagnostics.Debug.Assert(bHE.FHEId == hEId2);
                System.Diagnostics.Debug.Assert(bHE.ULId == uLId2);
            }

            addULId = UseLoopArray.GetFreeObjectId();
            uint parentULId;
            {
                System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId1));
                UseLoop uL1 = UseLoopArray.GetObject(uLId1);
                uL1.HEId = bHEId1;
                parentULId = uL1.ParentULId;
                if (parentULId != 0)
                {
                    System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(parentULId));
                    UseLoop parentUL = UseLoopArray.GetObject(parentULId);
                    System.Diagnostics.Debug.Assert(parentUL.ParentULId == parentULId);
                    uint uLId = uLId1;
                    while (true)
                    {
                        UseLoop uL = UseLoopArray.GetObject(uLId);
                        if (uL.ChildULId == 0)
                        {
                            uL.ChildULId = addULId;
                            break;
                        }
                        uLId = uL.ChildULId;
                    }
                }
            }
            {
                uint tmpId = UseLoopArray.AddObject(new UseLoop(addULId, fHEId1, 0, parentULId));
                System.Diagnostics.Debug.Assert(tmpId == addULId);
            }

            {
                System.Diagnostics.Debug.Assert(UseVertexArray.IsObjectId(uVId1));
                UseVertex uV = UseVertexArray.GetObject(uVId1);
                uV.HEId = fHEId2;
            }
            {
                System.Diagnostics.Debug.Assert(UseVertexArray.IsObjectId(uVId2));
                UseVertex uV = UseVertexArray.GetObject(uVId2);
                uV.HEId = fHEId1;
            }

            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId1));
                HalfEdge hE = HalfEdgeArray.GetObject(bHEId1);
                hE.FHEId = fHEId2;
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(fHEId2));
                HalfEdge hE = HalfEdgeArray.GetObject(fHEId2);
                hE.BHEId = bHEId1;
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(fHEId1));
                HalfEdge hE = HalfEdgeArray.GetObject(fHEId1);
                hE.BHEId = bHEId2;
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId2));
                HalfEdge hE = HalfEdgeArray.GetObject(bHEId2);
                hE.FHEId = fHEId1;
            }
            { 
                uint hEId = bHEId2;
                while (true)
                {
                    System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId));
                    HalfEdge hE = HalfEdgeArray.GetObject(hEId);
                    System.Diagnostics.Debug.Assert(hE.ULId == uLId1);
                    hE.ULId = addULId;
                    hEId = hE.FHEId;
                    if (hEId == bHEId2)
                    {
                        break;
                    }
                }
            }
            HalfEdgeArray.DeleteObject(hEId1);
            HalfEdgeArray.DeleteObject(hEId2);

            return true;
        }

        /// <summary>
        /// エッジを追加してループを削除する
        ///　 he1の起点uv1とhe2の起点uv2を結んで、２つのループをつなげる
        ///   he1は[uv1 - uv2]、he2は[uv2 - uv1]
        /// </summary>
        /// <param name="addHEId1"></param>
        /// <param name="addHEId2"></param>
        /// <param name="hEId1"></param>
        /// <param name="hEId2"></param>
        /// <returns></returns>
        public bool MEKL(out uint addHEId1, out uint addHEId2, uint hEId1, uint hEId2)
        {
            addHEId1 = 0;
            addHEId2 = 0;

            uint bHEId1;
            uint uVId1;
            uint uLId1;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId1));
                if (!HalfEdgeArray.IsObjectId(hEId1))
                {
                    return false;
                }
                HalfEdge hE = HalfEdgeArray.GetObject(hEId1);
                System.Diagnostics.Debug.Assert(hE.Id == hEId1);
                bHEId1 = hE.BHEId;
                uVId1 = hE.UVId;
                uLId1 = hE.ULId;
                System.Diagnostics.Debug.Assert(bHEId1 != hEId1);
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId1));
                HalfEdge bHE = HalfEdgeArray.GetObject(bHEId1);
                System.Diagnostics.Debug.Assert(bHE.Id == bHEId1);
                System.Diagnostics.Debug.Assert(bHE.FHEId == hEId1);
                System.Diagnostics.Debug.Assert(bHE.ULId == uLId1);
            }

            uint bHEId2;
            uint uVId2;
            uint uLId2;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId2));
                if (!HalfEdgeArray.IsObjectId(hEId2))
                {
                    return false;
                }
                HalfEdge hE = HalfEdgeArray.GetObject(hEId2);
                System.Diagnostics.Debug.Assert(hE.Id == hEId2);
                bHEId2 = hE.BHEId;
                uVId2 = hE.UVId;
                uLId2 = hE.ULId;
                System.Diagnostics.Debug.Assert(bHEId2 != hEId2);
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId2));
                HalfEdge bHE = HalfEdgeArray.GetObject(bHEId2);
                System.Diagnostics.Debug.Assert(bHE.Id == bHEId2);
                System.Diagnostics.Debug.Assert(bHE.FHEId == hEId2);
                System.Diagnostics.Debug.Assert(bHE.ULId == uLId2);
            }

            System.Diagnostics.Debug.Assert(uLId1 != uLId2);
            if (uLId1 == uLId2)
            {
                return false;
            }

            uint parentULId1;
            uint childULId1;
            {
                System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId1));
                UseLoop uL = UseLoopArray.GetObject(uLId1);
                parentULId1 = uL.ParentULId;
                childULId1 = uL.ChildULId;
            }

            uint parentULId2;
            uint childULId2;
            {
                System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId2));
                UseLoop uL = UseLoopArray.GetObject(uLId1);
                parentULId2 = uL.ParentULId;
                childULId2 = uL.ChildULId;
            }

            {
                IList<uint> freeEdgeIds = HalfEdgeArray.GetFreeObjectIds(2);
                System.Diagnostics.Debug.Assert(freeEdgeIds.Count == 2);
                addHEId1 = freeEdgeIds[0];
                addHEId2 = freeEdgeIds[1];
                System.Diagnostics.Debug.Assert(!HalfEdgeArray.IsObjectId(addHEId1));
                System.Diagnostics.Debug.Assert(!HalfEdgeArray.IsObjectId(addHEId2));
                System.Diagnostics.Debug.Assert(addHEId1 != addHEId2);
            }


            {
                HalfEdge tmpHE = new HalfEdge(addHEId1, uVId1, hEId2, bHEId1, addHEId2, uLId1);
                uint tmpId = HalfEdgeArray.AddObject(tmpHE);
                System.Diagnostics.Debug.Assert(tmpId == addHEId1);
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(addHEId1));
            }
            {
                HalfEdge tmpHE = new HalfEdge(addHEId2, uVId2, hEId1, bHEId2, addHEId1, uLId1);
                uint tmpId = HalfEdgeArray.AddObject(tmpHE);
                System.Diagnostics.Debug.Assert(tmpId == addHEId2);
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(addHEId2));
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId1));
                HalfEdge hE = HalfEdgeArray.GetObject(bHEId1);
                System.Diagnostics.Debug.Assert(hE.Id == bHEId1);
                hE.FHEId = addHEId1;
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId2));
                HalfEdge hE = HalfEdgeArray.GetObject(bHEId2);
                System.Diagnostics.Debug.Assert(hE.Id == bHEId2);
                hE.FHEId = addHEId2;
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId1));
                HalfEdge hE = HalfEdgeArray.GetObject(hEId1);
                System.Diagnostics.Debug.Assert(hE.Id == hEId1);
                hE.BHEId = addHEId2;
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId2));
                HalfEdge hE = HalfEdgeArray.GetObject(hEId2);
                System.Diagnostics.Debug.Assert(hE.Id == hEId2);
                hE.BHEId = addHEId1;
            }
            {
                // 途中で-1になるかもしれないからint型
                int hEId = (int)addHEId2;

                while (true)
                {
                    System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId((uint)hEId));
                    HalfEdge edge = HalfEdgeArray.GetObject((uint)hEId);
                    System.Diagnostics.Debug.Assert((int)edge.Id == hEId);
                    edge.ULId = uLId1;
                    hEId = (int)edge.FHEId;
                    if (hEId == (int)addHEId2)
                    {
                        break;
                    }
                }
            }

            UseLoopArray.DeleteObject(uLId2);
            System.Diagnostics.Debug.Assert(!UseLoopArray.IsObjectId(uLId2));
            if (parentULId2 != uLId2 && parentULId2 != 0)
            {
                uint uLId;
                if (parentULId1 != uLId1)
                {
                    System.Diagnostics.Debug.Assert(parentULId1 == parentULId2);
                    uLId = parentULId1;
                }
                else
                {
                    uLId = uLId1;
                }
                while (true)
                {
                    System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId));
                    UseLoop uL = UseLoopArray.GetObject(uLId);
                    if (uL.ChildULId == uLId2)
                    {
                        uL.ChildULId = childULId2;
                        break;
                    }
                    uLId = uL.ChildULId;
                    System.Diagnostics.Debug.Assert(uLId != 0);
                }
            }
            else if (parentULId1 != 0)
            {
                if (childULId2 == uLId1)
                {
                    {
                        System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId1));
                        UseLoop uL = UseLoopArray.GetObject(uLId1);
                        uL.ParentULId = uLId1;
                    }
                    uint uLId = childULId1;
                    while (true)
                    {
                        if (uLId == 0)
                        {
                            break;
                        }
                        System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId));
                        UseLoop uL = UseLoopArray.GetObject(uLId);
                        System.Diagnostics.Debug.Assert(uL.ParentULId == uLId2);
                        uL.ParentULId = uLId1;
                        uLId = uL.ChildULId;
                    }
                }
                else
                {
                    System.Diagnostics.Debug.Assert(childULId2 != 0);
                    uint uLId = childULId2;
                    while (true)
                    {
                        System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId));
                        UseLoop uL = UseLoopArray.GetObject(uLId);
                        System.Diagnostics.Debug.Assert(uL.ParentULId == uLId2);
                        uint childULId = uL.ChildULId;
                        System.Diagnostics.Debug.Assert(uL.ParentULId == uLId2);
                        if (uLId == uLId1)
                        {
                            uL.ParentULId = uLId;
                            uL.ChildULId = childULId2;
                        }
                        else
                        {
                            uL.ParentULId = uLId1;
                            if (uL.ChildULId == uLId1)
                            {
                                uL.ChildULId = childULId1;
                            }
                        }
                        uLId = childULId;
                        if (uLId == 0)
                        {
                            break;
                        }
                    }
                }
            }

            System.Diagnostics.Debug.Assert(AssertValidUse() == 0);

            return true;
        }

        /// <summary>
        /// 片方が端点である、HalfEdgeを削除する。
        /// he1の起点uv1は他の辺につながっていない端点である
        /// </summary>
        /// <param name="addULId"></param>
        /// <param name="hEId1"></param>
        /// <returns></returns>
        public bool KEMLOneFloatingVertex(out uint addULId, uint hEId1)
        {

            addULId = 0;

            uint fHEId1;
            uint hEId2;
            uint uVId1;
            uint uLId1;
            uint parentULId;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId1));
                if (!HalfEdgeArray.IsObjectId(hEId1))
                {
                    return false;
                }
                HalfEdge hE = HalfEdgeArray.GetObject(hEId1);
                System.Diagnostics.Debug.Assert(hE.Id == hEId1);
                fHEId1 = hE.FHEId;
                uVId1 = hE.UVId;
                hEId2 = hE.OHEId;
                System.Diagnostics.Debug.Assert(hE.BHEId == hEId2);
                System.Diagnostics.Debug.Assert(hE.FHEId != hEId2);
                System.Diagnostics.Debug.Assert(hEId1 != hEId2);
                uLId1 = hE.ULId;
                UseLoop uL = UseLoopArray.GetObject(uLId1);
                parentULId = uL.ParentULId;
                if (parentULId != 0)
                {
                    UseLoop parentUL = UseLoopArray.GetObject(parentULId);
                    System.Diagnostics.Debug.Assert(parentUL.ParentULId == parentULId);
                }
            }
            uint uVId2;
            uint bHEId2;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId2));
                HalfEdge hE = HalfEdgeArray.GetObject(hEId2);
                System.Diagnostics.Debug.Assert(hE.Id == hEId2);
                bHEId2 = hE.BHEId;
                uVId2 = hE.UVId;
                System.Diagnostics.Debug.Assert(bHEId2 != fHEId1);
                System.Diagnostics.Debug.Assert(hE.FHEId == hEId1);
                System.Diagnostics.Debug.Assert(hE.ULId == uLId1);
            }

            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(fHEId1));
                HalfEdge hE = HalfEdgeArray.GetObject(fHEId1);
                System.Diagnostics.Debug.Assert(hE.Id == fHEId1);
                System.Diagnostics.Debug.Assert(hE.UVId == uVId2);
                System.Diagnostics.Debug.Assert(hE.ULId == uLId1);
            }

            addULId = UseLoopArray.GetFreeObjectId();
            {
                uint tmpId = UseLoopArray.AddObject(new UseLoop(addULId, hEId1, 0, parentULId));
                System.Diagnostics.Debug.Assert(tmpId == addULId);
            }
            if (parentULId != 0)
            { 
                uint uLId = uLId1;
                while (true)
                {
                    UseLoop uL = UseLoopArray.GetObject(uLId);
                    System.Diagnostics.Debug.Assert(uL.ParentULId == parentULId);
                    if (uL.ChildULId == 0)
                    {
                        uL.ChildULId = addULId;
                        break;
                    }
                    uLId = uL.ChildULId;
                }
            }
            {
                UseLoop uL = UseLoopArray.GetObject(uLId1);
                uL.HEId = fHEId1;
            }
            {
                UseVertex uV = UseVertexArray.GetObject(uVId2);
                uV.HEId = fHEId1;
            }

            HalfEdgeArray.DeleteObject(hEId2);
            {
                HalfEdge hE = HalfEdgeArray.GetObject(hEId1);
                hE.BHEId = hEId1;
                hE.FHEId = hEId1;
                hE.OHEId = hEId1;
                hE.ULId = addULId;
                hE.UVId = uVId1;
                hE.EId = 0;
                hE.IsSameDir = true;
            }
            {
                HalfEdge hE = HalfEdgeArray.GetObject(fHEId1);
                hE.BHEId = bHEId2;
            }
            {
                HalfEdge hE = HalfEdgeArray.GetObject(bHEId2);
                hE.FHEId = fHEId1;
            }
            return true;
        }

        /// <summary>
        /// ループと浮遊点をつなげる,he1がLoop上のEdgeでhe2が浮遊点Edge
        /// he2は[uv2-uv1], he_add1は[uv1-uv2]のHalfEdgeとなる
        /// </summary>
        /// <param name="addHEId"></param>
        /// <param name="hEId1"></param>
        /// <param name="hEId2"></param>
        /// <returns></returns>
        public bool MEKLOneFloatingVertex(out uint addHEId, uint hEId1, uint hEId2)
        {
            addHEId = 0;

            uint bHEId1;
            uint uVId1;
            uint uLId1;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId1));
                if (!HalfEdgeArray.IsObjectId(hEId1))
                {
                    return false;
                }
                HalfEdge hE = HalfEdgeArray.GetObject(hEId1);
                System.Diagnostics.Debug.Assert(hE.Id == hEId1);
                bHEId1 = hE.BHEId;
                uVId1 = hE.UVId;
                uLId1 = hE.ULId;
                System.Diagnostics.Debug.Assert(bHEId1 != hEId1);
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId1));
                HalfEdge bHE = HalfEdgeArray.GetObject(bHEId1);
                System.Diagnostics.Debug.Assert(bHE.Id == bHEId1);
                System.Diagnostics.Debug.Assert(bHE.FHEId == hEId1);
                System.Diagnostics.Debug.Assert(bHE.ULId == uLId1);
            }
            uint parentULId1;
            {
                System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId1));
                UseLoop uL = UseLoopArray.GetObject(uLId1);
                parentULId1 = uL.ParentULId;
            }

            uint uVId2;
            uint uLId2;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId2));
                if (!HalfEdgeArray.IsObjectId(hEId2))
                {
                    return false;
                }
                HalfEdge hE = HalfEdgeArray.GetObject(hEId2);
                System.Diagnostics.Debug.Assert(hE.Id == hEId2);
                uVId2 = hE.UVId;
                uLId2 = hE.ULId;
                System.Diagnostics.Debug.Assert(hE.FHEId == hEId2);
                System.Diagnostics.Debug.Assert(hE.BHEId == hEId2);
            }

            System.Diagnostics.Debug.Assert(uLId1 != uLId2);
            if (uLId1 == uLId2)
            {
                return false;
            }

            {
                addHEId = HalfEdgeArray.GetFreeObjectId();
                System.Diagnostics.Debug.Assert(!HalfEdgeArray.IsObjectId(addHEId));
            }

            {
                HalfEdge tmpHE = new HalfEdge(addHEId, uVId1, hEId2, bHEId1, hEId2, uLId1);
                uint tmpId = HalfEdgeArray.AddObject(tmpHE);
                System.Diagnostics.Debug.Assert(tmpId == addHEId);
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(addHEId));
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId1));
                HalfEdge hE = HalfEdgeArray.GetObject(bHEId1);
                System.Diagnostics.Debug.Assert(hE.Id == bHEId1);
                hE.FHEId = addHEId;
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId1));
                HalfEdge hE = HalfEdgeArray.GetObject(hEId1);
                System.Diagnostics.Debug.Assert(hE.Id == hEId1);
                hE.BHEId = hEId2;
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId2));
                HalfEdge hE = HalfEdgeArray.GetObject(hEId2);
                System.Diagnostics.Debug.Assert(hE.Id == hEId2);
                hE.BHEId = addHEId;
                hE.FHEId = hEId1;
                hE.OHEId = addHEId;
                hE.ULId = uLId1;
                hE.UVId = uVId2;
            }
            {
                uint childULId0;
                {
                    System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId2));
                    UseLoop uL = UseLoopArray.GetObject(uLId2);
                    childULId0 = uL.ChildULId;
                }
                if (parentULId1 != 0)
                {
                    uint uLId = parentULId1;
                    while (true)
                    {
                        System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId));
                        UseLoop uL = UseLoopArray.GetObject(uLId);
                        uint childULId = uL.ChildULId;
                        if (childULId == uLId2)
                        {
                            uL.ChildULId = childULId0;
                            break;
                        }
                        if (childULId == 0)
                        {
                            break;
                        }
                        System.Diagnostics.Debug.Assert(childULId != childULId0);
                        uLId = childULId;
                    }
                }
            }
            UseLoopArray.DeleteObject(uLId2);
            System.Diagnostics.Debug.Assert(!UseLoopArray.IsObjectId(uLId2));
            System.Diagnostics.Debug.Assert(AssertValidUse() == 0);

            return true;
        }


        /// <summary>
        /// ループと浮遊点をつなげる,he1,he2が浮遊点Edge
        /// he1は[uv1-uv2],he2は[uv2-uv1]のHalfEdgeとなる
        /// </summary>
        /// <param name="hEId1"></param>
        /// <param name="hEId2"></param>
        /// <returns></returns>
        public bool MEKLTwoFloatingVertex(uint hEId1, uint hEId2)
        {
            uint uVId1;
            uint uLId1;
            uint parentULId1;
            uint childULId1;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId1));
                if (!HalfEdgeArray.IsObjectId(hEId1))
                {
                    return false;
                }
                HalfEdge hE = HalfEdgeArray.GetObject(hEId1);
                System.Diagnostics.Debug.Assert(hE.Id == hEId1);
                System.Diagnostics.Debug.Assert(hE.FHEId == hEId1);
                System.Diagnostics.Debug.Assert(hE.BHEId == hEId1);
                System.Diagnostics.Debug.Assert(hE.OHEId == hEId1);
                uVId1 = hE.UVId;
                uLId1 = hE.ULId;
                UseLoop uL = UseLoopArray.GetObject(uLId1);
                parentULId1 = uL.ParentULId;
                childULId1 = uL.ChildULId;
            }
            uint uVId2;
            uint uLId2;
            uint parentULId2;
            uint childULId2;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId2));
                if (!HalfEdgeArray.IsObjectId(hEId2))
                {
                    return false;
                }
                HalfEdge hE = HalfEdgeArray.GetObject(hEId2);
                System.Diagnostics.Debug.Assert(hE.Id == hEId2);
                System.Diagnostics.Debug.Assert(hE.FHEId == hEId2);
                System.Diagnostics.Debug.Assert(hE.BHEId == hEId2);
                System.Diagnostics.Debug.Assert(hE.OHEId == hEId2);
                uVId2 = hE.UVId;
                uLId2 = hE.ULId;
                UseLoop uL = UseLoopArray.GetObject(uLId2);
                parentULId2 = uL.ParentULId;
                childULId2 = uL.ChildULId;
            }
            System.Diagnostics.Debug.Assert(parentULId1 == parentULId2);

            {
                HalfEdge hE = HalfEdgeArray.GetObject(hEId1);
                hE.BHEId = hEId2;
                hE.FHEId = hEId2;
                hE.OHEId = hEId2;
                hE.ULId = uLId1;
            }
            {
                HalfEdge hE = HalfEdgeArray.GetObject(hEId2);
                hE.BHEId = hEId1;
                hE.FHEId = hEId1;
                hE.OHEId = hEId1;
                hE.ULId = uLId1;
            }

            UseLoopArray.DeleteObject(uLId2);
            if (parentULId1 != 0)
            {
                uint uLId = parentULId1;
                while (true)
                {
                    UseLoop uL = UseLoopArray.GetObject(uLId);
                    if (uL.ChildULId == uLId2)
                    {
                        uL.ChildULId = childULId2;
                        break;
                    }
                    uLId = uL.ChildULId;
                    System.Diagnostics.Debug.Assert(uLId != 0);
                }
            }

            return true;
        }

        /// <summary>
        /// 両方が端点であるEdgeを削除する
        /// </summary>
        /// <param name="addULId"></param>
        /// <param name="hEId1"></param>
        /// <returns></returns>
        public bool KEMLTwoFloatingVertex(out uint addULId, uint hEId1)
        {
            addULId = 0;

            uint hEId2;
            uint uVId1;
            uint uLIdl;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId1));
                if (!HalfEdgeArray.IsObjectId(hEId1))
                {
                    return false;
                }
                HalfEdge hE = HalfEdgeArray.GetObject(hEId1);
                System.Diagnostics.Debug.Assert(hE.Id == hEId1);
                System.Diagnostics.Debug.Assert(hE.OHEId == hE.FHEId);
                System.Diagnostics.Debug.Assert(hE.OHEId == hE.BHEId);
                uVId1 = hE.UVId;
                hEId2 = hE.OHEId;
                System.Diagnostics.Debug.Assert(hE.BHEId == hEId2);
                System.Diagnostics.Debug.Assert(hE.FHEId == hEId2);
                System.Diagnostics.Debug.Assert(hEId1 != hEId2);
                uLIdl = hE.ULId;
            }
            uint parentULId;
            uint childULId;
            {
                UseLoop uL = UseLoopArray.GetObject(uLIdl);
                parentULId = uL.ParentULId;
                childULId = uL.ChildULId;
            }
            uint uVId2;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId2));
                HalfEdge hE = HalfEdgeArray.GetObject(hEId2);
                System.Diagnostics.Debug.Assert(hE.Id == hEId2);
                System.Diagnostics.Debug.Assert(hE.BHEId == hEId1);
                System.Diagnostics.Debug.Assert(hE.FHEId == hEId1);
                System.Diagnostics.Debug.Assert(hE.OHEId == hEId1);
                uVId2 = hE.UVId;
                System.Diagnostics.Debug.Assert(hE.ULId == uLIdl);
            }

            addULId = UseLoopArray.GetFreeObjectId();
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId1));
                HalfEdge hE = HalfEdgeArray.GetObject(hEId1);
                hE.BHEId = hEId1;
                hE.FHEId = hEId1;
                hE.OHEId = hEId1;
                System.Diagnostics.Debug.Assert(hE.ULId == uLIdl);
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId2));
                HalfEdge hE = HalfEdgeArray.GetObject(hEId2);
                hE.BHEId = hEId2;
                hE.FHEId = hEId2;
                hE.OHEId = hEId2;
                hE.ULId = addULId;
            }
            {
                UseLoop uL = UseLoopArray.GetObject(uLIdl);
                uL.ChildULId = addULId;
                System.Diagnostics.Debug.Assert(uL.ParentULId == parentULId);
                // 2019-03-11 RemoveElement FIX
                //System.Diagnostics.Debug.Assert(uL.HEId == hEId1);
            }
            {
                uint tmpId = UseLoopArray.AddObject(new UseLoop(addULId, hEId2, childULId, parentULId));
                System.Diagnostics.Debug.Assert(tmpId == addULId);
            }

            return true;
        }

        /// <summary>
        /// ul2側のループが消去される
        /// </summary>
        /// <param name="hEId1"></param>
        /// <returns></returns>
        public bool KEL(uint hEId1)
        {
            uint uVId1;
            uint hEId2;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId1));
                if (!HalfEdgeArray.IsObjectId(hEId1))
                {
                    return false;
                }
                HalfEdge hE1 = HalfEdgeArray.GetObject(hEId1);
                System.Diagnostics.Debug.Assert(hE1.Id == hEId1);
                hEId2 = hE1.OHEId;
                uVId1 = hE1.UVId;
            }
            uint uVId2;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId2));
                HalfEdge hE2 = HalfEdgeArray.GetObject(hEId2);
                System.Diagnostics.Debug.Assert(hE2.Id == hEId2);
                System.Diagnostics.Debug.Assert(hE2.OHEId == hEId1);
                uVId2 = hE2.UVId;
            }

            uint bHEId1;
            uint fHEId1;
            uint uLId1;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId1));
                HalfEdge hE1 = HalfEdgeArray.GetObject(hEId1);
                System.Diagnostics.Debug.Assert(hE1.Id == hEId1);
                System.Diagnostics.Debug.Assert(hE1.UVId == uVId1);
                uLId1 = hE1.ULId;
                bHEId1 = hE1.BHEId;
                fHEId1 = hE1.FHEId;
                {
                    System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId1));
                    HalfEdge bHE1 = HalfEdgeArray.GetObject(bHEId1);
                    System.Diagnostics.Debug.Assert(bHE1.Id == bHEId1);
                    System.Diagnostics.Debug.Assert(bHE1.FHEId == hEId1);
                    System.Diagnostics.Debug.Assert(bHE1.ULId == uLId1);
                }
                {
                    System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(fHEId1));
                    HalfEdge he1f = HalfEdgeArray.GetObject(fHEId1);
                    System.Diagnostics.Debug.Assert(he1f.Id == fHEId1);
                    System.Diagnostics.Debug.Assert(he1f.BHEId == hEId1);
                    System.Diagnostics.Debug.Assert(he1f.ULId == uLId1);
                    System.Diagnostics.Debug.Assert(he1f.UVId == uVId2);
                }
                {
                    System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId1));
                }
            }

            uint bHEId2;
            uint fHEId2;
            uint uLId2;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId2));
                HalfEdge hE2 = HalfEdgeArray.GetObject(hEId2);
                System.Diagnostics.Debug.Assert(hE2.Id == hEId2);
                System.Diagnostics.Debug.Assert(hE2.UVId == uVId2);
                uLId2 = hE2.ULId;
                bHEId2 = hE2.BHEId;
                fHEId2 = hE2.FHEId;
                {
                    System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId2));
                    HalfEdge bHE2 = HalfEdgeArray.GetObject(bHEId2);
                    System.Diagnostics.Debug.Assert(bHE2.Id == bHEId2);
                    System.Diagnostics.Debug.Assert(bHE2.FHEId == hEId2);
                    System.Diagnostics.Debug.Assert(bHE2.ULId == uLId2);
                }
                {
                    System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(fHEId2));
                    HalfEdge fHE2 = HalfEdgeArray.GetObject(fHEId2);
                    System.Diagnostics.Debug.Assert(fHE2.Id == fHEId2);
                    System.Diagnostics.Debug.Assert(fHE2.BHEId == hEId2);
                    System.Diagnostics.Debug.Assert(fHE2.ULId == uLId2);
                    System.Diagnostics.Debug.Assert(fHE2.UVId == uVId1);
                }
                {
                    System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId2));
                }
            }

            System.Diagnostics.Debug.Assert(bHEId2 != fHEId1);
            System.Diagnostics.Debug.Assert(bHEId1 != fHEId2);
            System.Diagnostics.Debug.Assert(uLId1 != uLId2);

            HalfEdgeArray.DeleteObject(hEId1);
            HalfEdgeArray.DeleteObject(hEId2);
            {
                HalfEdge fHE1 = HalfEdgeArray.GetObject(fHEId1);
                fHE1.BHEId = bHEId2;
            }
            {
                HalfEdge bHE1 = HalfEdgeArray.GetObject(bHEId1);
                bHE1.FHEId = fHEId2;
            }
            {
                HalfEdge fHE2 = HalfEdgeArray.GetObject(fHEId2);
                fHE2.BHEId = bHEId1;
            }
            {
                HalfEdge bHE2 = HalfEdgeArray.GetObject(bHEId2);
                bHE2.FHEId = fHEId1;
            }
            {
                uint hEId = fHEId1;
                while (true)
                {
                    System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId));
                    HalfEdge hE = HalfEdgeArray.GetObject(hEId);
                    System.Diagnostics.Debug.Assert((int)hE.Id == hEId);
                    System.Diagnostics.Debug.Assert(hE.ULId == uLId1 || hE.ULId == uLId2);

                    hE.ULId = uLId1;

                    hEId = hE.FHEId;
                    if (hEId == fHEId1)
                    {
                        break;
                    }
                }
            }

            {
                UseVertex uV1 = UseVertexArray.GetObject(uVId1);
                uV1.HEId = fHEId2;
            }
            {
                UseVertex uV2 = UseVertexArray.GetObject(uVId2);
                uV2.HEId = fHEId1;
            }

            {
                UseLoop uL1 = UseLoopArray.GetObject(uLId1);
                uL1.HEId = fHEId1;
            }
            uint parentULId1;
            uint childULId1;
            {
                UseLoop uL1 = UseLoopArray.GetObject(uLId1);
                parentULId1 = uL1.ParentULId;
                childULId1 = uL1.ChildULId;
            }
            uint parentULId2;
            uint childULId2;
            {
                UseLoop uL2 = UseLoopArray.GetObject(uLId2);
                parentULId2 = uL2.ParentULId;
                childULId2 = uL2.ChildULId;
            }
            UseLoopArray.DeleteObject(uLId2);
            if (parentULId2 == uLId2)
            { 
                {
                    uint uLId = childULId2;
                    while (true)
                    {
                        if (uLId == 0)
                        {
                            break;
                        }
                        UseLoop uL = UseLoopArray.GetObject(uLId);
                        System.Diagnostics.Debug.Assert(uL.ParentULId == uLId2);

                        uL.ParentULId = parentULId1;

                        uLId = uL.ChildULId;
                    }
                }
                {
                    UseLoop uL1 = UseLoopArray.GetObject(uLId1);
                    uL1.HEId = fHEId1;
                    uint uLId = uLId1;
                    while (true)
                    {
                        UseLoop uL = UseLoopArray.GetObject(uLId);
                        System.Diagnostics.Debug.Assert(uL.ParentULId == parentULId1 ||
                            (uL.ParentULId == 0 && uLId == uLId1));
                        uLId = uL.ChildULId;
                        if (uLId == 0)
                        {
                            uL.ChildULId = childULId2;
                            break;
                        }
                    }
                }
            }
            else
            {
                // 2019-03-11 RemoveElement FIX
                //System.Diagnostics.Debug.Assert(parentULId1 == uLId1);
                { 
                    uint uLId = uLId1;
                    while (true)
                    {
                        System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId));
                        UseLoop uL = UseLoopArray.GetObject(uLId);
                        System.Diagnostics.Debug.Assert(uL.ParentULId == uLId1 ||
                            (uL.ParentULId == 0 && uLId == uLId1));

                        uL.ParentULId = parentULId2;

                        uLId = uL.ChildULId;
                        if (uLId == 0)
                        {
                            break;
                        }
                    }
                }
                if (parentULId2 != 0)
                {
                    uint uLId = parentULId2;
                    while (true)
                    {
                        System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId));
                        UseLoop uL = UseLoopArray.GetObject(uLId);
                        System.Diagnostics.Debug.Assert(uL.ParentULId == parentULId2 ||
                            (uL.ParentULId == 0 && uLId == parentULId2));
                        if (uL.ChildULId == uLId2)
                        {
                            uL.ChildULId = childULId2;
                        }
                        if (uL.ChildULId == 0)
                        {
                            uL.ChildULId = uLId1;
                            break;
                        }
                        uLId = uL.ChildULId;
                    }
                }
            }
            return true;
        }

        /// <summary>
        /// ループを２つに分ける
        /// </summary>
        /// <param name="addHEId1"></param>
        /// <param name="addHEId2"></param>
        /// <param name="addULId"></param>
        /// <param name="hEId1"></param>
        /// <param name="hEId2"></param>
        /// <returns></returns>
        public bool MEL(out uint addHEId1, out uint addHEId2, out uint addULId, uint hEId1, uint hEId2)
        {
            addHEId1 = 0;
            addHEId2 = 0;
            addULId = 0;

            uint uLId;
            uint bHEId1;
            uint uVId1;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId1));
                if (!HalfEdgeArray.IsObjectId(hEId1))
                {
                    return false;
                }
                HalfEdge hE1 = HalfEdgeArray.GetObject(hEId1);
                System.Diagnostics.Debug.Assert(hE1.Id == hEId1);
                uLId = hE1.ULId;
                uVId1 = hE1.UVId;
                bHEId1 = hE1.BHEId;
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId1));
                HalfEdge bHE1 = HalfEdgeArray.GetObject(bHEId1);
                System.Diagnostics.Debug.Assert(bHE1.Id == bHEId1);
                System.Diagnostics.Debug.Assert(bHE1.FHEId == hEId1);
                System.Diagnostics.Debug.Assert(bHE1.ULId == uLId);
            }

            uint bHEId2;
            uint uVId2;
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId2));
                if (!HalfEdgeArray.IsObjectId(hEId2))
                {
                    return false;
                }
                HalfEdge hE2 = HalfEdgeArray.GetObject(hEId2);
                System.Diagnostics.Debug.Assert(hE2.Id == hEId2);
                System.Diagnostics.Debug.Assert(hE2.ULId == uLId);
                uVId2 = hE2.UVId;
                bHEId2 = hE2.BHEId;
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId2));
                HalfEdge bHE2 = HalfEdgeArray.GetObject(bHEId2);
                System.Diagnostics.Debug.Assert(bHE2.Id == bHEId2);
                System.Diagnostics.Debug.Assert(bHE2.FHEId == hEId2);
                System.Diagnostics.Debug.Assert(bHE2.ULId == uLId);
            }

            addULId = UseLoopArray.GetFreeObjectId();

            {
                IList<uint> freeEdgeIds = HalfEdgeArray.GetFreeObjectIds(2);
                System.Diagnostics.Debug.Assert(freeEdgeIds.Count == 2);
                addHEId1 = freeEdgeIds[0];
                addHEId2 = freeEdgeIds[1];
                System.Diagnostics.Debug.Assert(!HalfEdgeArray.IsObjectId(addHEId1));
                System.Diagnostics.Debug.Assert(!HalfEdgeArray.IsObjectId(addHEId2));
                System.Diagnostics.Debug.Assert(addHEId1 != addHEId2);
            }

            {
                uint tmpId = UseLoopArray.AddObject(new UseLoop(addULId, addHEId2, 0, 0));
                System.Diagnostics.Debug.Assert(addULId == tmpId);
                System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(addULId));
            }
            {
                System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId));
                UseLoop loop = UseLoopArray.GetObject(uLId);
                System.Diagnostics.Debug.Assert(loop.Id == uLId);
                loop.HEId = addHEId1;
            }

            {
                HalfEdge tmpHE = new HalfEdge(addHEId1, uVId1, hEId2, bHEId1, addHEId2, uLId);
                uint tmpId = HalfEdgeArray.AddObject(tmpHE);
                System.Diagnostics.Debug.Assert(tmpId == addHEId1);
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(addHEId1));
            }
            {
                HalfEdge tmp_he = new HalfEdge(addHEId2, uVId2, hEId1, bHEId2, addHEId1, addULId);
                uint tmpId = HalfEdgeArray.AddObject(tmp_he);
                System.Diagnostics.Debug.Assert(tmpId == addHEId2);
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(addHEId2));
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId1));
                HalfEdge bHE1 = HalfEdgeArray.GetObject(bHEId1);
                System.Diagnostics.Debug.Assert(bHE1.Id == bHEId1);
                bHE1.FHEId = addHEId1;
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId2));
                HalfEdge hE2 = HalfEdgeArray.GetObject(hEId2);
                System.Diagnostics.Debug.Assert(hE2.Id == hEId2);
                hE2.BHEId = addHEId1;
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(bHEId2));
                HalfEdge bHE2 = HalfEdgeArray.GetObject(bHEId2);
                System.Diagnostics.Debug.Assert(bHE2.Id == bHEId2);
                bHE2.FHEId = addHEId2;
            }
            {
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId1));
                HalfEdge hE1 = HalfEdgeArray.GetObject(hEId1);
                System.Diagnostics.Debug.Assert(hE1.Id == hEId1);
                hE1.BHEId = addHEId2;
            }
            {
                uint hEId = addHEId2;
                while (true)
                {
                    System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId));
                    HalfEdge edge = HalfEdgeArray.GetObject(hEId);
                    System.Diagnostics.Debug.Assert(edge.Id == hEId);
                    edge.ULId = addULId;
                    hEId = edge.FHEId;
                    if (hEId == addHEId2)
                    {
                        break;
                    }
                }
            }

            return true;
        }

        public bool MoveUseLoop(uint uLId1, uint uLId2)
        {
            System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId1));
            System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId2));
            uint parentULId1;
            uint childULId1;
            {
                UseLoop tmpUL1 = UseLoopArray.GetObject(uLId1);
                parentULId1 = tmpUL1.ParentULId;
                System.Diagnostics.Debug.Assert(parentULId1 != uLId1);
                childULId1 = tmpUL1.ChildULId;
            }
            uint childULId2;
            uint parentULId2;
            {
                UseLoop uL2 = UseLoopArray.GetObject(uLId2);
                childULId2 = uL2.ChildULId;
                parentULId2 = uL2.ParentULId;
            }

            UseLoop uL1 = UseLoopArray.GetObject(uLId1);
            uL1.ParentULId = parentULId2;
            uL1.ChildULId = 0;

            {
                uint uLId = uLId2;
                while (true)
                {
                    UseLoop uL = UseLoopArray.GetObject(uLId);
                    System.Diagnostics.Debug.Assert(uL.ParentULId == parentULId2 || (uLId == parentULId2 && uL.ParentULId == 0));
                    if (uL.ChildULId == 0)
                    {
                        uL.ChildULId = uLId1;
                        break;
                    }
                    uLId = uL.ChildULId;
                }
            }
            {
                uint uLId = parentULId1;
                while (true)
                {
                    UseLoop uL = UseLoopArray.GetObject(uLId);
                    System.Diagnostics.Debug.Assert(uL.ParentULId == parentULId1 || (uLId == parentULId1 && uL.ParentULId == 0));
                    if (uL.ChildULId == uLId1)
                    {
                        uL.ChildULId = childULId1;
                        break;
                    }
                    System.Diagnostics.Debug.Assert(uL.ChildULId != 0);
                    uLId = uL.ChildULId;
                }
            }

            return true;
        }

        public bool SwapUseLoop(uint uLId1, uint uLId2)
        {
            System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId1));
            System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId2));
            uint parentULId1;
            {
                UseLoop uL1 = UseLoopArray.GetObject(uLId1);
                parentULId1 = uL1.ParentULId;
            }
            uint parentULId2;
            {
                UseLoop uL2 = UseLoopArray.GetObject(uLId2);
                parentULId2 = uL2.ParentULId;
            }
            if (parentULId1 == uLId1 && parentULId2 == uLId1)
            {
                return SwapChildUseLoopParentSameLoop(uLId2, uLId1);
            }
            else if (parentULId2 == uLId2 && parentULId1 == uLId2)
            {
                return SwapChildUseLoopParentSameLoop(uLId1, uLId2);
            }
            else if ((parentULId1 == uLId1 || parentULId1 == 0) && parentULId2 != uLId1)
            {
                return SwapChildUseLoopParentDifferentLoop(uLId2, uLId1);
            }
            else if ((parentULId2 == uLId2 || parentULId2 == 0) && parentULId1 != uLId2)
            {
                return SwapChildUseLoopParentDifferentLoop(uLId1, uLId2);
            }
            System.Diagnostics.Debug.Assert(false);
            return true;
        }

        private bool SwapChildUseLoopParentDifferentLoop(uint uLId1, uint uLId2)
        {
            System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId1));
            System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId2));
            uint parentULId1;
            uint childULId1;
            {
                UseLoop uL1 = UseLoopArray.GetObject(uLId1);
                parentULId1 = uL1.ParentULId;
                System.Diagnostics.Debug.Assert(parentULId1 != uLId1);
                System.Diagnostics.Debug.Assert(uLId2 != parentULId1);
                childULId1 = uL1.ChildULId;
            }
            uint childULId2;
            {
                System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId2));
                UseLoop uL2 = UseLoopArray.GetObject(uLId2);
                System.Diagnostics.Debug.Assert(uL2.ParentULId == uLId2 || uL2.ParentULId == 0);
                childULId2 = uL2.ChildULId;
            }
            uint uLId = parentULId1;
            while (true)
            {
                System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId));
                UseLoop uL = UseLoopArray.GetObject(uLId);
                System.Diagnostics.Debug.Assert(uL.Id == uLId);
                uint childULId = uL.ChildULId;

                if (uLId == uLId1)
                {
                    uL.ChildULId = childULId2;
                    uL.ParentULId = uLId;
                }
                else
                {
                    if (uL.ChildULId == uLId1)
                    {
                        uL.ChildULId = uLId2;
                    }
                    System.Diagnostics.Debug.Assert(uL.ParentULId == parentULId1 ||
                        (uLId == parentULId1 && uL.ParentULId == 0));
                }

                if (childULId == 0)
                {
                    break;
                }
                uLId = childULId;
            }

            uLId = uLId2;
            while (true)
            {
                System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId));
                UseLoop uL = UseLoopArray.GetObject(uLId);
                System.Diagnostics.Debug.Assert(uL.Id == uLId);
                uint childULId = uL.ChildULId;

                if (uLId == uLId2)
                {
                    uL.ChildULId = childULId1;
                    uL.ParentULId = parentULId1;
                }
                else
                {
                    uL.ParentULId = uLId1;
                }

                if (childULId == 0)
                {
                    break;
                }
                uLId = childULId;
            }

            return true;
        }

        private bool SwapChildUseLoopParentSameLoop(uint uLId1, uint uLId2)
        {
            System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId1));
            System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId2));
            uint chidlULId1;
            {
                UseLoop uL1 = UseLoopArray.GetObject(uLId1);
                System.Diagnostics.Debug.Assert(uLId2 == uL1.ParentULId);
                chidlULId1 = uL1.ChildULId;
            }
            uint childULId2;
            {
                System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId2));
                UseLoop uL2 = UseLoopArray.GetObject(uLId2);
                System.Diagnostics.Debug.Assert(uL2.ParentULId == uLId2);
                childULId2 = uL2.ChildULId;
                System.Diagnostics.Debug.Assert(childULId2 != 0);
            }

            uint uLId = uLId2;
            while (true)
            {
                System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId));
                UseLoop uL = UseLoopArray.GetObject(uLId);
                System.Diagnostics.Debug.Assert(uL.Id == uLId);
                uint childULId = uL.ChildULId;

                if (uLId == uLId2)
                {
                    System.Diagnostics.Debug.Assert(chidlULId1 != uLId2);
                    uL.ChildULId = chidlULId1;
                    uL.ParentULId = uLId1;
                }
                else if (uLId == uLId1)
                {
                    if (childULId2 != uLId1)
                    {
                        uL.ChildULId = childULId2;
                    }
                    else
                    {
                        uL.ChildULId = uLId2;
                    }
                    uL.ParentULId = uLId;
                }
                else
                {
                    System.Diagnostics.Debug.Assert(uL.ChildULId != uLId2);
                    if (uL.ChildULId == uLId1)
                    {
                        uL.ChildULId = uLId2;
                    }
                    uL.ParentULId = uLId1;
                }

                if (childULId == 0)
                {
                    break;
                }
                uLId = childULId;
            }

            return true;
        }

        /////////////////////////////////////////////////////////////////////////////
        // Radial Edge
        private void SplitEdgeToRadialEdges(uint eId)
        {
            IList<uint> hEIds = FindHalfEdgeByEdge(eId);
            if (hEIds.Count != 2)
            {
                // すでにRadial Edgeになっている?
                System.Diagnostics.Debug.Assert(false);
                return;
            }

            uint hEId1 = hEIds[0];
            uint hEId2 = hEIds[1];
            HalfEdge hE1 = HalfEdgeArray.GetObject(hEId1);
            HalfEdge hE2 = HalfEdgeArray.GetObject(hEId2);
            System.Diagnostics.Debug.Assert(hEId1 == hE2.OHEId);
            System.Diagnostics.Debug.Assert(hEId2 == hE1.OHEId);
            System.Diagnostics.Debug.Assert(hE1.RadialHEId == 0);
            System.Diagnostics.Debug.Assert(hE2.RadialHEId == 0);
            System.Diagnostics.Debug.Assert(hE1.EId == eId);
            System.Diagnostics.Debug.Assert(hE2.EId == eId);

            IList<uint> hEFreeIds = HalfEdgeArray.GetFreeObjectIds(2);
            System.Diagnostics.Debug.Assert(hEFreeIds.Count == 2);
            uint addHEId1 = hEFreeIds[0];
            uint addHEId2 = hEFreeIds[1];
            // hE2の隣接関係をaddHE2に移動する
            {
                HalfEdge addHE1 = new HalfEdge(addHEId1, hE1.UVId, addHEId1, addHEId1, addHEId2, 0);
                addHE1.RadialHEId = hEId2;
                addHE1.EId = eId;
                addHE1.IsSameDir = hE1.IsSameDir;
                uint tmpId = HalfEdgeArray.AddObject(addHE1);
                System.Diagnostics.Debug.Assert(tmpId == addHEId1);
            }
            {
                HalfEdge addHE2 = new HalfEdge(addHEId2, hE2.UVId, hE2.FHEId, hE2.BHEId, addHEId1, hE2.ULId);
                addHE2.RadialHEId = hEId1;
                addHE2.EId = eId;
                addHE2.IsSameDir = hE2.IsSameDir;
                uint tmpId = HalfEdgeArray.AddObject(addHE2);
                System.Diagnostics.Debug.Assert(tmpId == addHEId2);
            }
            {
                hE1.RadialHEId = addHEId2;
                hE2.RadialHEId = addHEId1;
            }
            // hE2の隣接関係をaddHE2に置換する
            {
                UseLoop uL = UseLoopArray.GetObject(hE2.ULId);
                if (uL.HEId == hEId2)
                {
                    uL.HEId = addHEId2;
                }
                else
                {
                    // そのまま
                }
            }
            {
                HalfEdge fHE = HalfEdgeArray.GetObject(hE2.FHEId);
                fHE.BHEId = addHEId2;
            }
            {
                HalfEdge bHE = HalfEdgeArray.GetObject(hE2.BHEId);
                bHE.FHEId = addHEId2;
            }
            // hE2の前後隣接関係をクリアする
            {
                hE2.ULId = 0;
                hE2.FHEId = hE2.Id;
                hE2.BHEId = hE2.Id;
            }
        }

        /*
        private void UnSplitRadialEdgesToEdge(uint eId)
        {
            // TODO:
        }
        */

        private IList<uint> AddRadialEdgeToEdge(uint eId)
        {
            IList<uint> hEIds0 = FindHalfEdgeByEdge(eId);
            if (hEIds0.Count == 2)
            {
                SplitEdgeToRadialEdges(eId);
            }

            IList<uint> hEIds = FindHalfEdgeByEdge(eId);
            System.Diagnostics.Debug.Assert(hEIds.Count >= 4);

            uint hEId1 = hEIds[0];
            uint hEId2 = hEIds[1];
            uint hEId3 = hEIds[hEIds.Count - 2];
            uint hEId4 = hEIds[hEIds.Count - 1];
            HalfEdge hE1 = HalfEdgeArray.GetObject(hEId1);
            HalfEdge hE2 = HalfEdgeArray.GetObject(hEId2);
            HalfEdge hE3 = HalfEdgeArray.GetObject(hEId3);
            HalfEdge hE4 = HalfEdgeArray.GetObject(hEId4);
            System.Diagnostics.Debug.Assert(hEId3 == hE4.OHEId);
            System.Diagnostics.Debug.Assert(hEId4 == hE3.OHEId);
            System.Diagnostics.Debug.Assert(hE1.RadialHEId == hEId4);
            System.Diagnostics.Debug.Assert(hE4.RadialHEId == hEId1);
            System.Diagnostics.Debug.Assert(hE1.EId == eId);
            System.Diagnostics.Debug.Assert(hE2.EId == eId);
            System.Diagnostics.Debug.Assert(hE3.EId == eId);
            System.Diagnostics.Debug.Assert(hE4.EId == eId);

            IList<uint> hEFreeIds = HalfEdgeArray.GetFreeObjectIds(2);
            System.Diagnostics.Debug.Assert(hEFreeIds.Count == 2);
            uint addHEId1 = hEFreeIds[0];
            uint addHEId2 = hEFreeIds[1];
            {
                HalfEdge addHE1 = new HalfEdge(addHEId1, hE1.UVId, addHEId1, addHEId1, addHEId2, 0);
                addHE1.RadialHEId = hEId4;
                addHE1.EId = eId;
                addHE1.IsSameDir = hE1.IsSameDir;
                uint tmpId = HalfEdgeArray.AddObject(addHE1);
                System.Diagnostics.Debug.Assert(tmpId == addHEId1);
            }
            {
                HalfEdge addHE2 = new HalfEdge(addHEId2, hE2.UVId, addHEId2, addHEId2, addHEId1, 0);
                addHE2.RadialHEId = hEId1;
                addHE2.EId = eId;
                addHE2.IsSameDir = hE2.IsSameDir;
                uint tmpId = HalfEdgeArray.AddObject(addHE2);
                System.Diagnostics.Debug.Assert(tmpId == addHEId2);
            }
            {
                hE1.RadialHEId = addHEId2;
                hE4.RadialHEId = addHEId1;
            }

            IList<uint> addRadialHEIds = new List<uint>();
            addRadialHEIds.Add(addHEId1);
            addRadialHEIds.Add(addHEId2);
            return addRadialHEIds;
        }

        public uint MakeRadialUseLoop(IList<uint> eIds, IList<bool> isSameDirs)
        {
            int edgeCnt = eIds.Count;

            IList<IList<uint>> radialHEIdss = new List<IList<uint>>();
            for (int i = 0; i < edgeCnt; i++)
            {
                uint eId = eIds[i];
                IList<uint> radialHEIds = AddRadialEdgeToEdge(eId);
                System.Diagnostics.Debug.Assert(radialHEIds.Count == 2);
                radialHEIdss.Add(radialHEIds);
            }

            IList<uint> hitRadialHEIds = new List<uint>();
            for (int i = 0; i < edgeCnt; i++)
            {
                bool inIsSameDir = isSameDirs[i];
                IList<uint> hEIds = radialHEIdss[i];
                uint hEId1 = hEIds[0];
                uint hEId2 = hEIds[1];

                HalfEdge hE1 = HalfEdgeArray.GetObject(hEId1);
                HalfEdge hE2 = HalfEdgeArray.GetObject(hEId2);
                uint uVId1 = hE1.UVId;
                uint uVId2 = hE2.UVId;
                bool isSameDir1 = hE1.IsSameDir;
                bool isSameDir2 = hE2.IsSameDir;

                uint hitHEId = 0; 
                if (inIsSameDir == isSameDir1)
                {
                    hitHEId = hEId1;
                }
                else if (inIsSameDir == isSameDir2)
                {
                    hitHEId = hEId2;
                }
                else
                {
                    // fail
                    System.Diagnostics.Debug.Assert(false);
                }
                hitRadialHEIds.Add(hitHEId);
            }

            for (int i = 0; i < edgeCnt; i++)
            {
                uint curHEId = hitRadialHEIds[i];
                uint prevHEId = hitRadialHEIds[(i - 1) >= 0 ? (i - 1) : (edgeCnt - 1)];
                uint nextHEId = hitRadialHEIds[(i + 1) % edgeCnt];

                HalfEdge curHE = HalfEdgeArray.GetObject(curHEId);
                curHE.BHEId = prevHEId;
                curHE.FHEId = nextHEId;
            }

            uint addULId = UseLoopArray.GetFreeObjectId();
            uint parentULId = 0; //!!!!!!
            uint hEId0 = hitRadialHEIds[0];
            {
                UseLoop addUL = new UseLoop(addULId, hEId0, 0, parentULId);
                uint tmpId = UseLoopArray.AddObject(addUL);
                System.Diagnostics.Debug.Assert(tmpId == addULId);
            }

            for (int i = 0; i < edgeCnt; i++)
            {
                uint curHEId = hitRadialHEIds[i];
                HalfEdge curHE = HalfEdgeArray.GetObject(curHEId);
                curHE.ULId = addULId;
            }

            return addULId;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////

        public int AssertValidUse()
        {
            IList<uint> uVIds = UseVertexArray.GetObjectIds();
            for (int i = 0; i < uVIds.Count; i++)
            {
                uint uVId = uVIds[i];
                System.Diagnostics.Debug.Assert(UseVertexArray.IsObjectId(uVId));
                UseVertex uV = UseVertexArray.GetObject(uVId);
                System.Diagnostics.Debug.Assert(uV.Id == uVId);
                uint hEId = uV.HEId;
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId));
                HalfEdge hEdge = HalfEdgeArray.GetObject(hEId);
                System.Diagnostics.Debug.Assert(hEdge.Id == hEId);
                System.Diagnostics.Debug.Assert(hEdge.UVId == uVId);
            }

            IList<uint> hEdgeIds = HalfEdgeArray.GetObjectIds();
            for (int i = 0; i < hEdgeIds.Count; i++)
            {
                uint hEId = hEdgeIds[i];
                System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId));
                HalfEdge hEdge = HalfEdgeArray.GetObject(hEId);
                System.Diagnostics.Debug.Assert(hEdge.Id == hEId);

                {
                    uint uVId1 = hEdge.UVId;
                    System.Diagnostics.Debug.Assert(UseVertexArray.IsObjectId(uVId1));
                    UseVertex uV = UseVertexArray.GetObject(uVId1);
                    System.Diagnostics.Debug.Assert(uV.Id == uVId1);
                }

                uint uVId2 = 0;
                if (hEdge.RadialHEId != 0)
                {
                    // FIXME:
                    // uVId2をOHEIdからとってくることは出来る
                }
                else
                {
                    // original
                    uint fHEId = hEdge.FHEId;
                    System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(fHEId));
                    HalfEdge cwEdge = HalfEdgeArray.GetObject(fHEId);
                    System.Diagnostics.Debug.Assert(cwEdge.Id == fHEId);
                    System.Diagnostics.Debug.Assert(cwEdge.BHEId == hEId);
                    System.Diagnostics.Debug.Assert(cwEdge.ULId == hEdge.ULId);
                    uVId2 = cwEdge.UVId;
                    System.Diagnostics.Debug.Assert(UseVertexArray.IsObjectId(uVId2));
                    UseVertex uV = UseVertexArray.GetObject(uVId2);
                    System.Diagnostics.Debug.Assert(uV.Id == uVId2);
                }

                {
                    uint ccwHEId = hEdge.BHEId;
                    System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(ccwHEId));
                    HalfEdge ccwEdge = HalfEdgeArray.GetObject(ccwHEId);
                    System.Diagnostics.Debug.Assert(ccwEdge.Id == ccwHEId);
                    System.Diagnostics.Debug.Assert(ccwEdge.FHEId == hEId);
                    System.Diagnostics.Debug.Assert(ccwEdge.ULId == hEdge.ULId);
                }

                if (hEdge.RadialHEId != 0)
                {
                    // FIXME:
                }
                else
                {
                    // original
                    uint oHEId = hEdge.OHEId;
                    System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(oHEId));
                    HalfEdge oEdge = HalfEdgeArray.GetObject(oHEId);
                    System.Diagnostics.Debug.Assert(oEdge.Id == oHEId);
                    System.Diagnostics.Debug.Assert(oEdge.OHEId == hEId);
                    System.Diagnostics.Debug.Assert(oEdge.UVId == uVId2);
                }
            }

            IList<uint> uLIds = UseLoopArray.GetObjectIds();
            for (int i = 0; i < uLIds.Count; i++)
            {
                uint uLId = uLIds[i];
                System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId));
                UseLoop uL = UseLoopArray.GetObject(uLId);
                System.Diagnostics.Debug.Assert(uL.Id == uLId);
                {
                    IList<uint> passedHEdge = new List<uint>();
                    uint hEId = uL.HEId;
                    while (true)
                    {
                        System.Diagnostics.Debug.Assert(passedHEdge.IndexOf(hEId) == -1);
                        passedHEdge.Add(hEId);
                        System.Diagnostics.Debug.Assert(HalfEdgeArray.IsObjectId(hEId));
                        HalfEdge hE = HalfEdgeArray.GetObject(hEId);
                        System.Diagnostics.Debug.Assert(hE.Id == hEId);
                        System.Diagnostics.Debug.Assert(hE.ULId == uLId);
                        uint nextHEId = hE.FHEId;
                        if (nextHEId == uL.HEId)
                        {
                            break;
                        }
                        hEId = nextHEId;
                    }
                }
                if (uL.ParentULId != uLId && uL.ParentULId != 0)
                {
                    uint parentULId = uL.ParentULId;
                    uint uLId2 = parentULId;
                    bool flag = false;
                    while (true)
                    {
                        System.Diagnostics.Debug.Assert(UseLoopArray.IsObjectId(uLId2));
                        UseLoop uL2 = UseLoopArray.GetObject(uLId2);
                        System.Diagnostics.Debug.Assert(uL2.ParentULId == parentULId);
                        if (uLId2 == uLId)
                        {
                            flag = true;
                        }
                        uLId2 = uL2.ChildULId;
                        if (uLId2 == 0)
                        {
                            break;
                        }
                    }
                    System.Diagnostics.Debug.Assert(flag == true);
                }
            }
            return 0;
        }
    }
}
