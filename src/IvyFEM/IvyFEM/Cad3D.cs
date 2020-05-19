using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Cad3D
    {
        protected ObjectArray<Loop3D> LoopArray = new ObjectArray<Loop3D>();
        protected ObjectArray<Edge3D> EdgeArray = new ObjectArray<Edge3D>();
        protected ObjectArray<Vertex3D> VertexArray = new ObjectArray<Vertex3D>();
        protected BRep2D BRep = new BRep2D();

        public Cad3D()
        {

        }

        public Cad3D(Cad3D src)
        {
            Copy(src);
        }

        public void Copy(Cad3D src)
        {
            Clear();

            // Vertex
            {
                IList<uint> srcIds = src.VertexArray.GetObjectIds();
                foreach (uint srcId in srcIds)
                {
                    Vertex3D srcV = src.VertexArray.GetObject(srcId);
                    uint id = srcId;
                    Vertex3D v = new Vertex3D(srcV);
                    uint tmpId = VertexArray.AddObject(id, v);
                    System.Diagnostics.Debug.Assert(tmpId == id);
                }
            }
            // Edge
            {
                IList<uint> srcIds = src.EdgeArray.GetObjectIds();
                foreach (uint srcId in srcIds)
                {
                    Edge3D srcE = src.EdgeArray.GetObject(srcId);
                    uint id = srcId;
                    Edge3D e = new Edge3D(srcE);
                    uint tmpId = EdgeArray.AddObject(id, e);
                    System.Diagnostics.Debug.Assert(tmpId == id);
                }
            }
            // Loop
            {
                IList<uint> srcIds = src.LoopArray.GetObjectIds();
                foreach (uint srcId in srcIds)
                {
                    Loop3D srcL = src.LoopArray.GetObject(srcId);
                    uint id = srcId;
                    Loop3D e = new Loop3D(srcL);
                    uint tmpId = LoopArray.AddObject(id, e);
                    System.Diagnostics.Debug.Assert(tmpId == id);
                }
            }

            BRep.Copy(src.BRep);

            AssertValid();
        }

        public void Clear()
        {
            LoopArray.Clear();
            EdgeArray.Clear();
            VertexArray.Clear();
            BRep.Clear();
        }

        public AddVertexRes AddVertex(CadElementType type, uint id, OpenTK.Vector3d vec)
        {
            AddVertexRes res = new AddVertexRes();

            if (type == CadElementType.NotSet || id == 0)
            {
                uint addVId = BRep.AddVertexToLoop(0);
                uint tmpId = VertexArray.AddObject(addVId, new Vertex3D(vec));
                System.Diagnostics.Debug.Assert(tmpId == addVId);
                res.AddVId = addVId;
                return res;
            }
            else if (type == CadElementType.Loop)
            {
                uint lId = id;
                System.Diagnostics.Debug.Assert(LoopArray.IsObjectId(lId));
                if (!LoopArray.IsObjectId(lId))
                {
                    return res;
                }

                throw new NotImplementedException();
            }
            else if (type == CadElementType.Edge)
            {
                uint eId = id;
                System.Diagnostics.Debug.Assert(EdgeArray.IsObjectId(eId));
                if (!EdgeArray.IsObjectId(eId))
                {
                    return res;
                }

                throw new NotImplementedException();
            }

            return res;
        }

        public ConnectVertexRes ConnectVertexLine(uint vId1, uint vId2)
        {
            ConnectVertexRes res = new ConnectVertexRes();
            res.VId1 = vId1;
            res.VId2 = vId2;
            uint maxId = 0;
            {
                IList<uint> lIds = BRep.GetElementIds(CadElementType.Loop);
                for (int i = 0; i < lIds.Count; i++)
                {
                    maxId = (lIds[i] > maxId) ? lIds[i] : maxId;
                }
            }
            IList<int> flgs = new List<int>();
            for (int i = 0; i < (maxId + 1); i++)
            {
                flgs.Add(0);
            }
            VertexEdgeItr vItr0 = BRep.GetVertexEdgeItr(vId1);
            for (; !vItr0.IsEnd(); vItr0.Next())
            {
                uint lId0 = vItr0.GetLoopId();
                System.Diagnostics.Debug.Assert(BRep.IsElementId(CadElementType.Loop, lId0));
                System.Diagnostics.Debug.Assert(flgs[(int)lId0] == 0);
                flgs[(int)lId0] = 1;
            }
            uint lId = 0;
            VertexEdgeItr vItr1 = BRep.GetVertexEdgeItr(vId2);
            for (; !vItr1.IsEnd(); vItr1.Next())
            {
                uint lId0 = vItr1.GetLoopId();
                System.Diagnostics.Debug.Assert(BRep.IsElementId(CadElementType.Loop, lId0));
                if (flgs[(int)lId0] == 1)
                {
                    lId = lId0;
                    break;
                }
            }
            res.LId = lId;
            if (lId == 0)
            {
                return res;
            }
            for (vItr0.Begin(); !vItr0.IsEnd(); vItr0.Next())
            {
                uint lId0 = vItr0.GetLoopId();
                System.Diagnostics.Debug.Assert(BRep.IsElementId(CadElementType.Loop, lId0));
                if (lId0 == lId)
                {
                    break;
                }
            }
            System.Diagnostics.Debug.Assert(vItr0.GetLoopId() == lId);
            System.Diagnostics.Debug.Assert(vItr1.GetLoopId() == lId);
            res = BRep.ConnectVertex(vItr0, vItr1, true);
            {
                EdgeArray.AddObject(res.AddEId, new Edge3D());
            }
            if (res.AddLId != lId && res.AddLId != 0)
            {
                if (LoopArray.IsObjectId(res.LId))
                {
                    Loop3D addLoop = LoopArray.GetObject(lId);
                    addLoop = new Loop3D(addLoop); // copy
                    LoopArray.AddObject(res.AddLId, addLoop);
                }
                else
                {
                    LoopArray.AddObject(res.AddLId, new Loop3D());
                }
            }
            System.Diagnostics.Debug.Assert(AssertValid() == 0);
            return res;
        }
        
        public AddPolygonRes AddPolygon(IList<OpenTK.Vector3d> points, uint lId = 0)
        {
            AddPolygonRes res = new AddPolygonRes();

            if (points.Count < 3)
            {
                return res;
            }
            int ptCnt = points.Count;
            IList<uint> vIds = new List<uint>();
            for (int i = 0; i < ptCnt; i++)
            {
                vIds.Add(0);
            }
            for (int i = 0; i < ptCnt; i++)
            {
                uint uvId = BRep.AddVertexToLoop(lId);
                vIds[i] = VertexArray.AddObject(uvId, new Vertex3D(points[i]));
                res.VIds.Add(vIds[i]);
            }
            IList<uint> eIds = new List<uint>();
            for (int i = 0; i < ptCnt; i++)
            {
                eIds.Add(0);
            }
            for (int i = 0; i < ptCnt; i++)
            {
                int j = (i == ptCnt - 1) ? 0 : i + 1;
                var vItr0 = new VertexEdgeItr(BRep, vIds[i]);
                System.Diagnostics.Debug.Assert(vItr0.GetLoopId() == lId);
                var vItr1 = new VertexEdgeItr(BRep, vIds[j]);
                System.Diagnostics.Debug.Assert(vItr1.GetLoopId() == lId);
                eIds[i] = BRep.ConnectVertex(vItr0, vItr1, true).AddEId;
                EdgeArray.AddObject(eIds[i], new Edge3D());
                res.EIds.Add(eIds[i]);
            }
            uint addLId = BRep.GetEdgeLoopId(eIds[ptCnt - 1], true);
            res.AddLId = addLId;
            if (LoopArray.IsObjectId(res.AddLId))
            {
                Loop3D addLoop = LoopArray.GetObject(res.AddLId);
                addLoop = new Loop3D(addLoop); // copy
                LoopArray.AddObject(addLId, addLoop);
            }
            else
            {
                LoopArray.AddObject(res.AddLId, new Loop3D());
            }
            System.Diagnostics.Debug.Assert(AssertValid() == 0);

            return res;
        }

        public bool IsElementId(CadElementType type, uint id)
        {
            if (type == CadElementType.NotSet)
            {
                return false;
            }
            else if (type == CadElementType.Vertex)
            {
                return VertexArray.IsObjectId(id);
            }
            else if (type == CadElementType.Edge)
            {
                return EdgeArray.IsObjectId(id);
            }
            else if (type == CadElementType.Loop)
            {
                return LoopArray.IsObjectId(id);
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return false;
        }

        public IList<uint> GetElementIds(CadElementType type)
        {
            if (type == CadElementType.Vertex)
            {
                return VertexArray.GetObjectIds();
            }
            else if (type == CadElementType.Edge)
            {
                return EdgeArray.GetObjectIds();
            }
            else if (type == CadElementType.Loop)
            {
                return LoopArray.GetObjectIds();
            }

            System.Diagnostics.Debug.Assert(false);
            IList<uint> nullVec = new List<uint>();
            return nullVec;
        }

        public OpenTK.Vector3d GetVertexCoord(uint vId)
        {
            System.Diagnostics.Debug.Assert(VertexArray.IsObjectId(vId));
            Vertex3D v = VertexArray.GetObject(vId);
            return v.Point;
        }

        public bool GetEdgeVertexId(uint eId, out uint sVId, out uint eVId)
        {
            System.Diagnostics.Debug.Assert(BRep.IsElementId(CadElementType.Edge, eId));
            return BRep.GetEdgeVertexId(eId, out sVId, out eVId);
        }

        public uint GetEdgeVertexId(uint eId, bool isRoot)
        {
            System.Diagnostics.Debug.Assert(BRep.IsElementId(CadElementType.Edge, eId));
            return BRep.GetEdgeVertexId(eId, isRoot);
        }

        public bool GetEdgeLoopId(uint eId, out uint lLId, out uint rLId)
        {
            return BRep.GetEdgeLoopId(eId, out lLId, out rLId);
        }

        public uint GetEdgeLoopId(uint eId, bool isLeft)
        {
            return BRep.GetEdgeLoopId(eId, isLeft);
        }

        public LoopEdgeItr GetLoopEdgeItr(uint lId)
        {
            System.Diagnostics.Debug.Assert(LoopArray.IsObjectId(lId));
            return new LoopEdgeItr(BRep, lId);
        }

        public Loop3D GetLoop(uint lId)
        {
            System.Diagnostics.Debug.Assert(LoopArray.IsObjectId(lId));
            Loop3D l = LoopArray.GetObject(lId);
            l.Edges.Clear();
            l.EdgeIndexs.Clear();

            l.EdgeIndexs.Add(0);
            for (LoopEdgeItr lItr = BRep.GetLoopEdgeItr(lId); !lItr.IsChildEnd; lItr.ShiftChildLoop())
            {
                for (; !lItr.IsEnd(); lItr.Next())
                {
                    uint eId0;
                    bool isSameDir0;
                    lItr.GetEdgeId(out eId0, out isSameDir0);
                    if (eId0 == 0)
                    {
                        continue;
                    }
                    Edge2D e2 = GetEdge2D(eId0, lId);
                    l.Edges.Add(new KeyValuePair<Edge2D, bool>(e2, isSameDir0));
                }
                l.EdgeIndexs.Add((uint)l.Edges.Count);
            }
            l.ClearBoundingBox();
            return l;
        }

        public Edge3D GetEdge(uint eId)
        {
            System.Diagnostics.Debug.Assert(BRep.IsElementId(CadElementType.Edge, eId));
            System.Diagnostics.Debug.Assert(EdgeArray.IsObjectId(eId));
            Edge3D e = EdgeArray.GetObject(eId);
            uint sVId;
            uint eVId;
            BRep.GetEdgeVertexId(eId, out sVId, out eVId);
            e.SetVertexIds(sVId, eVId);
            System.Diagnostics.Debug.Assert(BRep.IsElementId(CadElementType.Vertex, sVId));
            System.Diagnostics.Debug.Assert(BRep.IsElementId(CadElementType.Vertex, eVId));
            System.Diagnostics.Debug.Assert(VertexArray.IsObjectId(sVId));
            System.Diagnostics.Debug.Assert(VertexArray.IsObjectId(eVId));
            e.SetVertexCoords(GetVertexCoord(sVId), GetVertexCoord(eVId));
            return e;
        }

        private Edge2D GetEdge2D(uint eId, uint lId)
        {
            System.Diagnostics.Debug.Assert(IsElementId(CadElementType.Edge, eId));
            System.Diagnostics.Debug.Assert(IsElementId(CadElementType.Loop, lId));
            Loop3D l = LoopArray.GetObject(lId);
            uint sVId;
            uint eVId;
            BRep.GetEdgeVertexId(eId, out sVId, out eVId);
            Edge2D e = new Edge2D(sVId, eVId);
            e.SetVertexCoords(l.Project(GetVertexCoord(sVId)), l.Project(GetVertexCoord(eVId)));
            return e;
        }

        public bool GetCurveAsPolyline(uint eId, out IList<OpenTK.Vector3d> points, double elen = -1)
        {
            points = new List<OpenTK.Vector3d>();

            if (!EdgeArray.IsObjectId(eId))
            {
                return false;
            }
            Edge3D e = GetEdge(eId);
            double len = e.GetCurveLength();
            if (elen > 0)
            {
                uint ndiv = (uint)((len / elen) + 1);
                return e.GetCurveAsPolyline(out points, (int)ndiv);
            }
            return e.GetCurveAsPolyline(out points, -1);
        }

        private int AssertValid()
        {
            if (!BRep.AssertValid())
            {
                return 6;
            }
            {
                IList<uint> lIds = BRep.GetElementIds(CadElementType.Loop);
                for (int i = 0; i < lIds.Count; i++)
                {
                    if (!LoopArray.IsObjectId(lIds[i]))
                    {
                        return 7;
                    }
                }
            }
            {
                IList<uint> eIds = BRep.GetElementIds(CadElementType.Edge);
                for (int i = 0; i < eIds.Count; i++)
                {
                    if (!EdgeArray.IsObjectId(eIds[i]))
                    {
                        return 8;
                    }
                }
            }
            {
                IList<uint> vIds = BRep.GetElementIds(CadElementType.Vertex);
                for (int i = 0; i < vIds.Count; i++)
                {
                    if (!VertexArray.IsObjectId(vIds[i]))
                    {
                        return 9;
                    }
                }
            }
            return 0;
        }
    }
}
