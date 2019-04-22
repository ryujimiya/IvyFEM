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

    public class AddPolygonRes
    {
        public uint AddLId { get; set; } = 0;
        public IList<uint> VIds { get; } = new List<uint>();
        public IList<uint> EIds { get; } = new List<uint>();

        public AddPolygonRes()
        {

        }

        public AddPolygonRes(AddPolygonRes src)
        {
            this.AddLId = src.AddLId;
            this.VIds.Clear();
            foreach (var id in src.VIds)
            {
                this.VIds.Add(id);
            }
            this.EIds.Clear();
            foreach (var id in src.EIds)
            {
                this.EIds.Add(id);
            }
        }
    }

    public class CadObject2D
    {
        protected ObjectArray<Loop2D> LoopArray = new ObjectArray<Loop2D>();
        protected ObjectArray<Edge2D> EdgeArray = new ObjectArray<Edge2D>();
        protected ObjectArray<Vertex2D> VertexArray = new ObjectArray<Vertex2D>();
        protected BRep2D BRep = new BRep2D();
        protected double MinClearance = 1.0e-3;

        public CadObject2D()
        {

        }

        public CadObject2D(CadObject2D src)
        {
            Copy(src);
        }

        public void Clear()
        {
            LoopArray.Clear();
            EdgeArray.Clear();
            VertexArray.Clear();
            BRep.Clear();
        }

        public AddVertexRes AddVertex(CadElementType type, uint id, OpenTK.Vector2d vec)
        {
            AddVertexRes res = new AddVertexRes();
            if (type == CadElementType.NotSet || id == 0)
            {
                uint addVId = BRep.AddVertexToLoop(0);
                uint tmpId = VertexArray.AddObject(addVId, new Vertex2D(vec));
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
                {
                    double dist = SignedDistancePointLoop(lId, vec);
                    if (dist < this.MinClearance) { return res; }
                }
                uint addVId = BRep.AddVertexToLoop(lId);
                uint tmpId = VertexArray.AddObject(addVId, new Vertex2D(vec));
                System.Diagnostics.Debug.Assert(tmpId == (int)addVId);
                System.Diagnostics.Debug.Assert(AssertValid() == 0);
                res.AddVId = addVId;
                return res;
            }
            else if (type == CadElementType.Edge)
            {
                uint eId = id;
                if (!EdgeArray.IsObjectId(eId))
                {
                    return res;
                }
                Edge2D oldEdge = GetEdge(eId);
                OpenTK.Vector2d addVec = oldEdge.GetNearestPoint(vec);
                if (CadUtils.SquareLength(addVec - oldEdge.GetVertex(false)) < 1.0e-20 ||
                    CadUtils.SquareLength(addVec - oldEdge.GetVertex(true)) < 1.0e-20)
                {
                    return res;
                }

                uint addVId = BRep.AddVertexToEdge(eId);
                uint addEId;
                {
                    VertexEdgeItr vItr = BRep.GetVertexEdgeItr(addVId);
                    bool isSameDir0;
                    bool isSameDir1;
                    uint bEId;
                    uint aEId;
                    vItr.GetBehindEdgeId(out bEId, out isSameDir0);
                    vItr.GetAheadEdgeId(out aEId, out isSameDir1);
                    addEId = (bEId == eId) ? aEId : bEId;
                }
                {
                    uint tmpId = VertexArray.AddObject(addVId, new Vertex2D(addVec));
                    System.Diagnostics.Debug.Assert(tmpId == addVId);
                }
                {
                    uint tmpId = EdgeArray.AddObject(addEId, new Edge2D(addVId, oldEdge.GetVertexId(false)));
                    System.Diagnostics.Debug.Assert(tmpId == addEId);
                }
                {
                    Edge2D addEdge = GetEdge(addEId);

                    oldEdge.Split(addEdge, addVec);

                    Edge2D edge = GetEdge(eId);

                    System.Diagnostics.Debug.Assert(edge.GetVertexId(false) == addVId);
                    edge.Copy(oldEdge);
                }
                System.Diagnostics.Debug.Assert(AssertValid() == 0);
                res.AddVId = addVId;
                res.AddEId = addEId;
                return res;
            }
            return res;
        }

        public bool SetCurveLine(uint eId)
        {
            if (!EdgeArray.IsObjectId(eId))
            {
                System.Diagnostics.Debug.Assert(false);
                return false;
            }
            System.Diagnostics.Debug.Assert(EdgeArray.IsObjectId(eId));
            Edge2D e = EdgeArray.GetObject(eId);
            Edge2D oldE = new Edge2D(e);
            ////////////////
            e.SetCurveLine();
            ////////////////
            IList<uint> loopIds = LoopArray.GetObjectIds();
            for (int iLId = 0; iLId < loopIds.Count; iLId++)
            {
                uint id_l = loopIds[iLId];
                if (CheckLoop(id_l) != 0)
                {
                    e.Copy(oldE);
                    return false;
                }
            }
            return true;
        }

        // 円弧の中心は(ePt1 + ePt2) * 0.5 + normal * distanceRatio * edgeLen
        public bool SetCurveArc(uint eId, bool isLeftSide, double distanceRatio)
        {
            if (!EdgeArray.IsObjectId(eId))
            {
                System.Diagnostics.Debug.Assert(false);
                return false;
            }
            System.Diagnostics.Debug.Assert(EdgeArray.IsObjectId(eId));
            Edge2D e = GetEdge(eId);
            Edge2D oldE = new Edge2D(e);
            ////////////////////////////////
            // ここからを現在のCurveTypeによって決める,次の設定は直線の場合
            e.SetCurveArc(isLeftSide, distanceRatio);
            // ここまで
            ////////////////////////////////
            IList<uint> loopIds = LoopArray.GetObjectIds();
            for (int iLId = 0; iLId < loopIds.Count; iLId++)
            {
                uint lId = loopIds[iLId];
                if (CheckLoop(lId) != 0)
                {
                    e.Copy(oldE);
                    return false;
                }
            }
            return true;
        }

        public bool SetCurvePolyline(uint eId)
        {
            if (!EdgeArray.IsObjectId(eId))
            {
                System.Diagnostics.Debug.Assert(false);
                return false;
            }
            System.Diagnostics.Debug.Assert(EdgeArray.IsObjectId(eId));
            Edge2D e = GetEdge(eId);
            Edge2D oldE = new Edge2D(e);
            ////////////////////////////////
            // ここからを現在のCurveTypeによって決める,次の設定は直線の場合
            OpenTK.Vector2d sPt = oldE.GetVertex(true);
            OpenTK.Vector2d ePt = oldE.GetVertex(false);
            IList<OpenTK.Vector2d> pts;
            oldE.GetCurveAsPolyline(out pts, 20);
            double sqLen = CadUtils.SquareLength(ePt - sPt);
            OpenTK.Vector2d hE = (ePt - sPt) * (1 / sqLen);
            OpenTK.Vector2d vE = new OpenTK.Vector2d(-hE.Y, hE.X);
            {
                IList<double> relCos = new List<double>();
                for (int ico = 0; ico < pts.Count; ico++)
                {
                    double x1 = OpenTK.Vector2d.Dot(pts[ico] - sPt, hE);
                    double y1 = OpenTK.Vector2d.Dot(pts[ico] - sPt, vE);
                    relCos.Add(x1);
                    relCos.Add(y1);
                }
                e.SetCurvePolyline(relCos);
            }
            // ここまで
            ////////////////////////////////
            IList<uint> loopIds = LoopArray.GetObjectIds();
            for (int iLId = 0; iLId < loopIds.Count; iLId++)
            {
                uint lId = loopIds[iLId];
                if (CheckLoop(lId) != 0)
                {
                    e.Copy(oldE);
                    return false;
                }
            }
            return true;
        }

        public bool SetCurvePolyline(uint eId, IList<OpenTK.Vector2d> points)
        {
            if (!EdgeArray.IsObjectId(eId))
            {
                System.Diagnostics.Debug.Assert(false);
                return false;
            }
            System.Diagnostics.Debug.Assert(EdgeArray.IsObjectId(eId));
            Edge2D e = GetEdge(eId);
            Edge2D oldE = new Edge2D(e);
            ////////////////
            {
                // 相対座標を作る    
                int n = points.Count;
                IList<double> relCos = new List<double>();
                OpenTK.Vector2d sPt = e.GetVertex(true);
                OpenTK.Vector2d ePt = e.GetVertex(false);
                double sqlen = CadUtils.SquareLength(ePt - sPt);
                OpenTK.Vector2d hE = (ePt - sPt) * (1 / sqlen);
                OpenTK.Vector2d vE = new OpenTK.Vector2d(-hE.Y, hE.X);
                for (int i = 0; i < n; i++)
                {
                    double x0 = OpenTK.Vector2d.Dot(points[i] - sPt, hE);
                    double y0 = OpenTK.Vector2d.Dot(points[i] - sPt, vE);
                    relCos.Add(x0);
                    relCos.Add(y0);
                }
                System.Diagnostics.Debug.Assert(relCos.Count == n * 2);
                e.SetCurvePolyline(relCos);
            }
            ////////////////
            IList<uint> loopIds = LoopArray.GetObjectIds();
            for (int iLId = 0; iLId < loopIds.Count; iLId++)
            {
                uint lId = loopIds[iLId];
                if (CheckLoop(lId) != 0)
                {
                    e.Copy(oldE);
                    return false;
                }
            }
            return true;
        }

        public bool SetCurveBezier(uint eId, double cx0, double cy0, double cx1, double cy1)
        {
            if (!EdgeArray.IsObjectId(eId))
            {
                return false;
            }
            System.Diagnostics.Debug.Assert(EdgeArray.IsObjectId(eId));
            Edge2D e = EdgeArray.GetObject(eId);
            Edge2D oldE = new Edge2D(e);
            ////////////////
            e.SetCurveBezier(cx0, cy0, cx1, cy1);
            ////////////////
            IList<uint> loopIds = LoopArray.GetObjectIds();
            for (int iLId = 0; iLId < loopIds.Count; iLId++)
            {
                uint lId = loopIds[iLId];
                if (CheckLoop(lId) != 0)
                {
                    e.Copy(oldE);
                    return false;
                }
            }
            return true;
        }

        /*
        public ResAddPolygon AddLoop(IList<KeyValuePair<CurveType, IList<double>>> points,  uint lId, double scale)
        {
            ResAddPolygon res = new ResAddPolygon();

            try
            {
                int ptCnt = points.Count;
                IList<OpenTK.Vector2d> vecs = new List<OpenTK.Vector2d>();
                //for (int iPt = 0; iPt < ptCnt - 1; iPt++)
                for (int iPt = 0; iPt < ptCnt; iPt++)
                {
                    IList<double> point = points[iPt].Value;
                    OpenTK.Vector2d vec = new OpenTK.Vector2d(
                        point[point.Count - 2] * scale, point[point.Count - 1] * scale);
                    vecs.Add(vec);
                    uint vId0 = AddVertex(CadElementType.Loop, lId, vec).AddVId;
                    if (vId0 == 0)
                    {
                        throw new InvalidOperationException("FAIL_ADD_POLYGON_INSIDE_LOOP");
                    }
                    res.VIds.Add(vId0);
                }
                System.Diagnostics.Debug.Assert(res.VIds.Count == ptCnt);
                //for (int iEdge = 0; iEdge < ptCnt - 1; iEdge++)
                for (int iEdge = 0; iEdge < ptCnt; iEdge++)
                {
                    int isPt = iEdge;
                    //int iePt = (iEdge != ptCnt - 2) ? iEdge + 1 : 0;
                    int iePt = (iEdge != ptCnt - 1) ? iEdge + 1 : 0;
                    Edge2D e = new Edge2D(res.VIds[isPt], res.VIds[iePt]);
                    //System.Diagnostics.Debug.WriteLine(
                    //    iEdge + " " + ptCnt + "    " + res.VIds[isPt] + " " + res.VIds[iePt] + "  ");
                    if (iEdge != ptCnt - 1)
                    {
                        CurveType type = points[iEdge + 1].Key;
                        if (type == CurveType.CurveBezier)
                        {
                            IList<double> point = points[iEdge + 1].Value;
                            OpenTK.Vector2d vs = vecs[isPt];
                            OpenTK.Vector2d ve = vecs[iePt];
                            OpenTK.Vector2d vsc = new OpenTK.Vector2d(point[0] * scale, point[1] * scale);
                            OpenTK.Vector2d vec = new OpenTK.Vector2d(point[2] * scale, point[3] * scale);
                            OpenTK.Vector2d vh = ve - vs;
                            {
                                double len = vh.Length;
                                vh *= 1.0 / (len * len);
                            }
                            OpenTK.Vector2d vv = new OpenTK.Vector2d(-vh.Y, vh.X);
                            double t0 = OpenTK.Vector2d.Dot(vsc - vs, vh);
                            double t1 = OpenTK.Vector2d.Dot(vsc - vs, vv);
                            double t2 = OpenTK.Vector2d.Dot(vec - vs, vh);
                            double t3 = OpenTK.Vector2d.Dot(vec - vs, vv);
                            //System.Diagnostics.Debug.WriteLine(t0 + " " + t1 + " " + t2 + " " + t3);
                            e.SetCurveBezier(t0, t1, t2, t3);
                        }
                    }
                    uint eId0 = ConnectVertex(e).AddEId;
                    //System.Diagnostics.Debug.WriteLine("edge add " + eId0);
                    if (eId0 == 0)
                    {
                        throw new InvalidOperationException("FAIL_ADD_POLYGON_INSIDE_LOOP");
                    }
                    res.EIds.Add(eId0);
                }
                System.Diagnostics.Debug.Assert(res.EIds.Count == ptCnt);
                System.Diagnostics.Debug.Assert(AssertValid() == 0);
                // 新しく出来たループのIDを取得  
                {
                    // 辺の両側のループを調べる
                    uint eId0 = res.EIds[ptCnt - 1];
                    uint lId0;
                    uint lId1;
                    BRep.GetEdgeLoopId(eId0, out lId0, out lId1);
                    res.AddLId = (lId0 == lId) ? lId1 : lId0;
                }

            }
            catch (InvalidOperationException exception)
            {
                for (int iie = 0; iie < res.EIds.Count; iie++)
                {
                    uint id_e0 = res.EIds[iie];

                    RemoveElement(CadElementType.Edge, id_e0);

                }
                for (int iiv = 0; iiv < res.VIds.Count; iiv++)
                {
                    uint id_v0 = res.VIds[iiv];

                    RemoveElement(CadElementType.Vertex, id_v0);

                }
                System.Diagnostics.Debug.Assert(AssertValid() == 0);
                return new ResAddPolygon();
            }
            return res;
        }
        */

        public AddPolygonRes AddPolygon(IList<OpenTK.Vector2d> points, uint lId = 0)
        {
            AddPolygonRes res = new AddPolygonRes();

            int ptCnt = points.Count;
            if (ptCnt < 3)
            {
                return res;
            }

            try
            {
                IList<OpenTK.Vector2d> points1 = new List<OpenTK.Vector2d>(points);
                {
                    uint n = (uint)points.Count;
                    IList<Edge2D> edges = new List<Edge2D>();
                    for (uint i = 0; i < n - 1; i++)
                    {
                        Edge2D e = new Edge2D(i, i + 1);
                        e.SetVertexCoords(points[(int)i], points[(int)i + 1]);
                        edges.Add(e);
                    }
                    {
                        Edge2D e = new Edge2D(n - 1, 0);
                        e.SetVertexCoords(points[(int)n - 1], points[0]);
                        edges.Add(e);
                    }
                    if (CadUtils.CheckEdgeIntersection(edges) != 0)
                    {
                        return res;
                    }
                }
                for (uint i = 0; i < ptCnt; i++)
                {
                    uint vId0 = AddVertex(CadElementType.Loop, lId, points1[(int)i]).AddVId;
                    if (vId0 == 0)
                    {
                        throw new InvalidOperationException("FAIL_ADD_POLYGON_INSIDE_LOOP");
                    }
                    res.VIds.Add(vId0);
                }
                for (uint iEdge = 0; iEdge < ptCnt - 1; iEdge++)
                {
                    uint eId0 = ConnectVertexLine(res.VIds[(int)iEdge], res.VIds[(int)iEdge + 1]).AddEId;
                    if (eId0 == 0)
                    {
                        throw new InvalidOperationException("FAIL_ADD_POLYGON_INSIDE_LOOP");
                    }
                    res.EIds.Add(eId0);
                }
                {
                    uint eId0 = ConnectVertexLine(res.VIds[ptCnt - 1], res.VIds[0]).AddEId;
                    if (eId0 == 0)
                    {
                        throw new InvalidOperationException("FAIL_ADD_POLYGON_INSIDE_LOOP");
                    }
                    res.EIds.Add(eId0);
                }

                System.Diagnostics.Debug.Assert(AssertValid() == 0);

                {
                    uint eId0 = res.EIds[ptCnt - 1];
                    uint lId0 = 0;
                    uint lId1 = 0;
                    BRep.GetEdgeLoopId(eId0, out lId0, out lId1);
                    res.AddLId = (lId0 == lId) ? lId1 : lId0;
                }
                return res;
            }
            catch (InvalidOperationException exception)
            {
                System.Diagnostics.Debug.WriteLine(exception.ToString());
                for (uint iEId = 0; iEId < res.EIds.Count; iEId++)
                {
                    uint eId0 = res.EIds[(int)iEId];
                    RemoveElement(CadElementType.Edge, eId0);
                }
                for (uint iVId = 0; iVId < res.VIds.Count; iVId++)
                {
                    uint vId0 = res.VIds[(int)iVId];
                    RemoveElement(CadElementType.Vertex, vId0);
                }
                System.Diagnostics.Debug.Assert(AssertValid() == 0);
            }

            //　失敗したとき
            return new AddPolygonRes();

        }

        public AddPolygonRes AddCircle(OpenTK.Vector2d cPt, double r, uint lId)
        {
            AddPolygonRes res = new AddPolygonRes();
            OpenTK.Vector2d v1 = new OpenTK.Vector2d(cPt.X, cPt.Y + r);
            OpenTK.Vector2d v2 = new OpenTK.Vector2d(cPt.X - r, cPt.Y);
            OpenTK.Vector2d v3 = new OpenTK.Vector2d(cPt.X, cPt.Y - r);
            OpenTK.Vector2d v4 = new OpenTK.Vector2d(cPt.X + r, cPt.Y);
            var resV1 = AddVertex(CadElementType.Loop, lId, v1);
            var resV2 = AddVertex(CadElementType.Loop, lId, v2);
            var resV3 = AddVertex(CadElementType.Loop, lId, v3);
            var resV4 = AddVertex(CadElementType.Loop, lId, v4);
            var resE1 = ConnectVertexLine(resV1.AddVId, resV2.AddVId);
            var resE2 = ConnectVertexLine(resV2.AddVId, resV3.AddVId);
            var resE3 = ConnectVertexLine(resV3.AddVId, resV4.AddVId);
            var resE4 = ConnectVertexLine(resV4.AddVId, resV1.AddVId);
            double distanceRatio = 0.5;
            bool success1 = SetCurveArc(resE1.AddEId, false, distanceRatio);
            System.Diagnostics.Debug.Assert(success1);
            bool success2 = SetCurveArc(resE2.AddEId, false, distanceRatio);
            System.Diagnostics.Debug.Assert(success2);
            bool success3 = SetCurveArc(resE3.AddEId, false, distanceRatio);
            System.Diagnostics.Debug.Assert(success3);
            bool success4 = SetCurveArc(resE4.AddEId, false, distanceRatio);
            System.Diagnostics.Debug.Assert(success4);
            res.AddLId = resE4.AddLId;
            res.EIds.Add(resE1.AddEId);
            res.EIds.Add(resE2.AddEId);
            res.EIds.Add(resE3.AddEId);
            res.EIds.Add(resE4.AddEId);
            res.VIds.Add(resV1.AddVId);
            res.VIds.Add(resV2.AddVId);
            res.VIds.Add(resV3.AddVId);
            res.VIds.Add(resV4.AddVId);
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

        public Vertex2D GetVertex(uint vId)
        {
            if (!BRep.IsElementId(CadElementType.Vertex, vId))
            {
                return null;
            }
            if (!VertexArray.IsObjectId(vId))
            {
                return null;
            }
            Vertex2D v = VertexArray.GetObject(vId);
            return v;
        }

        public OpenTK.Vector2d GetVertexCoord(uint vId)
        {
            System.Diagnostics.Debug.Assert(VertexArray.IsObjectId(vId));
            Vertex2D v = VertexArray.GetObject(vId);
            return v.Point;
        }

        public double[] GetVertexColor(uint vId)
        {
            double[] color = new double[3];
            if (!BRep.IsElementId(CadElementType.Vertex, vId))
            {
                return color;
            }
            if (!VertexArray.IsObjectId(vId))
            {
                return color;
            }
            Vertex2D v = VertexArray.GetObject(vId);
            v.Color.CopyTo(color, 0);
            return color;
        }

        public bool SetVertexColor(uint vId, double[] color)
        {
            if (!BRep.IsElementId(CadElementType.Vertex, vId))
            {
                return false;
            }
            if (!VertexArray.IsObjectId(vId))
            {
                return false;
            }
            Vertex2D v = VertexArray.GetObject(vId);
            color.CopyTo(v.Color, 0);
            return true;
        }

        public Edge2D GetEdge(uint eId)
        {
            System.Diagnostics.Debug.Assert(BRep.IsElementId(CadElementType.Edge, eId));
            System.Diagnostics.Debug.Assert(EdgeArray.IsObjectId(eId));
            Edge2D e = EdgeArray.GetObject(eId);
            uint sVId;
            uint eVId;
            BRep.GetEdgeVertexIds(eId, out sVId, out eVId);
            e.SetVertexIds(sVId, eVId);
            System.Diagnostics.Debug.Assert(BRep.IsElementId(CadElementType.Vertex, sVId));
            System.Diagnostics.Debug.Assert(BRep.IsElementId(CadElementType.Vertex, eVId));
            System.Diagnostics.Debug.Assert(VertexArray.IsObjectId(sVId));
            System.Diagnostics.Debug.Assert(VertexArray.IsObjectId(eVId));
            e.SetVertexCoords(GetVertexCoord(sVId), GetVertexCoord(eVId));
            return e;
        }

        public double[] GetEdgeColor(uint eId)
        {
            double[] color = new double[3];
            if (!BRep.IsElementId(CadElementType.Edge, eId))
            {
                return color;
            }
            if (!EdgeArray.IsObjectId(eId))
            {
                return color;
            }
            Edge2D e = EdgeArray.GetObject(eId);
            e.Color.CopyTo(color, 0);
            return color;
        }

        public bool SetEdgeColor(uint eId, double[] color)
        {
            if (!BRep.IsElementId(CadElementType.Edge, eId))
            {
                return false;
            }
            if (!EdgeArray.IsObjectId(eId))
            {
                return false;
            }
            Edge2D e = EdgeArray.GetObject(eId);
            color.CopyTo(e.Color, 0);
            return true;
        }

        public bool GetEdgeVertexId(out uint sVId, out uint eVId, uint eId)
        {
            System.Diagnostics.Debug.Assert(BRep.IsElementId(CadElementType.Edge, eId));
            return BRep.GetEdgeVertexIds(eId, out sVId, out eVId);
        }

        public uint GetEdgeVertexId(uint eId, bool isS)
        {
            System.Diagnostics.Debug.Assert(BRep.IsElementId(CadElementType.Edge, eId));
            return BRep.GetEdgeVertexId(eId, isS);
        }

        public bool GetEdgeLoopId(out uint lLId, out uint rLId, uint eId)
        {
            return BRep.GetEdgeLoopId(eId, out lLId, out rLId);
        }

        public CurveType GetEdgeCurveType(uint eId)
        {
            System.Diagnostics.Debug.Assert(EdgeArray.IsObjectId(eId));
            Edge2D e = EdgeArray.GetObject(eId);
            return e.CurveType;
        }

        public bool GetCurveAsPolyline(uint eId, out IList<OpenTK.Vector2d> points, double elen = -1)
        {
            points = new List<OpenTK.Vector2d>();

            if (!EdgeArray.IsObjectId(eId))
            {
                return false;
            }
            Edge2D e = GetEdge(eId);
            double len = e.GetCurveLength();
            if (elen > 0)
            {
                uint ndiv = (uint)((len / elen) + 1);
                return e.GetCurveAsPolyline(out points, (int)ndiv);
            }
            return e.GetCurveAsPolyline(out points, -1);
        }

        public int GetLayer(CadElementType type, uint id)
        {
            if (type == CadElementType.Loop)
            {
                if (!LoopArray.IsObjectId(id))
                {
                    return 0;
                }
                Loop2D l = LoopArray.GetObject(id);
                return (int)l.Layer;
            }
            else if (type == CadElementType.Edge)
            {
                uint lLId;
                uint rLId;
                GetEdgeLoopId(out lLId, out rLId, id);

                bool bl = IsElementId(CadElementType.Loop, lLId);
                bool br = IsElementId(CadElementType.Loop, rLId);
                if (!bl && !br) { return 0; }
                if (bl && !br) { return GetLayer(CadElementType.Loop, lLId); }
                if (!bl && br) { return GetLayer(CadElementType.Loop, rLId); }
                int ilayer_l = GetLayer(CadElementType.Loop, lLId);
                int ilayer_r = GetLayer(CadElementType.Loop, rLId);
                return (ilayer_l > ilayer_r) ? ilayer_l : ilayer_r;
            }
            else if (type == CadElementType.Vertex)
            {
                int layer = 0;
                bool iflg = true;
                for (VertexEdgeItr vItr = BRep.GetVertexEdgeItr(id); !vItr.IsEnd(); vItr++)
                {
                    uint lId0 = vItr.GetLoopId();
                    if (!IsElementId(CadElementType.Loop, lId0))
                    {
                        continue;
                    }
                    int layer0 = GetLayer(CadElementType.Loop, lId0);
                    if (iflg == true)
                    {
                        layer = layer0;
                        iflg = false;
                    }
                    else
                    {
                        layer = (layer0 > layer) ? layer0 : layer;
                    }
                }
                return layer;
            }
            return 0;
        }

        public void GetLayerMinMax(out int minLayer, out int maxLayer)
        {
            IList<uint> lIds = GetElementIds(CadElementType.Loop);
            if (lIds.Count == 0)
            {
                minLayer = 0;
                maxLayer = 0;
                return;
            }

            {
                System.Diagnostics.Debug.Assert(lIds.Count > 0);
                uint lId0 = lIds[0];
                minLayer = GetLayer(CadElementType.Loop, lId0);
                maxLayer = minLayer;
            }
            for (int i = 0; i < lIds.Count; i++)
            {
                uint lId = lIds[i];
                int layer = GetLayer(CadElementType.Loop, lId);
                minLayer = (layer < minLayer) ? layer : minLayer;
                maxLayer = (layer > maxLayer) ? layer : maxLayer;
            }
        }

        public double[] GetLoopColor(uint lId)
        {
            double[] color = new double[3];
            if (!LoopArray.IsObjectId(lId))
            {
                return color;
            }
            Loop2D l = LoopArray.GetObject(lId);
            l.Color.CopyTo(color, 0);
            return color;
        }

        public bool SetLoopColor(uint lId, double[] color)
        {
            if (!LoopArray.IsObjectId(lId))
            {
                return false;
            }
            Loop2D l = LoopArray.GetObject(lId);
            color.CopyTo(l.Color, 0);
            return true;
        }

        private int AssertValid()
        {
            {
                IList<uint> lIds = LoopArray.GetObjectIds();
                for (uint i = 0; i < lIds.Count; i++)
                {
                    uint lId = lIds[(int)i];
                    int res = CheckLoop(lId);
                    if (res != 0)
                    {
                        if (res == 1)
                        {
                            System.Diagnostics.Debug.WriteLine("Intersectoin in the loop");
                        }
                        else if (res == 2)
                        {
                            System.Diagnostics.Debug.WriteLine("Check area parent plus, childe minus");
                        }
                        else if (res == 3)
                        {
                            System.Diagnostics.Debug.WriteLine("Check whether childe loop included in parent loop");
                        }
                        else if (res == 4)
                        {
                            System.Diagnostics.Debug.WriteLine("Check childe loop excluded from other child loop");
                        }
                        else if (res == 5)
                        {
                            System.Diagnostics.Debug.WriteLine("Check positive angle around vertex on the loop");
                        }
                        return res;
                    }
                }
            }
            if (!BRep.AssertValid())
            {
                return 6;
            }
            {
                IList<uint> lIds = BRep.GetElementIds(CadElementType.Loop);
                for (uint i = 0; i < lIds.Count; i++)
                {
                    if (!LoopArray.IsObjectId(lIds[(int)i]))
                    {
                        //System.Diagnostics.Debug.WriteLine(lIds[(int)i]);
                        return 7;
                    }
                }
            }
            {
                IList<uint> eIds = BRep.GetElementIds(CadElementType.Edge);
                for (uint i = 0; i < eIds.Count; i++)
                {
                    if (!EdgeArray.IsObjectId(eIds[(int)i]))
                    {
                        return 7;
                    }
                }
            }
            {
                IList<uint> vIds = BRep.GetElementIds(CadElementType.Vertex);
                for (uint i = 0; i < vIds.Count; i++)
                {
                    if (!VertexArray.IsObjectId(vIds[(int)i]))
                    {
                        return 7;
                    }
                }
            }
            return 0;
        }

        private bool CheckIsPointInsideLoopEdgeItr(LoopEdgeItr lItr, OpenTK.Vector2d point)
        {
            // 29 is handy prime number
            for (uint i = 1; i < 29; i++)
            {
                uint crossCounter = 0;
                bool iflg = true;
                OpenTK.Vector2d dir = new OpenTK.Vector2d(Math.Sin(6.28 * i / 29.0), Math.Cos(6.28 * i / 29.0));
                for (lItr.Begin(); !lItr.IsEnd(); lItr++)
                {
                    uint eId;
                    bool isSameDir;
                    lItr.GetEdgeId(out eId, out isSameDir);
                    if (eId == 0)
                    {
                        return false;
                    }
                    System.Diagnostics.Debug.Assert(EdgeArray.IsObjectId(eId));
                    Edge2D e = GetEdge(eId);
                    int ires = e.NumIntersectAgainstHalfLine(point, dir);
                    if (ires == -1)
                    {
                        iflg = false;
                        break;
                    }
                    crossCounter += (uint)ires;
                }
                if (iflg == true)
                {
                    if (crossCounter % 2 == 0) return false;
                    return true;
                }
            }
            System.Diagnostics.Debug.Assert(false);
            return false;
        }

        private uint CheckLoopEdgeItrPointInOutLoopEdgeItr(LoopEdgeItr lItr1, LoopEdgeItr lItr2)
        {
            uint outCount = 0;
            uint inCount = 0;
            for (lItr1.Begin(); !lItr1.IsEnd(); lItr1++)
            {
                uint vId = lItr1.GetVertexId();
                Vertex2D v = VertexArray.GetObject(vId);
                double dist = DistancePointLoopEdgeItr(lItr2, v.Point);
                if (Math.Abs(dist) < MinClearance)
                {
                    return 1;
                } 
                if (CheckIsPointInsideLoopEdgeItr(lItr2, v.Point))
                {
                    if (outCount != 0)
                    {
                        return 1;
                    }
                    inCount++;
                }
                else
                {
                    if (inCount != 0)
                    {
                        return 1;
                    }
                    outCount++;
                }
            }
            if (inCount == 0)
            {
                System.Diagnostics.Debug.Assert(outCount != 0);
                return 2;
            }
            System.Diagnostics.Debug.Assert(outCount == 0);
            System.Diagnostics.Debug.Assert(inCount != 0);
            return 0;
        }

        public bool CheckIsPointInsideLoop(uint lId1, OpenTK.Vector2d point)
        {
            System.Diagnostics.Debug.Assert(LoopArray.IsObjectId(lId1));
            for (LoopEdgeItr lItr = BRep.GetLoopEdgeItr(lId1); !lItr.IsChildEnd; lItr.ShiftChildLoop())
            {
                if (lItr.IsParent())
                {
                    if (!CheckIsPointInsideLoopEdgeItr(lItr, point))
                    {
                        return false;
                    }
                }
                else
                {
                    if (CheckIsPointInsideLoopEdgeItr(lItr, point))
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        public double GetLoopArea(uint lId)
        {
            System.Diagnostics.Debug.Assert(LoopArray.IsObjectId(lId));
            double area = 0.0;
            for (LoopEdgeItr lItr = BRep.GetLoopEdgeItr(lId); !lItr.IsChildEnd; lItr.ShiftChildLoop())
            {
                area += GetLoopEdgeItrArea(lItr);
            }
            return area;
        }

        private double GetLoopEdgeItrArea(LoopEdgeItr lItr)
        {
            double area = 0.0;
            for (lItr.Begin(); !lItr.IsEnd(); lItr++)
            {
                uint eId;
                bool isSameDir;
                lItr.GetEdgeId(out eId, out isSameDir);
                if (!IsElementId(CadElementType.Edge, eId))
                {
                    return 0;
                }
                System.Diagnostics.Debug.Assert(IsElementId(CadElementType.Edge, eId));
                Edge2D e = GetEdge(eId);
                System.Diagnostics.Debug.Assert(e.GetVertexId(isSameDir) == lItr.GetVertexId());
                System.Diagnostics.Debug.Assert(e.GetVertexId(!isSameDir) == lItr.GetAheadVertexId());
                double earea = CadUtils.TriArea(
                    e.GetVertex(true), e.GetVertex(false), new OpenTK.Vector2d(0, 0)) + e.EdgeArea();
                if (isSameDir)
                {
                    area += earea;
                }
                else
                {
                    area -= earea;
                }
            }
            return area;
        }

        private double DistancePointLoopEdgeItr(LoopEdgeItr lItr, OpenTK.Vector2d point)
        {
            double minDist = -1;
            for (lItr.Begin(); !lItr.IsEnd(); lItr++)
            {
                uint eId;
                bool isSameDir;
                lItr.GetEdgeId(out eId, out isSameDir);
                if (eId == 0)
                {
                    uint vId0 = lItr.GetVertexId();
                    System.Diagnostics.Debug.Assert(IsElementId(CadElementType.Vertex, vId0));
                    OpenTK.Vector2d p1 = GetVertexCoord(vId0);
                    return OpenTK.Vector2d.Distance(point, p1);
                }
                System.Diagnostics.Debug.Assert(EdgeArray.IsObjectId(eId));
                Edge2D e = GetEdge(eId);
                OpenTK.Vector2d v = e.GetNearestPoint(point);
                double d0 = OpenTK.Vector2d.Distance(v, point);
                if (minDist < 0 || d0 < minDist)
                {
                    minDist = d0;
                }
            }
            return minDist;
        }

        public double SignedDistancePointLoop(uint lId1, OpenTK.Vector2d point, uint ignoreVId = 0)
        {
            double minSd = 0;
            System.Diagnostics.Debug.Assert(LoopArray.IsObjectId(lId1));
            for (LoopEdgeItr lItr = BRep.GetLoopEdgeItr(lId1); !lItr.IsChildEnd; lItr.ShiftChildLoop())
            {
                if (lItr.IsParent())
                {
                    minSd = DistancePointLoopEdgeItr(lItr, point);
                    System.Diagnostics.Debug.Assert(minSd >= 0);
                    if (!CheckIsPointInsideLoopEdgeItr(lItr, point))
                    {
                        minSd = -minSd;
                    }
                }
                else
                {
                    if (lItr.GetVertexId() == lItr.GetAheadVertexId())
                    {
                        uint vId = lItr.GetVertexId();
                        if (vId == ignoreVId)
                        {
                            continue;
                        }
                    }
                    double sd0 = DistancePointLoopEdgeItr(lItr, point);
                    if (sd0 < 0)
                    {
                        continue;
                    }
                    if (CheckIsPointInsideLoopEdgeItr(lItr, point))
                    {
                        sd0 = -sd0;
                    }
                    if (Math.Abs(sd0) < Math.Abs(minSd))
                    {
                        minSd = sd0;
                    }
                }
            }
            return minSd;
        }

        public ConnectVertexRes ConnectVertexLine(uint vId1, uint vId2)
        {
            Edge2D e = new Edge2D(vId1, vId2);
            return ConnectVertex(e);
        }

        public ConnectVertexRes ConnectVertex(Edge2D edge)
        {
            uint vId1 = edge.GetVertexId(true);
            uint vId2 = edge.GetVertexId(false);
            ConnectVertexRes res = new ConnectVertexRes();
            res.VId1 = vId1;
            res.VId2 = vId2;
            ////  
            if (!VertexArray.IsObjectId(vId1))
            {
                return res;
            }
            if (!VertexArray.IsObjectId(vId2))
            {
                return res;
            }
            if (vId1 == vId2)
            {
                return res;
            }

            if (edge.CurveType == CurveType.CurveLine)
            {
                IList<uint> eIds = EdgeArray.GetObjectIds();
                for (uint i = 0; i < eIds.Count; i++)
                {
                    uint eId = eIds[(int)i];
                    System.Diagnostics.Debug.Assert(EdgeArray.IsObjectId(eId));
                    Edge2D e = GetEdge(eId);
                    if (e.CurveType != CurveType.CurveLine)
                    {
                        continue;
                    }
                    uint sVId = e.GetVertexId(true);
                    uint eVId = e.GetVertexId(false);
                    if (sVId != vId1 && sVId != vId2)
                    {
                        continue;
                    }
                    if (eVId != vId1 && eVId != vId2)
                    {
                        continue;
                    }
                    return res;
                }
            }
            edge.SetVertexCoords(VertexArray.GetObject(vId1).Point, VertexArray.GetObject(vId2).Point);
            if (edge.IsCrossEdgeSelf())
            {
                return res;
            }

            {
                VertexEdgeItr vItr1 = FindCornerHalfLine(vId1, edge.GetTangentEdge(true));
                VertexEdgeItr vItr2 = FindCornerHalfLine(vId2, edge.GetTangentEdge(false));
                if (vItr1.GetLoopId() != vItr2.GetLoopId())
                {
                    return res;
                }
                uint lId = vItr1.GetLoopId();
                if (CheckEdgeAgainstLoopIntersection(edge, lId))
                {
                    return res;
                }
                bool isLeftAddL = false;
                if (vItr1.IsSameUseLoop(vItr2) && (!vItr1.IsParent() || lId == 0))
                {
                    IList<KeyValuePair<uint, bool>> eId2Dir = BRep.GetEdgesForConnectVertex(vItr1, vItr2);
                    System.Diagnostics.Debug.Assert(eId2Dir.Count != 0);
                    double area = CadUtils.TriArea(
                        edge.GetVertex(true), edge.GetVertex(false), new OpenTK.Vector2d(0, 0)) + edge.EdgeArea();
                    for (uint i = 0; i < eId2Dir.Count; i++)
                    {
                        uint eId = eId2Dir[(int)i].Key;
                        System.Diagnostics.Debug.Assert(IsElementId(CadElementType.Edge, eId));
                        Edge2D e = GetEdge(eId);
                        double earea = e.EdgeArea() +
                            CadUtils.TriArea(e.GetVertex(true), e.GetVertex(false), new OpenTK.Vector2d(0, 0));
                        if (eId2Dir[(int)i].Value)
                        {
                            area += earea;
                        }
                        else
                        {
                            area -= earea;
                        }
                    }
                    isLeftAddL = (area > 0);
                }
                res = BRep.ConnectVertex(vItr1, vItr2, isLeftAddL);
            }
            {
                uint tmpId = EdgeArray.AddObject(res.AddEId, edge);
                System.Diagnostics.Debug.Assert(tmpId == (int)res.AddEId);
            }

            {
                uint lLId;
                uint rLId;
                BRep.GetEdgeLoopId(res.AddEId, out lLId, out rLId);
                if (lLId == rLId)
                {
                    System.Diagnostics.Debug.Assert(AssertValid() == 0);
                    return res;
                }
                res.AddLId = (res.IsLeftAddL) ? lLId : rLId;
                System.Diagnostics.Debug.Assert(res.AddLId != res.LId || res.LId == 0);
                System.Diagnostics.Debug.Assert(((res.IsLeftAddL) ? lLId : rLId) == res.AddLId);
            }

            if (!BRep.IsElementId(CadElementType.Loop, res.LId))
            {
                if (CheckLoopIntersection(res.AddLId))
                {
                    BRep.MakeHoleFromLoop(res.AddLId);
                    System.Diagnostics.Debug.Assert(AssertValid() == 0);
                    return res;
                }
            }

            if (BRep.IsElementId(CadElementType.Loop, res.LId) && BRep.IsElementId(CadElementType.Loop, res.AddLId))
            {
                while (true)
                {
                    bool iflg = true;
                    LoopEdgeItr lItrAddInner = BRep.GetLoopEdgeItrSideEdge(res.AddEId, res.IsLeftAddL);
                    LoopEdgeItr lItrAddOuter = BRep.GetLoopEdgeItrSideEdge(res.AddEId, !res.IsLeftAddL);
                    for (LoopEdgeItr cLItr = BRep.GetLoopEdgeItr(res.LId); !cLItr.IsChildEnd; cLItr.ShiftChildLoop())
                    {
                        if (cLItr.IsParent())
                        {
                            continue;
                        }
                        if (cLItr.IsSameUseLoop(lItrAddOuter))
                        {
                            continue;
                        }
                        uint ires = CheckLoopEdgeItrPointInOutLoopEdgeItr(cLItr, lItrAddInner);
                        System.Diagnostics.Debug.Assert(ires == 0 || ires == 2);
                        if (ires == 0)
                        {
                            BRep.SwapLoopEdgeItr(cLItr, res.AddLId);
                            iflg = false;
                            break;
                        }
                    }
                    if (iflg)
                    {
                        break;
                    }
                }
            }

            if (LoopArray.IsObjectId(res.LId))
            {
                Loop2D addLoop = LoopArray.GetObject(res.LId);
                // 2019-03-20 BUGFIX
                addLoop = new Loop2D(addLoop);
                LoopArray.AddObject(res.AddLId, addLoop);
            }
            else
            {
                LoopArray.AddObject(res.AddLId, new Loop2D());
            }

            System.Diagnostics.Debug.Assert(AssertValid() == 0);
            return res;
        }

        private VertexEdgeItr FindCornerHalfLine(uint vId, OpenTK.Vector2d dir1) 
        {
            System.Diagnostics.Debug.Assert(VertexArray.IsObjectId(vId));
            OpenTK.Vector2d dir = dir1;
            dir = CadUtils.Normalize(dir);
            OpenTK.Vector2d zeroVec = new OpenTK.Vector2d(0, 0);
            VertexEdgeItr vItr = BRep.GetVertexEdgeItr(vId);
            if (vItr.GetEdgeCount() < 2)
            {
                return vItr;
            }
            for (; !vItr.IsEnd(); vItr++)
            {
                uint eId0;
                bool isSameDir0;
                vItr.GetBehindEdgeId(out eId0, out isSameDir0);
                System.Diagnostics.Debug.Assert(EdgeArray.IsObjectId(eId0));
                Edge2D e0 = GetEdge(eId0);
                System.Diagnostics.Debug.Assert(e0.GetVertexId(isSameDir0) == vId);
                OpenTK.Vector2d tan0 = e0.GetTangentEdge(isSameDir0);

                uint eId1;
                bool isSameDir1;
                vItr.GetAheadEdgeId(out eId1, out isSameDir1);
                System.Diagnostics.Debug.Assert(EdgeArray.IsObjectId(eId1));
                Edge2D e1 = GetEdge(eId1);
                System.Diagnostics.Debug.Assert(e1.GetVertexId(isSameDir1) == vId);
                OpenTK.Vector2d tan1 = e1.GetTangentEdge(isSameDir1);
                System.Diagnostics.Debug.Assert(eId0 != eId1);

                double area0 = CadUtils.TriArea(tan1, zeroVec, tan0);
                double area1 = CadUtils.TriArea(tan1, zeroVec, dir);
                double area2 = CadUtils.TriArea(dir, zeroVec, tan0);
                if (area0 > 0.0)
                {
                    if (area1 > 0.0 && area2 > 0.0)
                    {
                        return vItr;
                    }
                }
                else
                {
                    if (area1 > 0.0 || area2 > 0.0)
                    {
                        return vItr;
                    }
                }
            }
            return vItr;
        }

        public bool RemoveElement(CadElementType type, uint id)
        {
            if (!IsElementId(type, id))
            {
                return false;
            }
            if (type == CadElementType.Edge)
            {
                LoopEdgeItr lItrL = BRep.GetLoopEdgeItrSideEdge(id, true);
                LoopEdgeItr lItrR = BRep.GetLoopEdgeItrSideEdge(id, false);
                uint lLId = lItrL.GetLoopId();
                uint rLId = lItrR.GetLoopId();
                uint vId1;
                uint vId2;
                BRep.GetEdgeVertexIds(id, out vId1, out vId2);
                VertexEdgeItr vItr1 = BRep.GetVertexEdgeItr(vId1);
                VertexEdgeItr vItr2 = BRep.GetVertexEdgeItr(vId2);
                bool isDelCP = false;
                if (lItrL.IsSameUseLoop(lItrR) && lItrL.IsParent() && lItrL.GetLoopId() != 0 &&
                    vItr1.GetEdgeCount() > 1 && vItr2.GetEdgeCount() > 1)
                {
                    IList<KeyValuePair<uint, bool>> eId2Dir = BRep.GetEdgesForRemoveEdge(id);
                    System.Diagnostics.Debug.Assert(eId2Dir.Count != 0);
                    {
                        int iE = 0;
                        for (; iE < eId2Dir.Count; iE++)
                        {
                            int jE = 0;
                            for (; jE < eId2Dir.Count; jE++)
                            {
                                if (iE == jE)
                                {
                                    continue;
                                }
                                if (eId2Dir[iE].Key == eId2Dir[jE].Key)
                                {
                                    System.Diagnostics.Debug.Assert(eId2Dir[iE].Value != eId2Dir[jE].Value);
                                    break;
                                }
                            }
                            if (jE == eId2Dir.Count)
                            {
                                break;
                            }
                        }
                        isDelCP = (iE == eId2Dir.Count);
                    }
                    if (!isDelCP)
                    {
                        double area = 0.0;
                        for (int iE = 0; iE < eId2Dir.Count; iE++)
                        {
                            uint eId = eId2Dir[iE].Key;
                            System.Diagnostics.Debug.Assert(IsElementId(CadElementType.Edge, eId));
                            Edge2D e = GetEdge(eId);
                            double earea = e.EdgeArea() + 
                                CadUtils.TriArea(e.GetVertex(true), e.GetVertex(false), new OpenTK.Vector2d(0, 0));
                            if (eId2Dir[iE].Value)
                            {
                                area += earea;
                            }
                            else
                            {
                                area -= earea;
                            }
                        }
                        if (area < 0)
                        {
                            isDelCP = true;
                        }
                    }
                }
                if (!BRep.RemoveEdge(id, isDelCP))
                {
                    System.Diagnostics.Debug.WriteLine( "Remove Edge B-Rep unsuccessfull : " + id);
                    return false;
                }
                EdgeArray.DeleteObject(id);
                if (!BRep.IsElementId(CadElementType.Loop, lLId))
                {
                    LoopArray.DeleteObject(lLId);
                }
                if (!BRep.IsElementId(CadElementType.Loop, rLId))
                {
                    LoopArray.DeleteObject(rLId);
                }
                System.Diagnostics.Debug.Assert(AssertValid() == 0);
                return true;
            }
            else if (type == CadElementType.Vertex)
            {
                VertexEdgeItr vItr = BRep.GetVertexEdgeItr(id);
                if (vItr.GetEdgeCount() == 2)
                {
                    uint eId1;
                    uint eId2;
                    bool isSameDir1;
                    bool isSameDir2;
                    vItr.GetAheadEdgeId(out eId1, out isSameDir1);
                    vItr.GetBehindEdgeId(out eId2, out isSameDir2);
                    Edge2D tmpEdge = GetEdge(eId1);
                    {
                        uint vId2 = BRep.GetEdgeVertexId(eId2, !isSameDir2);
                        System.Diagnostics.Debug.Assert(BRep.GetEdgeVertexId(eId1, isSameDir1) == id);
                        System.Diagnostics.Debug.Assert(BRep.GetEdgeVertexId(eId2, isSameDir2) == id);
                        tmpEdge.ConnectEdge(GetEdge(eId2), !isSameDir1, isSameDir1 != isSameDir2);
                        if (isSameDir1)
                        {
                            tmpEdge.SetVertexIds(vId2, tmpEdge.GetVertexId(false));
                            tmpEdge.SetVertexCoords(GetVertexCoord(vId2), tmpEdge.GetVertex(false));
                        }
                        else
                        {
                            tmpEdge.SetVertexIds(tmpEdge.GetVertexId(true), vId2);
                            tmpEdge.SetVertexCoords(tmpEdge.GetVertex(true), GetVertexCoord(vId2));
                        }
                    }
                    {
                        uint iPt0 = tmpEdge.GetVertexId(true);
                        uint iPt1 = tmpEdge.GetVertexId(false);
                        BoundingBox2D iBB = tmpEdge.GetBoundingBox();
                        IList<uint> eIds = BRep.GetElementIds(CadElementType.Edge);
                        for (int ijE = 0; ijE < eIds.Count; ijE++)
                        {
                            uint jEId = eIds[ijE];
                            if (jEId == eId2 || jEId == eId1)
                            {
                                continue;
                            }
                            Edge2D jEdge = GetEdge(jEId);
                            uint jPt0 = jEdge.GetVertexId(true);
                            uint jPt1 = jEdge.GetVertexId(false);
                            if ((iPt0 - jPt0) * (iPt0 - jPt1) * (iPt1 - jPt0) * (iPt1 - jPt1) != 0)
                            {
                                BoundingBox2D jBB = jEdge.GetBoundingBox();
                                if (!iBB.IsIntersect(jBB, MinClearance))
                                {
                                    continue;
                                }
                                double dist = tmpEdge.Distance(jEdge);
                                if (dist > MinClearance)
                                {
                                    continue;
                                }
                                return true;
                            }
                            else if (iPt0 == jPt0 && iPt1 == jPt1)
                            {
                                if (tmpEdge.IsCrossEdgeShareBothPoints(jEdge, true))
                                {
                                    return false;
                                }
                            }
                            else if (iPt0 == jPt1 && iPt1 == jPt0)
                            {
                                if (tmpEdge.IsCrossEdgeShareBothPoints(jEdge, false))
                                {
                                    return false;
                                }
                            }
                            else if (iPt0 == jPt0)
                            {
                                if (tmpEdge.IsCrossEdgeShareOnePoint(jEdge, true, true))
                                {
                                    return false;
                                }
                            }
                            else if (iPt0 == jPt1)
                            {
                                if (tmpEdge.IsCrossEdgeShareOnePoint(jEdge, true, false))
                                {
                                    return false;
                                }
                            }
                            else if (iPt1 == jPt0)
                            {
                                if (tmpEdge.IsCrossEdgeShareOnePoint(jEdge, false, true))
                                {
                                    return false;
                                }
                            }
                            else if (iPt1 == jPt1)
                            {
                                if (tmpEdge.IsCrossEdgeShareOnePoint(jEdge, false, false))
                                {
                                    return false;
                                }
                            }
                        }
                    }
                    if (!BRep.RemoveVertex(id))
                    {
                        return false;
                    }
                    System.Diagnostics.Debug.Assert(BRep.IsElementId(CadElementType.Edge, eId1));
                    System.Diagnostics.Debug.Assert(!BRep.IsElementId(CadElementType.Edge, eId2));
                    EdgeArray.DeleteObject(eId2);
                    Edge2D e1 = GetEdge(eId1);
                    e1.Copy(tmpEdge);
                    VertexArray.DeleteObject(id);
                    System.Diagnostics.Debug.Assert(AssertValid() == 0);
                    return true;
                }
                else if (vItr.GetEdgeCount() == 0)
                {
                    if (!BRep.RemoveVertex(id)) { return false; }
                    VertexArray.DeleteObject(id);
                    System.Diagnostics.Debug.Assert(AssertValid() == 0);
                    return true;
                }
            }
            return false;
        }

        protected int CheckLoop(uint lId)
        {
            {
                if (CheckLoopIntersection(lId))
                {
                    return 1;
                }
            }
            {
                for (LoopEdgeItr lItr = BRep.GetLoopEdgeItr(lId); !lItr.IsChildEnd; lItr.ShiftChildLoop())
                {
                    if (lItr.IsParent())
                    {
                        if (lItr.GetUseLoopType() != 2)
                        {
                            return 2;
                        }
                        if (GetLoopEdgeItrArea(lItr) < 0)
                        {
                            return 2;
                        }
                    }
                    else if (lItr.GetUseLoopType() == 2)
                    {
                        if (GetLoopEdgeItrArea(lItr) > 0)
                        {
                            return 2;
                        }
                    }
                }
            }
            {
                LoopEdgeItr pLItr = BRep.GetLoopEdgeItr(lId);
                for (LoopEdgeItr cLItr = BRep.GetLoopEdgeItr(lId); !cLItr.IsChildEnd; cLItr.ShiftChildLoop())
                {
                    if (cLItr.IsParent())
                    {
                        continue;
                    }
                    if (CheckLoopEdgeItrPointInOutLoopEdgeItr(cLItr, pLItr) != 0)
                    {
                        return 3;
                    }
                }
            }
            {
                for (LoopEdgeItr lItr1 = BRep.GetLoopEdgeItr(lId); !lItr1.IsChildEnd; lItr1.ShiftChildLoop())
                {
                    if (lItr1.IsParent())
                    {
                        continue;
                    }
                    for (LoopEdgeItr lItr2 = BRep.GetLoopEdgeItr(lId); !lItr2.IsChildEnd; lItr2.ShiftChildLoop())
                    {
                        if (lItr2.IsParent())
                        {
                            continue;
                        }
                        if (lItr1.IsSameUseLoop(lItr2))
                        {
                            continue;
                        }
                        if (CheckLoopEdgeItrPointInOutLoopEdgeItr(lItr1, lItr2) != 2)
                        {
                            return 4;
                        }
                    }
                }
            }
            {
                OpenTK.Vector2d zeroVec = new OpenTK.Vector2d(0, 0);
                for (LoopEdgeItr lItr = BRep.GetLoopEdgeItr(lId); !lItr.IsChildEnd; lItr.ShiftChildLoop())
                {
                    for (lItr.Begin(); !lItr.IsEnd(); lItr++)
                    {
                        uint vmId;
                        uint vfId;
                        OpenTK.Vector2d dir;
                        {
                            vmId = lItr.GetVertexId();
                            vfId = lItr.GetAheadVertexId();
                            uint eId0;
                            bool isSameDir0;
                            lItr.GetEdgeId(out eId0, out isSameDir0);
                            if (!EdgeArray.IsObjectId(eId0))
                            {
                                continue;
                            }
                            Edge2D e0 = GetEdge(eId0);
                            dir = e0.GetTangentEdge(isSameDir0);
                            System.Diagnostics.Debug.Assert(e0.GetVertexId(isSameDir0) == vmId);
                            System.Diagnostics.Debug.Assert(e0.GetVertexId(!isSameDir0) == vfId);
                        }
                        for (VertexEdgeItr vItr = BRep.GetVertexEdgeItr(vmId); !vItr.IsEnd(); vItr++)
                        {   // 点周りの辺をめぐる
                            uint eId0;
                            bool isSameDir0;
                            vItr.GetBehindEdgeId(out eId0, out isSameDir0);
                            if (!EdgeArray.IsObjectId(eId0))
                            {
                                continue;
                            }
                            Edge2D e0 = GetEdge(eId0);
                            System.Diagnostics.Debug.Assert(e0.GetVertexId(isSameDir0) == vmId);
                            if (e0.GetVertexId(!isSameDir0) == vfId)
                            {
                                continue;
                            }
                            OpenTK.Vector2d tan0 = e0.GetTangentEdge(isSameDir0);
                            uint eId1;
                            bool isSameDir1;
                            vItr.GetAheadEdgeId(out eId1, out isSameDir1);
                            if (!EdgeArray.IsObjectId(eId1))
                            {
                                continue;
                            }
                            Edge2D e1 = GetEdge(eId1);
                            System.Diagnostics.Debug.Assert(e1.GetVertexId(isSameDir1) == vmId);
                            if (e1.GetVertexId(!isSameDir1) == vfId)
                            {
                                continue;
                            }
                            OpenTK.Vector2d tan1 = e1.GetTangentEdge(isSameDir1);
                            double area0 = CadUtils.TriArea(tan1, zeroVec, tan0);
                            double area1 = CadUtils.TriArea(tan1, zeroVec, dir);
                            double area2 = CadUtils.TriArea(dir, zeroVec, tan0);
                            if ((area0 > 0.0 && area1 > 0.0 && area2 > 0.0) ||
                                   (area0 < 0.0 && (area1 > 0.0 || area2 > 0.0)))
                            {
                                return 5;
                            }
                        }
                    }
                }
            }
            return 0;
        }

        private bool CheckLoopIntersection(uint lId)
        {
            IList<uint> eIds = new List<uint>();
            if (BRep.IsElementId(CadElementType.Loop, lId))
            {
                for (LoopEdgeItr lItr = BRep.GetLoopEdgeItr(lId); !lItr.IsChildEnd; lItr.ShiftChildLoop())
                {
                    for (lItr.Begin(); !lItr.IsEnd(); lItr++)
                    {
                        uint eId;
                        bool isSameDir;
                        if (!lItr.GetEdgeId(out eId, out isSameDir))
                        {
                            continue;
                        }
                        if (lItr.IsEdgeBothSideSameLoop() && !isSameDir)
                        {
                            continue;
                        }
                        eIds.Add(eId);
                    }
                }
            }
            else
            {
                eIds = GetElementIds(CadElementType.Edge);
            }

            int eIdCnt = eIds.Count;
            for (int iE = 0; iE < eIdCnt; iE++)
            {
                Edge2D edge = GetEdge(eIds[iE]);
                if (edge.IsCrossEdgeSelf())
                {
                    return true;
                }
                uint iPt0 = edge.GetVertexId(true);
                uint iPt1 = edge.GetVertexId(false);
                BoundingBox2D iBB = edge.GetBoundingBox();
                for (int jE = iE + 1; jE < eIdCnt; jE++)
                {
                    Edge2D jEdge = GetEdge(eIds[jE]);
                    uint jPt0 = jEdge.GetVertexId(true);
                    uint jPt1 = jEdge.GetVertexId(false);
                    if ((iPt0 - jPt0) * (iPt0 - jPt1) * (iPt1 - jPt0) * (iPt1 - jPt1) != 0)
                    {
                        BoundingBox2D jBB = jEdge.GetBoundingBox();
                        if (!iBB.IsIntersect(jBB, MinClearance))
                        {
                            continue;
                        }
                        double dist = edge.Distance(jEdge);
                        if (dist < MinClearance)
                        {
                            return true;
                        }
                        continue;
                    }
                    else if (iPt0 == jPt0 && iPt1 == jPt1)
                    {
                        if (edge.IsCrossEdgeShareBothPoints(jEdge, true))
                        {
                            return true;
                        }
                    }
                    else if (iPt0 == jPt1 && iPt1 == jPt0)
                    {
                        if (edge.IsCrossEdgeShareBothPoints(jEdge, false))
                        {
                            return true;
                        }
                    }
                    else if (iPt0 == jPt0)
                    {
                        if (edge.IsCrossEdgeShareOnePoint(jEdge, true, true))
                        {
                            return true;
                        }
                    }
                    else if (iPt0 == jPt1)
                    {
                        if (edge.IsCrossEdgeShareOnePoint(jEdge, true, false))
                        {
                            return true;
                        }
                    }
                    else if (iPt1 == jPt0)
                    {
                        if (edge.IsCrossEdgeShareOnePoint(jEdge, false, true))
                        {
                            return true;
                        }
                    }
                    else if (iPt1 == jPt1)
                    {
                        if (edge.IsCrossEdgeShareOnePoint(jEdge, false, false))
                        {
                            return true;
                        }
                    }
                }
            }
            return false;
        }

        private bool CheckEdgeAgainstLoopIntersection(Edge2D edge, uint lId)
        {
            IList<uint> edgeIds = new List<uint>();
            if (BRep.IsElementId(CadElementType.Loop, lId))
            {
                for (LoopEdgeItr lItr = BRep.GetLoopEdgeItr(lId); !lItr.IsChildEnd; lItr.ShiftChildLoop())
                {
                    for (lItr.Begin(); !lItr.IsEnd(); lItr++)
                    {
                        uint eId;
                        bool isSameDir;
                        if (!lItr.GetEdgeId(out eId, out isSameDir))
                        {
                            continue;
                        }
                        if (lItr.IsEdgeBothSideSameLoop() && !isSameDir)
                        {
                            continue;
                        }
                        edgeIds.Add(eId);
                    }
                }
            }
            else
            {
                edgeIds = GetElementIds(CadElementType.Edge);
            }

            uint edgeIdCnt = (uint)edgeIds.Count;
            if (edge.IsCrossEdgeSelf())
            {
                return true;
            }
            uint iPt0 = edge.GetVertexId(true);
            uint iPt1 = edge.GetVertexId(false);
            BoundingBox2D iBB = edge.GetBoundingBox();
            for (int jE = 0; jE < edgeIdCnt; jE++)
            {
                Edge2D jEdge = GetEdge(edgeIds[jE]);
                uint jPt0 = jEdge.GetVertexId(true);
                uint jPt1 = jEdge.GetVertexId(false);
                if ((iPt0 - jPt0) * (iPt0 - jPt1) * (iPt1 - jPt0) * (iPt1 - jPt1) != 0)
                {
                    BoundingBox2D jBB = jEdge.GetBoundingBox();
                    if (!iBB.IsIntersect(jBB, MinClearance))
                    {
                        continue;
                    }
                    double dist = edge.Distance(jEdge);
                    if (dist < MinClearance)
                    {
                        return true;
                    }
                    continue;
                }
                else if (iPt0 == jPt0 && iPt1 == jPt1)
                {
                    if (edge.IsCrossEdgeShareBothPoints(jEdge, true))
                    {
                        return true;
                    }
                }
                else if (iPt0 == jPt1 && iPt1 == jPt0)
                {
                    if (edge.IsCrossEdgeShareBothPoints(jEdge, false))
                    {
                        return true;
                    }
                }
                else if (iPt0 == jPt0)
                {
                    if (edge.IsCrossEdgeShareOnePoint(jEdge, true, true))
                    {
                        return true;
                    }
                }
                else if (iPt0 == jPt1)
                {
                    if (edge.IsCrossEdgeShareOnePoint(jEdge, true, false))
                    {
                        return true;
                    }
                }
                else if (iPt1 == jPt0)
                {
                    if (edge.IsCrossEdgeShareOnePoint(jEdge, false, true))
                    {
                        return true;
                    }
                }
                else if (iPt1 == jPt1)
                {
                    if (edge.IsCrossEdgeShareOnePoint(jEdge, false, false))
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        public LoopEdgeItr GetLoopEdgeItr(uint lId)
        {
            return new LoopEdgeItr(BRep, lId);
        }

        public VertexEdgeItr GetVertexEdgeItr(uint vId)
        {
            return new VertexEdgeItr(BRep, vId);
        }

        public bool Serialize(Serializer arch)
        {
            if (arch.IsLoading)
            {
                // load
                Clear();

                ////////////////
                string className;
                string[] values;
                
                className = arch.ReadDepthClassName();
                if (className != "CadObj2D")
                {
                    return false;
                }
                int nv;
                int ne;
                int nl;
                {
                    values = arch.GetValues();
                    nv = int.Parse(values[0]);
                    ne = int.Parse(values[1]);
                    nl = int.Parse(values[2]);
                    System.Diagnostics.Debug.Assert(nv > 0);
                    System.Diagnostics.Debug.Assert(ne > 0);
                    System.Diagnostics.Debug.Assert(nl > 0);
                    /*
                    m_VertexSet.Reserve(nv * 2);
                    m_EdgeSet.Reserve(ne * 2);
                    m_LoopSet.Reserve(nl * 2);
                    */
                }
                arch.ShiftDepth(true);
                ////////////////////////////////////////////////
                for (int iv = 0; iv < nv; iv++)
                {
                    className = arch.ReadDepthClassName();
                    System.Diagnostics.Debug.Assert(className == "CVertex2D");
                    int id;
                    values = arch.GetValues();
                    id = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(id > 0);

                    double x;
                    double y;
                    values = arch.GetValues();
                    x = double.Parse(values[0]);
                    y = double.Parse(values[1]);

                    // added by ryujimiya
                    double[] c = new double[3];
                    values = arch.GetValues();
                    for (int i = 0; i < 3; i++)
                    {
                        c[i] = double.Parse(values[i]);
                    }
                    // added by ryujimiya

                    Vertex2D v = new Vertex2D(new OpenTK.Vector2d(x, y));
                    c.CopyTo(v.Color, 0); // added by ryujimiya
                    uint tmpId = VertexArray.AddObject((uint)id, v);
                    System.Diagnostics.Debug.Assert(tmpId == id);
                }
                for (int ie = 0; ie < ne; ie++)
                {
                    className = arch.ReadDepthClassName();
                    System.Diagnostics.Debug.Assert(className == "CEdge2D");
                    int id;
                    values = arch.GetValues();
                    id = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(id > 0);

                    int sVId;
                    int eVId;
                    values = arch.GetValues();
                    sVId = int.Parse(values[0]);
                    eVId = int.Parse(values[1]);
                    System.Diagnostics.Debug.Assert(sVId > 0 && eVId > 0);

                    CurveType type = CurveType.CurveLine; // dummy
                    {
                        int itype;
                        values = arch.GetValues();
                        itype = int.Parse(values[0]);
                        if (itype == 0) { type = CurveType.CurveLine; }
                        else if (itype == 1) { type = CurveType.CurveArc; }
                        else if (itype == 2) { type = CurveType.CurvePolyline; }
                        else if (itype == 3) { type = CurveType.CurveBezier; }
                        else { System.Diagnostics.Debug.Assert(false); }
                    }

                    int iIsLeftSide;
                    double dist;
                    values = arch.GetValues();
                    iIsLeftSide = int.Parse(values[0]);
                    dist = double.Parse(values[1]);
                    System.Diagnostics.Debug.Assert(iIsLeftSide == 0 || iIsLeftSide == 1);
                    bool isLeftSide = (iIsLeftSide != 0);

                    IList<double> relCos = new List<double>();
                    int npo;
                    values = arch.GetValues();
                    npo = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(npo >= 0);
                    for (uint ipo = 0; ipo < npo; ipo++)
                    {
                        double x;
                        double y;
                        values = arch.GetValues();
                        x = double.Parse(values[0]);
                        y = double.Parse(values[1]);
                        relCos.Add(x);
                        relCos.Add(y);
                    }

                    // added by ryujimiya
                    double[] c = new double[3];
                    values = arch.GetValues();
                    for (int i = 0; i < 3; i++)
                    {
                        c[i] = double.Parse(values[i]);
                    }
                    // added by ryujimiya

                    Edge2D e = new Edge2D((uint)sVId, (uint)eVId);
                    if (type == CurveType.CurveLine)
                    {
                        e.SetCurveLine();
                    }
                    else if (type == CurveType.CurveArc)
                    {
                        e.SetCurveArc(isLeftSide, dist);
                    }
                    else if (type == CurveType.CurvePolyline)
                    {
                        e.SetCurvePolyline(relCos);
                    }
                    else if (type == CurveType.CurveBezier)
                    {
                        // TODO:
                    }
                    c.CopyTo(e.Color, 0); // added by ryujimiya

                    uint tmpId = EdgeArray.AddObject((uint)id, e);
                    System.Diagnostics.Debug.Assert(tmpId == id);
                }

                for (int il = 0; il < nl; il++)
                {
                    className = arch.ReadDepthClassName();
                    System.Diagnostics.Debug.Assert(className == "CLoop2D");
                    int id;
                    values = arch.GetValues();
                    id = int.Parse(values[0]);
                    System.Diagnostics.Debug.Assert(id > 0);

                    int layer;
                    values = arch.GetValues();
                    layer = int.Parse(values[0]);

                    double[] c = new double[3];
                    values = arch.GetValues();
                    for (int i = 0; i < 3; i++)
                    {
                        c[i] = double.Parse(values[i]);
                    }
                    Loop2D l = new Loop2D();
                    l.Layer = (uint)layer;
                    c.CopyTo(l.Color, 0);

                    uint tmpId = LoopArray.AddObject((uint)id, l);
                    System.Diagnostics.Debug.Assert(tmpId == id);
                }

                BRep.Serialize(arch);

                AssertValid();

                arch.ShiftDepth(false);

                return true;
            }
            else
            {
                // write
                string line;

                // class name and size
                arch.WriteDepthClassName("CadObj2D");
                line = string.Format("{0} {1} {2}",
                    VertexArray.GetObjectIds().Count,
                    EdgeArray.GetObjectIds().Count,
                    LoopArray.GetObjectIds().Count);
                arch.WriteLine(line);
                arch.ShiftDepth(true);

                {
                    // print Vertex2D
                    IList<uint> ids = VertexArray.GetObjectIds();
                    foreach (uint vId in ids)
                    {
                        System.Diagnostics.Debug.Assert(VertexArray.IsObjectId(vId));
                        Vertex2D v = VertexArray.GetObject(vId);

                        arch.WriteDepthClassName("CVertex2D");

                        line = string.Format("{0}", vId);
                        arch.WriteLine(line);

                        line = string.Format("{0} {1}", v.Point.X, v.Point.Y);
                        arch.WriteLine(line);

                        // added by ryujimiya
                        line = string.Format("{0} {1} {2}", v.Color[0], v.Color[1], v.Color[2]);
                        arch.WriteLine(line);
                    }
                }
                {
                    // print Edge2D
                    IList<uint> ids = EdgeArray.GetObjectIds();
                    foreach (uint eId in ids)
                    {
                        System.Diagnostics.Debug.Assert(EdgeArray.IsObjectId(eId));
                        Edge2D e = GetEdge(eId);

                        arch.WriteDepthClassName("CEdge2D");

                        line = string.Format("{0}", eId);
                        arch.WriteLine(line);

                        line = string.Format("{0} {1}", e.GetVertexId(true), e.GetVertexId(false));
                        arch.WriteLine(line);

                        {
                            CurveType type = e.CurveType;
                            int itype = 0;
                            if (type == CurveType.CurveLine) { itype = 0; }
                            else if (type == CurveType.CurveArc) { itype = 1; }
                            else if (type == CurveType.CurvePolyline) { itype = 2; }
                            else if (type == CurveType.CurveBezier) { itype = 3; }
                            else { System.Diagnostics.Debug.Assert(false); }
                            line = string.Format("{0}", itype);
                            arch.WriteLine(line);
                        }
                        {
                            bool isLeftSide;
                            double dist;
                            e.GetCurveArc(out isLeftSide, out dist);
                            int iIsLeftSide = (isLeftSide) ? 1 : 0;
                            line = string.Format("{0} {1}", iIsLeftSide, dist);
                            arch.WriteLine(line);
                        }
                        IList<double> relCos = e.GetCurveRelPoint();
                        uint n = (uint)(relCos.Count / 2);
                        line = string.Format("{0}", n);
                        arch.WriteLine(line);
                        for (int i = 0; i < n; i++)
                        {
                            line = string.Format("{0} {1}", relCos[i * 2 + 0], relCos[i * 2 + 1]);
                            arch.WriteLine(line);
                        }
                        // added by ryujimiya
                        {
                            line = string.Format("{0} {1} {2}", e.Color[0], e.Color[1], e.Color[2]);
                            arch.WriteLine(line);
                        }
                        // endof added
                    }
                }
                { 
                    // print Loop2D
                    IList<uint> ids = LoopArray.GetObjectIds();
                    foreach (uint lId in ids)
                    {
                        System.Diagnostics.Debug.Assert(LoopArray.IsObjectId(lId));
                        Loop2D l = LoopArray.GetObject(lId);

                        arch.WriteDepthClassName("CLoop2D");

                        line = string.Format("{0}", lId);
                        arch.WriteLine(line);

                        line = string.Format("{0}", l.Layer);
                        arch.WriteLine(line);

                        line = string.Format("{0} {1} {2}", l.Color[0], l.Color[1], l.Color[2]);
                        arch.WriteLine(line);
                    }
                }

                BRep.Serialize(arch);

                arch.ShiftDepth(false);
            }
            return true;
        }

        public void Copy(CadObject2D src)
        {
            Clear();

            // Vertex
            {
                IList<uint> srcIds = src.VertexArray.GetObjectIds();
                foreach (uint srcId in srcIds)
                {
                    Vertex2D srcV = src.VertexArray.GetObject(srcId);
                    uint id = srcId;
                    Vertex2D v = new Vertex2D(srcV);
                    uint tmpId = VertexArray.AddObject(id, v);
                    System.Diagnostics.Debug.Assert(tmpId == id);
                }
            }
            // Edge
            {
                IList<uint> srcIds = src.EdgeArray.GetObjectIds();
                foreach (uint srcId in srcIds)
                {
                    Edge2D srcE = src.EdgeArray.GetObject(srcId);
                    uint id = srcId;
                    Edge2D e = new Edge2D(srcE);
                    uint tmpId = EdgeArray.AddObject(id, e);
                    System.Diagnostics.Debug.Assert(tmpId == id);
                }
            }
            // Loop
            {
                IList<uint> srcIds = src.LoopArray.GetObjectIds();
                foreach (uint srcId in srcIds)
                {
                    Loop2D srcL = src.LoopArray.GetObject(srcId);
                    uint id = srcId;
                    Loop2D e = new Loop2D(srcL);
                    uint tmpId = LoopArray.AddObject(id, e);
                    System.Diagnostics.Debug.Assert(tmpId == id);
                }
            }

            BRep.Copy(src.BRep);

            AssertValid();
        }

    }
}
