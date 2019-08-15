using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class CadObject2DMove : CadObject2D
    {
        private CadEdge2DPolyline Polyline = new CadEdge2DPolyline();

        public bool MoveEdgeCtrl(uint eId, uint iCtrl, OpenTK.Vector2d deltaVec)
        {
            if (!IsElementId(CadElementType.Edge, eId))
            {
                return false;
            }

            IList<double> oldRelCos = null;
            {
                Edge2D e = GetEdge(eId);
                oldRelCos = e.GetCurveRelPoint();
            }

            Edge2D edge = GetEdge(eId);
            if (iCtrl >= edge.GetCurvePoint().Count)
            {
                return true;
            }
            {
                OpenTK.Vector2d sPt = edge.GetVertexCoord(true);
                OpenTK.Vector2d ePt = edge.GetVertexCoord(false);
                double sqLen = CadUtils.SquareLength(ePt - sPt);
                var eh = (ePt - sPt) * (1 / sqLen);
                OpenTK.Vector2d ev = new OpenTK.Vector2d(-eh.Y, eh.X);
                var p0 = edge.GetCurvePoint()[(int)iCtrl] + deltaVec;
                IList<double> aRelCo = edge.GetCurveRelPoint();
                aRelCo[(int)(iCtrl * 2 + 0)] = OpenTK.Vector2d.Dot(p0 - sPt, eh);
                aRelCo[(int)(iCtrl * 2 + 1)] = OpenTK.Vector2d.Dot(p0 - sPt, ev);
                edge.SetCurveRelPoint(aRelCo);
            }

            HashSet<uint> setLId = new HashSet<uint>();
            for (VertexEdgeItr vItr = this.BRep.GetVertexEdgeItr(edge.GetVertexId(true)); !vItr.IsEnd(); vItr.Next())
            {
                uint tmpLId = vItr.GetLoopId();
                if (!setLId.Contains(tmpLId))
                {
                    setLId.Add(tmpLId);
                }
            }
            for (VertexEdgeItr vItr = this.BRep.GetVertexEdgeItr(edge.GetVertexId(false)); !vItr.IsEnd(); vItr.Next())
            {
                uint tmpLId = vItr.GetLoopId();
                if (!setLId.Contains(tmpLId))
                {
                    setLId.Add(tmpLId);
                }
            }

            foreach (uint lId in setLId)
            {
                if (!this.BRep.IsElementId(CadElementType.Loop, lId))
                {
                    continue;
                }
                if (CheckLoop(lId) != 0)
                {
                    // fail
                    edge.SetCurveRelPoint(oldRelCos);
                    return true;
                }
            }
            return true;
        }

        public bool MoveLoop(uint lId, OpenTK.Vector2d vecDelta)
        {
            if (!IsElementId(CadElementType.Loop, lId))
            {
                return false;
            }
            Dictionary<uint, OpenTK.Vector2d> mapOldVec = new Dictionary<uint, OpenTK.Vector2d>();
            HashSet<uint> setLId = new HashSet<uint>();  // check these loop for intersection detection
            for (LoopEdgeItr lItr = this.BRep.GetLoopEdgeItr(lId); !lItr.IsChildEnd; lItr.ShiftChildLoop())
            {
                for (lItr.Begin(); !lItr.IsEnd(); lItr.Next())
                {
                    uint vId = lItr.GetVertexId();
                    if (mapOldVec.ContainsKey(vId))
                    {
                        continue; // this point is already moved
                    }
                    mapOldVec.Add(vId, GetVertexCoord(vId));
                    Vertex2D vertex = VertexArray.GetObject(vId);
                    {
                        double oldX = vertex.Point.X;
                        double oldY = vertex.Point.Y;
                        vertex.Point = new OpenTK.Vector2d(oldX + vecDelta.X, oldY + vecDelta.Y);
                    }
                    uint eId0;
                    bool isSameDir0;
                    lItr.GetEdgeId(out eId0, out isSameDir0);
                    if (!IsElementId(CadElementType.Edge, eId0))
                    {
                        continue;  // this is point
                    }
                    uint lLId;
                    uint rLId;
                    this.GetEdgeLoopId(out lLId, out rLId, eId0);
                    if (!setLId.Contains(lLId))
                    {
                        setLId.Add(lLId);
                    }
                    if (!setLId.Contains(rLId))
                    {
                        setLId.Add(rLId);
                    }
                }
            }
            foreach (uint lIdl0 in setLId)
            {
                if (!IsElementId(CadElementType.Loop, lIdl0))
                {
                    continue;  // id0 can be 0 if this loop(id_l) have open boundary
                }
                if (CheckLoop(lIdl0) != 0)
                {
                    // fail
                    foreach (var oldVec in mapOldVec)
                    {
                        uint vId = oldVec.Key;
                        Vertex2D vertex = VertexArray.GetObject(vId);
                        vertex.Point = oldVec.Value;
                    }
                    return false;
                }
            }
            return true;
        }

        public bool MoveEdge(uint eId, OpenTK.Vector2d vecDelta)
        {
            System.Diagnostics.Debug.Assert(IsElementId(CadElementType.Edge, eId));
            uint sVId = GetEdgeVertexId(eId, true);
            uint eVId = GetEdgeVertexId(eId, false);

            OpenTK.Vector2d preSVec;
            {
                // 点を動かす→駄目だったら元に戻す
                Vertex2D vertex = VertexArray.GetObject(sVId);
                preSVec = vertex.Point;
                {
                    double oldX = vertex.Point.X;
                    double oldY = vertex.Point.Y;
                    vertex.Point = new OpenTK.Vector2d(oldX + vecDelta.X, oldY + vecDelta.Y);
                }
            }

            OpenTK.Vector2d preEVec;
            {   // 点を動かす→駄目だったら元に戻す
                Vertex2D ver = VertexArray.GetObject(eVId);
                preEVec = ver.Point;
                {
                    double oldX = ver.Point.X;
                    double oldY = ver.Point.Y;
                    ver.Point = new OpenTK.Vector2d(oldX + vecDelta.X, oldY + vecDelta.Y);
                }
            }

            try
            {
                // Check Interfarance
                HashSet<uint> lIds = new HashSet<uint>();
                for (VertexEdgeItr vItr = this.BRep.GetVertexEdgeItr(sVId); !vItr.IsEnd(); vItr.Next())
                {
                    uint lId = vItr.GetLoopId();
                    if (IsElementId(CadElementType.Loop, lId))
                    {
                        if (!lIds.Contains(lId))
                        {
                            lIds.Add(lId);
                        }
                        else
                        {
                            continue;
                        }
                        if (CheckLoop(lId) != 0)
                        {
                            throw new InvalidOperationException();
                        }
                    }
                }
                for (VertexEdgeItr vItr = this.BRep.GetVertexEdgeItr(eVId); !vItr.IsEnd(); vItr.Next())
                {
                    uint lId = vItr.GetLoopId();
                    if (IsElementId(CadElementType.Loop, lId))
                    {
                        if (!lIds.Contains(lId))
                        {
                            lIds.Add(lId);
                        }
                        else
                        {
                            continue;
                        }
                        if (CheckLoop(lId) != 0)
                        {
                            throw new InvalidOperationException();
                        }
                    }
                }
            }
            catch(Exception exception)
            {
                ////////////////////////////////  
                //if the operation fails
                {   // 動かした点を元に戻す
                    Vertex2D vertex = VertexArray.GetObject(sVId);
                    vertex.Point = preSVec;
                }
                {   // 動かした点を元に戻す
                    Vertex2D vertex = VertexArray.GetObject(eVId);
                    vertex.Point = preEVec;
                }
                return false;
            }

            return true;
        }

        public bool MoveVertex(uint vId, OpenTK.Vector2d vec)
        {
            if (!VertexArray.IsObjectId(vId))
            {
                return false;
            }

            OpenTK.Vector2d preVec;
            {
                // store point to move point back in case it fails
                Vertex2D vertex = VertexArray.GetObject(vId);
                preVec = vertex.Point;
            }

            OpenTK.Vector2d dist = vec;
            {
                // move point
                Vertex2D vertex = VertexArray.GetObject(vId);
                vertex.Point = dist;
            }

            try
            {
                VertexEdgeItr vItr = this.BRep.GetVertexEdgeItr(vId);
                if (vItr.GetEdgeCount() == 0)
                {
                    // move point inside loop
                    uint lId = vItr.GetLoopId();
                    if (IsElementId(CadElementType.Loop, lId))
                    {
                        // ignore vtx(id_v) in the signd distance computation
                        double distance = SignedDistancePointLoop(lId, vec, vId);
                        if (distance < this.MinClearance)
                        {
                            throw new InvalidOperationException();
                        }
                    }
                }
                else
                {   
                    // move point adjacent to loop
                    HashSet<uint> lIds = new HashSet<uint>();
                    for (; !vItr.IsEnd(); vItr.Next())
                    {
                        uint lId = vItr.GetLoopId();
                        if (IsElementId(CadElementType.Loop, lId))
                        {
                            if (!lIds.Contains(lId))
                            {
                                lIds.Add(lId);
                            }
                            else
                            {
                                continue;
                            }
                            if (CheckLoop(lId) != 0)
                            {
                                throw new InvalidOperationException();
                            }
                        }
                    }
                }
            }
            catch(Exception exception)
            {
                ////////////////////////////////
                // fail
                {   // reset the moved point
                    Vertex2D vertex = VertexArray.GetObject(vId);
                    vertex.Point = preVec;
                }
                return false;
            }
            return true;
        }

        public bool MoveVertex(IList<KeyValuePair<uint, OpenTK.Vector2d>> idVecs)
        {
            IList<OpenTK.Vector2d> oldVecs = new List<OpenTK.Vector2d>();
            try
            {
                foreach (var pair in idVecs)
                {
                    uint id_v = pair.Key;
                    if (!VertexArray.IsObjectId(id_v))
                    {
                        throw new InvalidOperationException();
                    }
                    Vertex2D ver = VertexArray.GetObject(id_v);
                    oldVecs.Add(ver.Point);
                    ver.Point = pair.Value;
                }
                {
                    IList<uint> loopIds = GetElementIds(CadElementType.Loop);
                    foreach (uint lId in loopIds)
                    {
                        if (CheckLoop(lId) != 0)
                        {
                            throw new InvalidOperationException();
                        }
                    }
                }
            }
            catch (Exception exception)
            {
                ////////////////////////////////
                // fail
                // 動かした点を元に戻す
                for (int i = 0; i < oldVecs.Count; i++)
                {
                    uint id_v = idVecs[i].Key;
                    Vertex2D ver = VertexArray.GetObject(id_v);
                    ver.Point = oldVecs[i];
                }
                return false;
            }
            return true;
        }

        public bool DragArc(uint eId, OpenTK.Vector2d vec)
        {
            if (!IsElementId(CadElementType.Edge, eId))
            {
                return false;
            }
            if (GetEdgeCurveType(eId) != CurveType.CurveArc)
            {
                return true;
            }

            bool oldIsLeftSide;
            double oldDist;
            {
                Edge2D e = GetEdge(eId);
                e.GetCurveArc(out oldIsLeftSide, out oldDist);
            }

            Edge2D edge = GetEdge(eId);
            System.Diagnostics.Debug.Assert(edge.CurveType == CurveType.CurveArc);
            {
                OpenTK.Vector2d sPt = edge.GetVertexCoord(true);
                OpenTK.Vector2d ePt = edge.GetVertexCoord(false);
                double edgeLen = Math.Sqrt(CadUtils.SquareLength(sPt, ePt));
                if (Math.Abs(CadUtils.TriHeight(vec, sPt, ePt)) > edgeLen * 0.02)
                {
                    OpenTK.Vector2d cPt;
                    bool ret = CadUtils.CenterCircumcircle(sPt, ePt, vec, out cPt);
                    System.Diagnostics.Debug.Assert(ret);
                    double dist = CadUtils.TriHeight(cPt, sPt, ePt);
                    double distRatio = dist / edgeLen;
                    edge.SetCurveArc(CadUtils.TriArea(sPt, ePt, vec) > 0, distRatio);
                }
                else
                {
                    return true;
                }
            }

            HashSet<uint> setLId = new HashSet<uint>();
            for (VertexEdgeItr vItr = this.BRep.GetVertexEdgeItr(edge.GetVertexId(true)); !vItr.IsEnd(); vItr.Next())
            {
                uint tmpLId = vItr.GetLoopId();
                if (!setLId.Contains(tmpLId))
                {
                    setLId.Add(tmpLId);
                }
            }
            for (VertexEdgeItr vItr = this.BRep.GetVertexEdgeItr(edge.GetVertexId(false)); !vItr.IsEnd(); vItr.Next())
            {
                uint tmpLId = vItr.GetLoopId();
                if (!setLId.Contains(tmpLId))
                {
                    setLId.Add(tmpLId);
                }
            }

            foreach (uint lId in setLId)
            {
                if (!this.BRep.IsElementId(CadElementType.Loop, lId))
                {
                    continue;
                }
                if (CheckLoop(lId) != 0)
                {
                    // fail
                    edge.SetCurveArc(oldIsLeftSide, oldDist);
                    return true;
                }
            }
            return true;
        }

        public bool PreCompDragPolyline(uint eId, OpenTK.Vector2d pickPos)
        {
            if (!IsElementId(CadElementType.Edge, eId))
            {
                return false;
            }
            if (GetEdgeCurveType(eId) != CurveType.CurveArc) // 2が指定されている 
            {
                return true;
            }

            Polyline.SetCadEdge(this, eId, pickPos);

            return true;
        }

        public bool DragPolyline(uint eId, OpenTK.Vector2d dist)
        {
            if (!IsElementId(CadElementType.Edge, eId))
            {
                return false;
            }
            if (GetEdgeCurveType(eId) != CurveType.CurveArc) // 2が指定されている 
            {
                return true;
            }

            Polyline.Drag(this, dist);

            return true;
        }

        public bool SmoothingPolylineEdge(uint eId, uint nIter, OpenTK.Vector2d pos, double radius)
        {
            if (!this.IsElementId(CadElementType.Edge, eId))
            {
                return false;
            }
            if (GetEdgeCurveType(eId) != CurveType.CurvePolyline)
            {
                return true;
            }
            if (nIter == 0)
            {
                return true;
            }
            Edge2D edge = GetEdge(eId);
            IList<double> axys = edge.GetCurveRelPoint();
            System.Diagnostics.Debug.Assert(axys.Count % 2 == 0);
            uint nno = (uint)axys.Count / 2;
            IList<uint> aIndNo = new List<uint>();
            if (radius > 0)
            {
                uint[] tmp = new uint[nno];
                for (int i = 0; i < nno; i++)
                {
                    if (i < aIndNo.Count)
                    {
                        tmp[i] = aIndNo[i];
                    }
                    else
                    {
                        tmp[i] = 0;
                    }
                }
                aIndNo.Clear();
                aIndNo = tmp.ToList();

                OpenTK.Vector2d sPt = edge.GetVertexCoord(true);
                OpenTK.Vector2d ePt = edge.GetVertexCoord(false);
                var v0 = ePt - sPt;
                OpenTK.Vector2d v1 = new OpenTK.Vector2d(-v0.Y, v0.X);
                for (int ino = 0; ino < nno; ino++)
                {
                    var p = sPt + v0 * axys[ino * 2 + 0] + v1 * axys[ino * 2 + 1];
                    if (CadUtils.SquareLength(p - pos) < radius * radius)
                    {
                        aIndNo.Add((uint)ino);
                    }
                }
            }
            else
            {
                aIndNo = new List<uint>();
                for (int ino = 0; ino < nno; ino++)
                {
                    aIndNo[ino] = (uint)ino;
                }
            }

            // strength of smoothing(0から1)
            double w = 0.8;

            for (uint iiter = 0; iiter < nIter; iiter++)
            {
                for (int iIndNo = 0; iIndNo < aIndNo.Count; iIndNo++)
                {
                    int ino = (int)aIndNo[iIndNo];
                    double[] po0 = new double[2];
                    double[] po1 = new double[2];
                    double[] po2 = new double[2];
                    if (ino == 0)
                    {
                        po0[0] = 0;
                        po0[1] = 0;
                    }
                    else
                    {
                        po0[0] = axys[ino * 2 - 2];
                        po0[1] = axys[ino * 2 - 1];
                    }
                    po1[0] = axys[ino * 2 + 0];
                    po1[1] = axys[ino * 2 + 1];
                    if (ino == nno - 1)
                    {
                        po2[0] = 1; po2[1] = 0;
                    }
                    else
                    {
                        po2[0] = axys[ino * 2 + 2];
                        po2[1] = axys[ino * 2 + 3];
                    }
                    ////////////////
                    double[] mPt = { (po0[0] + po2[0]) * 0.5, (po0[1] + po2[1]) * 0.5 };
                    double[] dPt = { (1.0 - w) * po1[0] + w * mPt[0], (1.0 - w) * po1[1] + w * mPt[1] };
                    axys[ino * 2 + 0] = dPt[0];
                    axys[ino * 2 + 1] = dPt[1];
                }
            }
            edge.SetCurveRelPoint(axys);
            return true;
        }
    }
}
