using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class LiftEdgeProp
    {
        public uint OutVId { get; set; }
        public uint InVId { get; set; }
        public uint EId { get; set; }
        public bool IsSameDir { get; set; }
        public bool IsVert { get; set; }

        public LiftEdgeProp(uint vId, uint eId, bool isSameDir, bool isVert)
        {
            OutVId = vId;
            EId = eId;
            IsSameDir = isSameDir;
            IsVert = isVert;
        }
    }

    public class Cad3D
    {
        protected ObjectArray<Solid> SolidArray = new ObjectArray<Solid>();
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
                    Loop3D l = new Loop3D(srcL);
                    uint tmpId = LoopArray.AddObject(id, l);
                    System.Diagnostics.Debug.Assert(tmpId == id);
                }
            }
            // Solid
            {
                IList<uint> srcIds = src.SolidArray.GetObjectIds();
                foreach (uint srcId in srcIds)
                {
                    Solid srcS = src.SolidArray.GetObject(srcId);
                    uint id = srcId;
                    Solid s = new Solid(srcS);
                    uint tmpId = SolidArray.AddObject(id, s);
                    System.Diagnostics.Debug.Assert(tmpId == id);
                }
            }

            BRep.Copy(src.BRep);

            AssertValid();
        }

        public void Clear()
        {
            SolidArray.Clear();
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
            else if (type == CadElementType.Solid)
            {
                throw new NotImplementedException();
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }

            return res;
        }

        public ConnectVertexRes ConnectVertexLine(uint vId1, uint vId2)
        {
            ConnectVertexRes res = new ConnectVertexRes();
            res.VId1 = vId1;
            res.VId2 = vId2;

            IList<uint> chkLIds = new List<uint>();
            uint lId0 = 0;
            VertexEdgeItr vItr1 = BRep.GetVertexEdgeItr(vId1);
            for (vItr1.Begin(); !vItr1.IsEnd(); vItr1.Next())
            {
                uint lId1 = vItr1.GetLoopId();
                if (!chkLIds.Contains(lId1))
                {
                    chkLIds.Add(lId1);
                }
            }
            VertexEdgeItr vItr2 = BRep.GetVertexEdgeItr(vId2);
            for (vItr2.Begin(); !vItr2.IsEnd(); vItr2.Next())
            {
                uint lId2 = vItr2.GetLoopId();
                if (chkLIds.Contains(lId2))
                {
                    lId0 = lId2;
                    break;
                }
            }
            for (vItr1.Begin(); !vItr1.IsEnd(); vItr1.Next())
            {
                uint lId1 = vItr1.GetLoopId();
                if (lId1 == lId0)
                {
                    break;
                }
            }
            System.Diagnostics.Debug.Assert(vItr1.GetLoopId() == lId0);
            System.Diagnostics.Debug.Assert(vItr2.GetLoopId() == lId0);

            res = BRep.ConnectVertex(vItr1, vItr2, true);
            if (res.AddEId != 0)
            {
                EdgeArray.AddObject(res.AddEId, new Edge3D());
            }
            if (res.AddLId != lId0 && res.AddLId != 0)
            {
                OpenTK.Vector3d origin;
                OpenTK.Vector3d normal;
                OpenTK.Vector3d xdir;
                CalcLoopNormal(res.AddLId, out origin, out normal, out xdir);

                LoopArray.AddObject(res.AddLId, new Loop3D(origin, normal, xdir));
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
            {
                OpenTK.Vector3d origin;
                OpenTK.Vector3d normal;
                OpenTK.Vector3d xdir;
                CalcLoopNormal(res.AddLId, out origin, out normal, out xdir);

                LoopArray.AddObject(res.AddLId, new Loop3D(origin, normal, xdir));
            }
            System.Diagnostics.Debug.Assert(AssertValid() == 0);

            return res;
        }

        public AddSurfaceRes AddCube(IList<OpenTK.Vector3d> pts)
        {
            AddSurfaceRes res = new AddSurfaceRes();

            System.Diagnostics.Debug.Assert(pts.Count == 8);
            uint[] vIds = new uint[8];
            uint[] eIds = new uint[12];
            uint[] lIds = new uint[6];

            IList<OpenTK.Vector3d> pts1 = new List<OpenTK.Vector3d>();
            for (int i = 0; i < 4; i++)
            {
                pts1.Add(pts[i]);
            }
            for (int i = 0; i < pts1.Count; i++)
            {
                var pt = pts1[i];
                vIds[i] = AddVertex(CadElementType.Solid, 0, pt).AddVId;
            }

            IList<OpenTK.Vector3d> pts2 = new List<OpenTK.Vector3d>();
            for (int i = 0; i < 4; i++)
            {
                pts2.Add(pts[i + 4]);
            }
            for (int i = 0; i < pts2.Count; i++)
            {
                var pt = pts2[i];
                vIds[i + 4] = AddVertex(CadElementType.Solid, 0, pt).AddVId;
            }

            eIds[0] = ConnectVertexLine(vIds[0], vIds[1]).AddEId;
            eIds[1] = ConnectVertexLine(vIds[1], vIds[2]).AddEId;
            eIds[2] = ConnectVertexLine(vIds[2], vIds[3]).AddEId;
            var res1 = ConnectVertexLine(vIds[3], vIds[0]);
            eIds[3] = res1.AddEId;
            lIds[0] = res1.AddLId;

            eIds[4] = ConnectVertexLine(vIds[0], vIds[4]).AddEId;
            eIds[5] = ConnectVertexLine(vIds[1], vIds[5]).AddEId;
            eIds[6] = ConnectVertexLine(vIds[2], vIds[6]).AddEId;
            eIds[7] = ConnectVertexLine(vIds[3], vIds[7]).AddEId;

            var res2 = ConnectVertexLine(vIds[4], vIds[5]);
            eIds[8] = res2.AddEId;
            lIds[1] = res2.AddLId;
            var res3 = ConnectVertexLine(vIds[5], vIds[6]);
            eIds[9] = res3.AddEId;
            lIds[2] = res3.AddLId;
            var res4 = ConnectVertexLine(vIds[6], vIds[7]);
            eIds[10] = res4.AddEId;
            lIds[3] = res4.AddLId;
            var res5 = ConnectVertexLine(vIds[7], vIds[4]);
            eIds[11] = res5.AddEId;
            lIds[4] = res5.AddLId;
            lIds[5] = SealHole(eIds[11], false);

            res.AddVIds = vIds;
            res.AddEIds = eIds;
            res.AddLIds = lIds;
            return res;
        }

        public AddSurfaceRes AddSphere(OpenTK.Vector3d cPt, double r, int divCnt)
        {
            AddSurfaceRes res = new AddSurfaceRes();

            IList<uint> retVIds = new List<uint>();
            IList<uint> retEIds = new List<uint>();
            IList<uint> retLIds = new List<uint>();

            OpenTK.Vector3d ptU = new OpenTK.Vector3d(cPt.X, cPt.Y, cPt.Z + r);
            OpenTK.Vector3d ptD = new OpenTK.Vector3d(cPt.X, cPt.Y, cPt.Z - r);
            uint vIdU = AddVertex(CadElementType.Solid, 0, ptU).AddVId;
            uint vIdD = AddVertex(CadElementType.Solid, 0, ptD).AddVId;
            System.Diagnostics.Debug.Assert(vIdU != 0);
            System.Diagnostics.Debug.Assert(vIdD != 0);
            retVIds.Add(vIdU);
            retVIds.Add(vIdD);

            IList<IList<uint>> vIdss = new List<IList<uint>>();
            int divCntXY = divCnt;
            int divCntZY = divCnt / 2;
            for (int iXY = 0; iXY < divCntXY; iXY++)
            {
                double thetaXY = (iXY / (double)divCntXY) * 360.0 * (Math.PI / 180.0);
                IList<uint> vIds = new List<uint>();
                vIdss.Add(vIds);

                vIds.Add(vIdU);
                for (int iZY = 1; iZY <= (divCntZY - 1); iZY++)
                {
                    double thetaZY = (iZY / (double)divCntZY) * 180.0 * (Math.PI / 180);
                    double x = cPt.X + r * Math.Sin(thetaZY) * Math.Cos(thetaXY);
                    double y = cPt.Y + r * Math.Sin(thetaZY) * Math.Sin(thetaXY);
                    double z = cPt.Z + r * Math.Cos(thetaZY);
                    OpenTK.Vector3d pt = new OpenTK.Vector3d(x, y, z);
                    uint vId = AddVertex(CadElementType.Solid, 0, pt).AddVId;
                    System.Diagnostics.Debug.Assert(vId != 0);
                    vIds.Add(vId);
                    retVIds.Add(vId);
                }
                vIds.Add(vIdD);
            }
            {
                uint eId0 = 0; 
                for (int iXY = 0; iXY < vIdss.Count; iXY++)
                {
                    IList<uint> vIds = vIdss[iXY];
                    if (iXY != 0)
                    {
                        IList<uint> prevVIds = vIdss[iXY - 1];
                        for (int iZY = 1; iZY < (vIds.Count - 1); iZY++)
                        {
                            uint prevVId = prevVIds[iZY];
                            uint vId = vIds[iZY];
                            var res1 = ConnectVertexLine(prevVId, vId);
                            uint addEId = res1.AddEId;
                            uint addLId = res1.AddLId;
                            System.Diagnostics.Debug.Assert(addEId != 0);
                            System.Diagnostics.Debug.Assert(addLId == 0);
                            retEIds.Add(addEId);
                        }
                    }

                    for (int iZY = 0; iZY < (vIds.Count - 1); iZY++)
                    {
                        uint vId1 = vIds[iZY];
                        uint vId2 = vIds[iZY + 1];
                        var res1 = ConnectVertexLine(vId1, vId2);
                        uint addEId = res1.AddEId;
                        uint addLId = res1.AddLId;
                        System.Diagnostics.Debug.Assert(addEId != 0);
                        retEIds.Add(addEId);
                        if (iXY ==0)
                        {
                            System.Diagnostics.Debug.Assert(addLId == 0);
                        }
                        else
                        {
                            System.Diagnostics.Debug.Assert(addLId != 0);
                            retLIds.Add(addLId);
                        }
                    }

                    if (iXY == (vIdss.Count - 1))
                    {
                        IList<uint> nextVIds = vIdss[0];
                        for (int iZY = 1; iZY < (vIds.Count - 1); iZY++)
                        {
                            uint vId = vIds[iZY];
                            uint nextVId = nextVIds[iZY];
                            var res1 = ConnectVertexLine(vId, nextVId);
                            uint addEId = res1.AddEId;
                            uint addLId = res1.AddLId;
                            System.Diagnostics.Debug.Assert(addEId != 0);
                            retEIds.Add(addEId);
                            System.Diagnostics.Debug.Assert(addLId != 0);
                            retLIds.Add(addLId);

                            if (iZY == 1)
                            {
                                eId0 = addEId;
                            }
                        }
                    }
                }

                {
                    uint addLId = SealHole(eId0, false);
                    System.Diagnostics.Debug.Assert(addLId != 0);
                    retLIds.Add(addLId);
                }
            }

            res.AddVIds = retVIds;
            res.AddEIds = retEIds;
            res.AddLIds = retLIds;
            return res;
        }

        public AddPolygonRes AddRectLoop(uint lId, OpenTK.Vector2d sPt, OpenTK.Vector2d ePt)
        {
            var res = new AddPolygonRes();

            var p0 = new OpenTK.Vector2d(sPt.X, sPt.Y);
            var p1 = new OpenTK.Vector2d(ePt.X, sPt.Y);
            var p2 = new OpenTK.Vector2d(ePt.X, ePt.Y);
            var p3 = new OpenTK.Vector2d(sPt.X, ePt.Y);
            if ((sPt - ePt).X * (sPt - ePt).Y < 0)
            {
                p1 = new OpenTK.Vector2d(sPt.X, ePt.Y);
                p3 = new OpenTK.Vector2d(ePt.X, sPt.Y);
            }
            System.Diagnostics.Debug.Assert(IsElementId(CadElementType.Loop, lId));
            LoopEdgeItr lItr = BRep.GetLoopEdgeItr(lId);
            bool in0 = CheckIsPointInsideLoopEdgeItr(lItr, p0);
            bool in1 = CheckIsPointInsideLoopEdgeItr(lItr, p1);
            bool in2 = CheckIsPointInsideLoopEdgeItr(lItr, p2);
            bool in3 = CheckIsPointInsideLoopEdgeItr(lItr, p3);
            if (!in0)
            {
                System.Diagnostics.Debug.Assert(false);
                return res;
            }
            Loop3D l = GetLoop(lId);
            if (in1 && in2 && in3)
            {
                var pts = new List<OpenTK.Vector3d>();
                pts.Add(l.UnProject(p0));
                pts.Add(l.UnProject(p1));
                pts.Add(l.UnProject(p2));
                pts.Add(l.UnProject(p3));
                return AddPolygon(pts, lId);
            }
            if (in1 && !in2 && !in3)
            {
                uint eId7;
                bool isSameDir7;
                OpenTK.Vector2d p7;
                FindIntersectionEdge(lId, p0, p3, out eId7, out isSameDir7, out p7);
                uint vId7 = AddVertex(CadElementType.Edge, eId7, l.UnProject(p7)).AddVId;
                uint eId5;
                bool isSameDir5;
                OpenTK.Vector2d p5;
                FindIntersectionEdge(lId, p1, p2, out eId5, out isSameDir5, out p5);
                uint vId5 = AddVertex(CadElementType.Edge, eId5, l.UnProject(p5)).AddVId;
                uint vId0 = AddVertex(CadElementType.Loop, lId, l.UnProject(p0)).AddVId;
                uint vId1 = AddVertex(CadElementType.Loop, lId, l.UnProject(p1)).AddVId;
                ConnectVertexLine(vId7, vId0);
                ConnectVertexLine(vId0, vId1);
                ConnectVertexLine(vId1, vId5);
            }
            if (!in1 && !in2 && in3)
            {
                uint eId6;
                bool isSameDir6;
                OpenTK.Vector2d p6;
                FindIntersectionEdge(lId, p3, p2, out eId6, out isSameDir6, out p6);
                uint vId6 = AddVertex(CadElementType.Edge, eId6, l.UnProject(p6)).AddVId;
                uint eId4;
                bool isSameDir4;
                OpenTK.Vector2d p4;
                FindIntersectionEdge(lId, p0, p1, out eId4, out isSameDir4, out p4);
                uint vId4 = AddVertex(CadElementType.Edge, eId4, l.UnProject(p4)).AddVId;
                uint vId3 = AddVertex(CadElementType.Loop, lId, l.UnProject(p3)).AddVId;
                uint vId0 = AddVertex(CadElementType.Loop, lId, l.UnProject(p0)).AddVId;
                ConnectVertexLine(vId6, vId3);
                ConnectVertexLine(vId3, vId0);
                ConnectVertexLine(vId0, vId4);
            }
            if (!in1 && !in3)
            {
                uint eId7;
                bool isSameDir7;
                OpenTK.Vector2d p7;
                FindIntersectionEdge(lId, p0, p3, out eId7, out isSameDir7, out p7);
                uint vId7 = AddVertex(CadElementType.Edge, eId7, l.UnProject(p7)).AddVId;
                uint eId4;
                bool isSameDir4;
                OpenTK.Vector2d p4;
                FindIntersectionEdge(lId, p0, p1, out eId4, out isSameDir4, out p4);
                uint vId4 = AddVertex(CadElementType.Edge, eId4, l.UnProject(p4)).AddVId;
                uint vId0 = AddVertex(CadElementType.Loop, lId, l.UnProject(p0)).AddVId;
                ConnectVertexLine(vId7, vId0);
                ConnectVertexLine(vId0, vId4);
            }
            return res;
        }

        public AddSurfaceRes LiftLoop(uint lId, OpenTK.Vector3d dir)
        {
            AddSurfaceRes res = new AddSurfaceRes();

            if (!BRep.IsElementId(CadElementType.Loop, lId))
            {
                return res;
            }
            var edgePropss = new List<List<LiftEdgeProp>>();
            {
                OpenTK.Vector3d udir = new OpenTK.Vector3d(dir);
                udir.Normalize();
                for (LoopEdgeItr itr = BRep.GetLoopEdgeItr(lId); !itr.IsChildEnd; itr.ShiftChildLoop())
                {
                    int iiul = edgePropss.Count;
                    edgePropss.Add(new List<LiftEdgeProp>());
                    for (itr.Begin(); !itr.IsEnd(); itr.Next())
                    {
                        uint vId = itr.GetVertexId();
                        uint eId0;
                        bool isSameDir0;
                        itr.GetEdgeId(out eId0, out isSameDir0);
                        uint lId0 = BRep.GetEdgeLoopId(eId0, !isSameDir0);
                        System.Diagnostics.Debug.Assert(lId0 != lId);
                        System.Diagnostics.Debug.Assert(LoopArray.IsObjectId(lId0));
                        OpenTK.Vector3d n0 = LoopArray.GetObject(lId0).Normal;
                        bool isVert0 = Math.Abs(OpenTK.Vector3d.Dot(n0, dir)) < 1.0e-5;
                        edgePropss[iiul].Add(new LiftEdgeProp(vId, eId0, isSameDir0, isVert0));
                    }
                }
            }
            for (int iiul = 0; iiul < edgePropss.Count; iiul++)
            {
                IList<LiftEdgeProp> edgeProps = edgePropss[iiul];
                int npo = edgeProps.Count;
                for (int iiv = 0; iiv < npo; iiv++)
                {
                    int jiv = (iiv == 0) ? npo - 1 : iiv - 1;
                    uint outVId = edgeProps[iiv].OutVId;
                    if (!edgeProps[iiv].IsVert && !edgeProps[jiv].IsVert)
                    {
                        uint addVId = BRep.AddVertexToLoop(lId);
                        OpenTK.Vector3d outV = GetVertexCoord(outVId);
                        edgeProps[iiv].InVId = VertexArray.AddObject(addVId, new Vertex3D(outV + dir));
                        res.AddVIds.Add(addVId);
                        ////
                        var vItr0 = new VertexEdgeItr(BRep, outVId);
                        for (; vItr0.GetLoopId() != lId && !vItr0.IsEnd(); vItr0.Next())
                        {
                        }
                        System.Diagnostics.Debug.Assert(vItr0.GetLoopId() == lId);
                        var vItr1 = new VertexEdgeItr(BRep, edgeProps[iiv].InVId);
                        System.Diagnostics.Debug.Assert(vItr1.GetLoopId() == lId);
                        uint addEId0 = BRep.ConnectVertex(vItr0, vItr1, true).AddEId;
                        EdgeArray.AddObject(addEId0, new Edge3D());
                        res.AddEIds.Add(addEId0);
                    }
                    else if (edgeProps[iiv].IsVert && edgeProps[jiv].IsVert)
                    {
                        edgeProps[iiv].InVId = outVId;
                        Vertex3D ver = VertexArray.GetObject(edgeProps[iiv].InVId);
                        ver.Point += dir;
                    }
                    else
                    {
                        uint eIdInsert;
                        {
                            var vItr = new VertexEdgeItr(BRep, edgeProps[iiv].OutVId);
                            for (; vItr.GetLoopId() != lId && !vItr.IsEnd(); vItr.Next())
                            {
                            }
                            System.Diagnostics.Debug.Assert(vItr.GetLoopId() == lId);
                            uint eId0;
                            bool isSameDir0;
                            vItr.GetAheadEdgeId(out eId0, out isSameDir0);
                            uint eId1;
                            bool isSameDir1;
                            vItr.GetBehindEdgeId(out eId1, out isSameDir1);
                            if (edgeProps[iiv].IsVert)
                            {
                                eIdInsert = eId1;
                            }
                            else
                            {
                                System.Diagnostics.Debug.Assert(edgeProps[jiv].IsVert);
                                eIdInsert = eId0;
                            }
                        }
                        uint addVId = BRep.AddVertexToEdge(eIdInsert);
                        OpenTK.Vector3d outV = GetVertexCoord(outVId);
                        edgeProps[iiv].InVId = VertexArray.AddObject(addVId, new Vertex3D(outV + dir));
                        res.AddVIds.Add(addVId);
                        /////
                        uint addEId;
                        {
                            var vItr = new VertexEdgeItr(BRep, edgeProps[iiv].InVId);
                            uint eId0;
                            bool isSameDir0;
                            vItr.GetAheadEdgeId(out eId0, out isSameDir0);
                            uint eId1;
                            bool isSameDir1;
                            vItr.GetBehindEdgeId(out eId1, out isSameDir1);
                            addEId = (eId0 == eIdInsert) ? eId1 : eId0;
                            System.Diagnostics.Debug.Assert(addEId != eIdInsert);
                        }
                        EdgeArray.AddObject(addEId, new Edge3D());
                        res.AddEIds.Add(addEId);
                    }
                }
            }

            for (int iiul = 0; iiul < edgePropss.Count; iiul++)
            {
                IList<LiftEdgeProp> edgeProps = edgePropss[iiul];
                int npo = edgeProps.Count;
                for (int iiv = 0; iiv < npo; iiv++)
                {
                    if (edgeProps[iiv].IsVert)
                    {
                        continue;
                    }
                    int jiv = (iiv == npo - 1) ? 0 : iiv + 1;
                    var vItr0 = new VertexEdgeItr(BRep, edgeProps[iiv].InVId);
                    for (; vItr0.GetLoopId() != lId && !vItr0.IsEnd(); vItr0.Next())
                    {
                    }
                    System.Diagnostics.Debug.Assert(vItr0.GetLoopId() == lId);
                    var vItr1 = new VertexEdgeItr(BRep, edgeProps[jiv].InVId);
                    for (; vItr1.GetLoopId() != lId && !vItr1.IsEnd(); vItr1.Next())
                    {
                    }
                    System.Diagnostics.Debug.Assert(vItr1.GetLoopId() == lId);
                    var cvRes = BRep.ConnectVertex(vItr0, vItr1, iiul == 0);
                    uint addEId0 = cvRes.AddEId;
                    EdgeArray.AddObject(addEId0, new Edge3D());
                    res.AddEIds.Add(addEId0);
                    uint addLId = BRep.GetEdgeLoopId(addEId0, cvRes.IsLeftAddL);
                    System.Diagnostics.Debug.Assert(addLId != lId);
                    {
                        OpenTK.Vector3d v0 = GetVertexCoord(edgeProps[iiv].OutVId);
                        OpenTK.Vector3d v1 = GetVertexCoord(edgeProps[jiv].OutVId);
                        var x = v1 - v0;
                        x.Normalize();
                        var n = OpenTK.Vector3d.Cross(x, dir);
                        n.Normalize();
                        LoopArray.AddObject(addLId, new Loop3D(v0, n, x));
                        res.AddLIds.Add(addLId);
                    }
                }
            }
            {
                Loop3D l = LoopArray.GetObject(lId);
                l.Origin += dir;
            }
            System.Diagnostics.Debug.Assert(AssertValid() == 0);
            return res;
        }

        public uint SealHole(uint eId, bool isLeft)
        {
            uint addLId = BRep.SealHole(eId, isLeft);

            if (addLId != 0)
            {
                OpenTK.Vector3d origin;
                OpenTK.Vector3d normal;
                OpenTK.Vector3d xdir;
                CalcLoopNormal(addLId, out origin, out normal, out xdir);

                LoopArray.AddObject(addLId, new Loop3D(origin, normal, xdir));
            }
            System.Diagnostics.Debug.Assert(AssertValid() == 0);
            return addLId;
        }

        // TODO:
        // Radial EdgeでLoopを作ったとき、
        // LoopEdgeItrは特に変更点なし
        // VertexEdgeItrはなんらかの対応が必要か？(edgeが増えるわけではないので問題ない？)
        public uint MakeRadialLoop(IList<uint> eIds)
        {
            //-------------------------------------------------------------
            // HalfEdgeの方向を決める
            IList<bool> isSameDirs = new List<bool>();
            {
                int edgeCnt = eIds.Count;
                uint firstVId = 0;
                uint prevVId = 0;
                for (int i = 0; i < edgeCnt; i++)
                {
                    uint eId = eIds[i];
                    uint sVId;
                    uint eVId;
                    GetEdgeVertexId(eId, out sVId, out eVId);
                    System.Diagnostics.Debug.Assert(sVId != 0);
                    System.Diagnostics.Debug.Assert(eVId != 0);
                    bool isSameDir = true;
                    if (firstVId == 0)
                    {
                        firstVId = sVId;
                        prevVId = eVId;
                        isSameDir = true;
                    }
                    else
                    {
                        if (prevVId == sVId)
                        {
                            prevVId = eVId;
                            isSameDir = true;
                        }
                        else if (prevVId == eVId)
                        {
                            prevVId = sVId;
                            isSameDir = false;
                        }
                        else
                        {
                            // eIdsの指定順番がループの周回になっていない
                            System.Diagnostics.Debug.Assert(false);
                            return 0;
                        }
                    }
                    isSameDirs.Add(isSameDir);
                }
                System.Diagnostics.Debug.Assert(prevVId == firstVId);
            }
            //-------------------------------------------------------------

            uint addLId = BRep.MakeRadialLoop(eIds, isSameDirs);

            if (addLId != 0)
            {
                OpenTK.Vector3d origin;
                OpenTK.Vector3d normal;
                OpenTK.Vector3d xdir;
                CalcLoopNormal(addLId, out origin, out normal, out xdir);

                LoopArray.AddObject(addLId, new Loop3D(origin, normal, xdir));
            }
            System.Diagnostics.Debug.Assert(AssertValid() == 0);
            return addLId;
        }

        public uint AddSolid(IList<uint> lIds, IList<OpenTK.Vector3d> holes, IList<uint> insideVIds)
        {
            Solid addSolid = new Solid(lIds, holes, insideVIds);
            uint addSId = SolidArray.GetFreeObjectId();

            SolidArray.AddObject(addSId, addSolid);
            return addSId;
        }

        private void CalcLoopNormal(
            uint lId, out OpenTK.Vector3d origin, out OpenTK.Vector3d normal, out OpenTK.Vector3d xdir)
        {
            origin = new OpenTK.Vector3d();
            normal = new OpenTK.Vector3d();
            xdir = new OpenTK.Vector3d();

            IList<uint> vIds = new List<uint>();
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
                    //Edge2D e2D = GetEdge2D(eId0, lId);
                    uint sVId = GetEdgeVertexId(eId0, true);
                    uint eVId = GetEdgeVertexId(eId0, false);
                    if (!vIds.Contains(sVId))
                    {
                        vIds.Add(sVId);
                    }
                    if (!vIds.Contains(eVId))
                    {
                        vIds.Add(eVId);
                    }
                }

                uint vId1 = vIds[0];
                uint vId2 = vIds[1];
                uint vId3 = vIds[2];
                OpenTK.Vector3d pt1 = GetVertexCoord(vId1);
                OpenTK.Vector3d pt2 = GetVertexCoord(vId2);
                OpenTK.Vector3d pt3 = GetVertexCoord(vId3);
                origin = pt1;
                xdir = pt2 - pt1;
                normal = CadUtils3D.TriNormal(pt1, pt2, pt3);
            }
        }

        private bool FindIntersectionEdge(
            uint lId, OpenTK.Vector2d sPt, OpenTK.Vector2d ePt,
            out uint eIdNearest, out bool isSameDirNearest, out OpenTK.Vector2d ptNearest)
        {
            eIdNearest = 0;
            isSameDirNearest = false;
            ptNearest = new OpenTK.Vector2d();

            OpenTK.Vector2d dir = ePt - sPt;
            dir.Normalize();
            double sqdist = -1;
            for (LoopEdgeItr pItr = BRep.GetLoopEdgeItr(lId); !pItr.IsChildEnd; pItr.ShiftChildLoop())
            {
                for (; !pItr.IsEnd(); pItr.Next())
                {
                    uint eId;
                    bool isSameDir;
                    pItr.GetEdgeId(out eId, out isSameDir);
                    Edge2D edge = GetEdge2D(eId, lId);
                    OpenTK.Vector2d sec;
                    if (!edge.GetNearestIntersectionPointAgainstHalfLine(out sec, sPt, dir))
                    {
                        continue;
                    }
                    if (sqdist < 0 || (sec - sPt).LengthSquared < sqdist)
                    {
                        eIdNearest = eId;
                        isSameDirNearest = isSameDir;
                        ptNearest = sec;
                        sqdist = (sec - sPt).LengthSquared;
                    }
                }
            }
            if (!IsElementId(CadElementType.Edge, eIdNearest))
            {
                return false;
            }
            return true;
        }

        private bool CheckIsPointInsideLoopEdgeItr(LoopEdgeItr lItr, OpenTK.Vector2d point)
        {
            for (uint i = 1; i < 29; i++)
            { // 29 is handy prim number
                int crossCounter = 0;
                bool iflg = true;
                var dir = new OpenTK.Vector2d(Math.Sin(6.28 * i / 29.0), Math.Cos(6.28 * i / 29.0));
                for (lItr.Begin(); !lItr.IsEnd(); lItr.Next())
                {
                    uint eId;
                    bool isSameDir;
                    lItr.GetEdgeId(out eId, out isSameDir);
                    if (eId == 0)
                    {
                        return false;
                    }
                    System.Diagnostics.Debug.Assert(EdgeArray.IsObjectId(eId));
                    Edge2D e = GetEdge2D(eId, lItr.GetLoopId());
                    int ires = e.NumIntersectAgainstHalfLine(point, dir);
                    if (ires == -1)
                    {
                        iflg = false;
                        break;
                    }
                    crossCounter += ires;
                }
                if (iflg == true)
                {
                    if (crossCounter % 2 == 0) return false;
                    return true;
                }
            }

            // fail
            System.Diagnostics.Debug.Assert(false);
            return false;
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
            else if (type == CadElementType.Solid)
            {
                return SolidArray.IsObjectId(id);
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
            else if (type == CadElementType.Solid)
            {
                return SolidArray.GetObjectIds();
            }

            System.Diagnostics.Debug.Assert(false);
            IList<uint> nullVec = new List<uint>();
            return nullVec;
        }

        public Vertex3D GetVertex(uint vId)
        {
            if (!BRep.IsElementId(CadElementType.Vertex, vId))
            {
                return null;
            }
            if (!VertexArray.IsObjectId(vId))
            {
                return null;
            }
            Vertex3D v = VertexArray.GetObject(vId);
            return v;
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
                    Edge2D e2D = GetEdge2D(eId0, lId);
                    l.Edges.Add(new KeyValuePair<Edge2D, bool>(e2D, isSameDir0));
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

        public double[] GetVertexColor(uint vId)
        {
            double[] color = new double[3];
            if (!BRep.IsElementId(CadElementType.Vertex, vId))
            {
                System.Diagnostics.Debug.Assert(false);
                return color;
            }
            if (!VertexArray.IsObjectId(vId))
            {
                System.Diagnostics.Debug.Assert(false);
                return color;
            }
            Vertex3D v = VertexArray.GetObject(vId);
            v.Color.CopyTo(color, 0);
            return color;
        }

        public Solid GetSolid(uint sId)
        {
            System.Diagnostics.Debug.Assert(IsElementId(CadElementType.Solid, sId));
            Solid s = SolidArray.GetObject(sId);
            return s;
        }

        public bool SetVertexColor(uint vId, double[] color)
        {
            if (!BRep.IsElementId(CadElementType.Vertex, vId))
            {
                System.Diagnostics.Debug.Assert(false);
                return false;
            }
            if (!VertexArray.IsObjectId(vId))
            {
                System.Diagnostics.Debug.Assert(false);
                return false;
            }
            Vertex3D v = VertexArray.GetObject(vId);
            color.CopyTo(v.Color, 0);
            return true;
        }

        public double[] GetEdgeColor(uint eId)
        {
            double[] color = new double[3];
            if (!BRep.IsElementId(CadElementType.Edge, eId))
            {
                System.Diagnostics.Debug.Assert(false);
                return color;
            }
            if (!EdgeArray.IsObjectId(eId))
            {
                System.Diagnostics.Debug.Assert(false);
                return color;
            }
            Edge3D e = EdgeArray.GetObject(eId);
            e.Color.CopyTo(color, 0);
            return color;
        }

        public bool SetEdgeColor(uint eId, double[] color)
        {
            if (!BRep.IsElementId(CadElementType.Edge, eId))
            {
                System.Diagnostics.Debug.Assert(false);
                return false;
            }
            if (!EdgeArray.IsObjectId(eId))
            {
                System.Diagnostics.Debug.Assert(false);
                return false;
            }
            Edge3D e = EdgeArray.GetObject(eId);
            color.CopyTo(e.Color, 0);
            return true;
        }

        public double[] GetLoopColor(uint lId)
        {
            double[] color = new double[3];
            if (!LoopArray.IsObjectId(lId))
            {
                System.Diagnostics.Debug.Assert(false);
                return color;
            }
            Loop3D l = LoopArray.GetObject(lId);
            color = new double[l.Color.Length];
            l.Color.CopyTo(color, 0);
            return color;
        }

        public bool SetLoopColor(uint lId, double[] color)
        {
            if (!LoopArray.IsObjectId(lId))
            {
                System.Diagnostics.Debug.Assert(false);
                return false;
            }
            Loop3D l = LoopArray.GetObject(lId);
            l.Color = new double[color.Length];
            color.CopyTo(l.Color, 0);
            return true;
        }

        private int AssertValid()
        {
            if (!BRep.AssertValid())
            {
                System.Diagnostics.Debug.Assert(false);
                return 6;
            }
            {
                IList<uint> lIds = BRep.GetElementIds(CadElementType.Loop);
                for (int i = 0; i < lIds.Count; i++)
                {
                    if (!LoopArray.IsObjectId(lIds[i]))
                    {
                        System.Diagnostics.Debug.Assert(false);
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
                        System.Diagnostics.Debug.Assert(false);
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
                        System.Diagnostics.Debug.Assert(false);
                        return 9;
                    }
                }
            }
            return 0;
        }
    }
}
