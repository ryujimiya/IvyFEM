using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Cad3DToXY : Cad2D
    {
        public OpenTK.Vector3d Origin { get; set; }
        public OpenTK.Vector3d Normal { get; set; }
        public OpenTK.Vector3d XDir { get; set; }
        private Dictionary<uint, uint> VId3D2D  = new Dictionary<uint, uint>();
        private Dictionary<uint, uint> EId3D2D  = new Dictionary<uint, uint>();
        private Dictionary<uint, uint> VId2D3D = new Dictionary<uint, uint>();
        private Dictionary<uint, uint> EId2D3D = new Dictionary<uint, uint>();

        public Cad3DToXY() : base()
        {
            Origin = new OpenTK.Vector3d(0.0, 0.0, 0.0);
            Normal = new OpenTK.Vector3d(0.0, 0.0, 1.0);
            XDir = new OpenTK.Vector3d(1.0, 0.0, 0.0);
        }

        public Cad3DToXY(Cad3D cad, CadElementType type, uint id) : base()
        {
            ToXY(cad, type, id);
        }

        public Cad3DToXY(Cad3DToXY src)
        {
            Copy(src);
        }

        public void Copy(Cad3DToXY src)
        {
            base.Copy(src);

            Origin = new OpenTK.Vector3d(src.Origin.X, src.Origin.Y, src.Origin.Z);
            Normal = new OpenTK.Vector3d(src.Normal.X, src.Normal.Y, src.Normal.Z);
            XDir = new OpenTK.Vector3d(src.XDir.X, src.XDir.Y, src.XDir.Z);

            VId3D2D  = new Dictionary<uint, uint>(src.VId3D2D);
            EId3D2D = new Dictionary<uint, uint>(src.EId3D2D);
            VId2D3D = new Dictionary<uint, uint>(src.VId2D3D);
            EId2D3D = new Dictionary<uint, uint>(src.EId2D3D);
        }

        public override void Clear()
        {
            base.Clear();

            VId3D2D = new Dictionary<uint, uint>();
            EId3D2D = new Dictionary<uint, uint>();
            VId2D3D = new Dictionary<uint, uint>();
            EId2D3D = new Dictionary<uint, uint>();
        }

        public OpenTK.Vector2d Project(OpenTK.Vector3d p)
        {
            double x = OpenTK.Vector3d.Dot((p - Origin), XDir);
            double y = OpenTK.Vector3d.Dot((p - Origin), OpenTK.Vector3d.Cross(Normal, XDir));
            return new OpenTK.Vector2d(x, y);
        }

        public OpenTK.Vector3d UnProject(OpenTK.Vector2d p)
        {
            return Origin + XDir * p.X + OpenTK.Vector3d.Cross(Normal, XDir) * p.Y;
        }

        public uint GetVertexId3DFrom2D(uint vId)
        {
            return VId2D3D[vId];
        }

        public uint GetEdgeId3DFrom2D(uint eId)
        {
            return EId2D3D[eId];
        }

        public uint GetVertexId2DFrom3D(uint vId)
        {
            return VId3D2D[vId];
        }

        public uint GetEdgeId2DFrom3D(uint eId)
        {
            if (!EId3D2D.ContainsKey(eId))
            {
                return 0;
            }
            return EId3D2D[eId];
        }

        public void ToXY(Cad3D cad, CadElementType type, uint id)
        {
            Clear();

            if (type == CadElementType.Loop)
            {
                Loop3DToXY(cad, id);
            }
            else if (type == CadElementType.Edge)
            {
                Edge3DToXY(cad, id);
            }
            else
            {
                // 対応不要
                throw new NotImplementedException();
            }
        }

        private OpenTK.Vector3d GetYDirFromXDir(OpenTK.Vector3d xDir)
        {
            OpenTK.Vector3d yDir;
            double[] co = { XDir.X, XDir.Y, XDir.Z };
            double min = double.MaxValue;
            int minIndex = -1;
            for (int i = 0; i < 3; i++)
            {
                double abs = Math.Abs(co[i]);
                if (abs < min)
                {
                    min = abs;
                    minIndex = i;
                }
            }
            yDir = new OpenTK.Vector3d();
            if (minIndex == 2)
            {
                yDir = new OpenTK.Vector3d(-XDir.Y, XDir.X, 0.0);
                yDir = OpenTK.Vector3d.Normalize(yDir);
            }
            else if (minIndex == 1)
            {
                yDir = new OpenTK.Vector3d(-XDir.Z, 0.0, XDir.X);
                yDir = OpenTK.Vector3d.Normalize(yDir);
            }
            else if (minIndex == 0)
            {
                yDir = new OpenTK.Vector3d(0.0, XDir.Z, -XDir.Y);
                yDir = OpenTK.Vector3d.Normalize(yDir);
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            // check
            double dotXY = OpenTK.Vector3d.Dot(XDir, yDir);
            System.Diagnostics.Debug.Assert(Math.Abs(dotXY) < Constants.PrecisionLowerLimit);
            return yDir;
        }

        private void Edge3DToXY(Cad3D cad, uint eId)
        {
            System.Diagnostics.Debug.Assert(cad.IsElementId(CadElementType.Edge, eId));
            uint sVId;
            uint eVId;
            cad.GetEdgeVertexId(eId, out sVId, out eVId);
            OpenTK.Vector3d sPt = cad.GetVertexCoord(sVId);
            OpenTK.Vector3d ePt = cad.GetVertexCoord(eVId);

            //Origin = new OpenTK.Vector3d(0.0, 0.0, 0.0);
            //XDir = new OpenTK.Vector3d(1.0, 0.0, 0.0);
            //Normal = new OpenTK.Vector3d(0.0, 0.0, 1.0);
            //
            // 辺を含む平面をXY平面とする直角座標系を求める
            Origin = new OpenTK.Vector3d(sPt.X, sPt.Y, sPt.Z);
            XDir = OpenTK.Vector3d.Normalize(ePt - sPt);
            OpenTK.Vector3d yDir = GetYDirFromXDir(XDir);
            Normal = OpenTK.Vector3d.Cross(XDir, yDir);

            uint[] vIds = { sVId, eVId };
            uint[] vId2Ds = { 0, 0 };

            for (int i = 0; i < 2; i++)
            {
                uint vId = vIds[0];
                OpenTK.Vector3d vec = cad.GetVertexCoord(vId);
                OpenTK.Vector2d vec2D = Project(vec);
                uint vId2D;
                if (VId3D2D.ContainsKey(vId))
                {
                    vId2D = VId3D2D[vId];
                }
                else
                {
                    vId2D = AddVertex(CadElementType.Loop, 0, vec2D).AddVId;
                    System.Diagnostics.Debug.Assert(vId2D != 0);
                    VId3D2D[vId] = vId2D;
                    VId2D3D[vId2D] = vId;
                }
                vId2Ds[i] = vId;
            }
            uint eId2D;
            if (EId3D2D.ContainsKey(eId))
            {
                eId2D = EId3D2D[eId];
            }
            else
            {
                eId2D = ConnectVertexLine(vId2Ds[0], vId2Ds[1]).AddEId;
                System.Diagnostics.Debug.Assert(eId2D != 0);
                EId3D2D[eId] = eId2D;
                EId2D3D[eId2D] = eId;
            }
        }

        private void Loop3DToXY(Cad3D cad, uint lId)
        {
            System.Diagnostics.Debug.Assert(cad.IsElementId(CadElementType.Loop, lId));
            Loop3D loop = cad.GetLoop(lId);
            Origin = new OpenTK.Vector3d(loop.Origin.X, loop.Origin.Y, loop.Origin.Z);
            Normal = new OpenTK.Vector3d(loop.Normal.X, loop.Normal.Y, loop.Normal.Z);
            XDir = new OpenTK.Vector3d(loop.XDir.X, loop.XDir.Y, loop.XDir.Z);

            int loopCounter = 0;
            uint parentLId = 0;
            for (LoopEdgeItr lItr = cad.GetLoopEdgeItr(lId); !lItr.IsChildEnd; lItr.ShiftChildLoop(), loopCounter++)
            {
                for (lItr.Begin(); !lItr.IsEnd(); lItr.Next())
                {
                    //uint vId = lItr.GetVertexId();

                    uint eId;
                    bool isSameDir;
                    lItr.GetEdgeId(out eId, out isSameDir);
                    if (!cad.IsElementId(CadElementType.Edge, eId))
                    {
                        continue;
                    }
                    uint sVId;
                    uint eVId;
                    cad.GetEdgeVertexId(eId, out sVId, out eVId);
                    uint[] vIds = { sVId, eVId };
                    uint[] vId2Ds = new uint[2];
                    OpenTK.Vector2d[] vec2Ds = new OpenTK.Vector2d[2];

                    for (int i = 0; i < 2; i++)
                    {
                        uint vId = vIds[i];
                        OpenTK.Vector3d vec = cad.GetVertexCoord(vId);
                        OpenTK.Vector2d vec2D = Project(vec);
                        uint vId2D;
                        if (VId3D2D.ContainsKey(vId))
                        {
                            vId2D = VId3D2D[vId];
                        }
                        else
                        {
                            vId2D = AddVertex(CadElementType.Loop, parentLId, vec2D).AddVId;
                            System.Diagnostics.Debug.Assert(vId2D != 0);
                            VId3D2D[vId] = vId2D;
                            VId2D3D[vId2D] = vId;
                        }
                        vId2Ds[i] = vId2D;
                        vec2Ds[i] = vec2D;
                    }

                    uint eId2D;
                    if (EId3D2D.ContainsKey(eId))
                    {
                        eId2D = EId3D2D[eId];
                    }
                    else
                    {
                        var res = ConnectVertexLine(vId2Ds[0], vId2Ds[1]);
                        eId2D = res.AddEId;
                        if (loopCounter == 0 && res.AddLId != 0)
                        {
                            parentLId = res.AddLId;
                        }
                        //System.Diagnostics.Debug.Assert(eId2D != 0);
                        if (eId2D != 0)
                        {
                            EId3D2D[eId] = eId2D;
                            EId2D3D[eId2D] = eId;
                        }
                        else
                        {
                            // 失敗
                            System.Diagnostics.Debug.Assert(false);
                        }
                    }
                }
            }
        }
    }
}
