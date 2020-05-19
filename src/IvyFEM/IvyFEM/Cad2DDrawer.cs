using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    public class Cad2DDrawer : IDrawer
    {
        public RotMode SutableRotMode { get; private set; } = RotMode.RotMode2D;
        public bool IsAntiAliasing { get; set; } = false;

        private double[] SelectedColor = { 1.0, 0.5, 1.0 }; //{ 1.0, 1.0, 0.0 };
        private byte[] Mask = new byte[128];
        private IList<Cad2DDrawPart> DrawParts = new List<Cad2DDrawPart>();
        private VertexArray VertexArray = new VertexArray();
        public uint LineWidth { get; set; } = 3;
        private uint PointSize = 5;
        private double TexCentX = 0;
        private double TexCentY = 0;
        private double TexScale = 1;

        public Cad2DDrawer()
        {
            SetupMask();
        }

        public Cad2DDrawer(Cad2D cad)
        {
            SetupMask();
            UpdateCadTopologyGeometry(cad);
        }

        private void SetupMask()
        {
            for (uint j = 0; j < 4; j++)
            {
                for (uint i = 0; i < 8; i++)
                {
                    Mask[j * 32 + i] = 0x33;
                }
                for (uint i = 0; i < 8; i++)
                {
                    Mask[j * 32 + 8 + i] = 0xcc;
                }
                for (uint i = 0; i < 8; i++)
                {
                    Mask[j * 32 + 16 + i] = 0x33;
                }
                for (uint i = 0; i < 8; i++)
                {
                    Mask[j * 32 + 24 + i] = 0xcc;
                }
            }
        }

        public bool UpdateCadTopologyGeometry(Cad2D cad)
        {
            SutableRotMode = RotMode.RotMode2D;
            IList<Cad2DDrawPart> oldDrawParts = new List<Cad2DDrawPart>();
            for (int i = 0; i < DrawParts.Count; i++)
            {
                oldDrawParts.Add(new Cad2DDrawPart(DrawParts[i]));
            }

            for (int idp = 0; idp < oldDrawParts.Count; idp++)
            {
                oldDrawParts[idp].MeshId = 0;
                oldDrawParts[idp].ShowMode = 0;
            }
            DrawParts.Clear();

            int minLayer;
            int maxLayer;
            cad.GetLayerMinMax(out minLayer, out maxLayer);
            double layerHeight = 1.0 / (maxLayer - minLayer + 1);

            {
                // 面をセット
                IList<uint> lIds = cad.GetElementIds(CadElementType.Loop);
                for (int iLId = 0; iLId < lIds.Count; iLId++)
                {
                    uint lId = lIds[iLId];
                    double height = 0;
                    {
                        int layer = cad.GetLayer(CadElementType.Loop, lId);
                        height = (layer - minLayer) * layerHeight;
                    }
                    int idp0 = 0;
                    for (; idp0 < oldDrawParts.Count; idp0++)
                    {
                        Cad2DDrawPart olddp = oldDrawParts[idp0];
                        if (olddp.Type == CadElementType.Loop && olddp.CadId == lId)
                        {
                            olddp.MeshId = 1;
                            olddp.Height = height;
                            double[] color = cad.GetLoopColor(lId);
                            color.CopyTo(olddp.Color, 0);
                            DrawParts.Add(oldDrawParts[idp0]);
                            break;
                        }
                    }
                    if (idp0 == oldDrawParts.Count)
                    {
                        Cad2DDrawPart dp = new Cad2DDrawPart();
                        dp.CadId = lId;
                        dp.Type = CadElementType.Loop;
                        dp.Height = height;
                        double[] color = cad.GetLoopColor(lId);
                        color.CopyTo(dp.Color, 0);
                        DrawParts.Add(dp);
                    }
                }
            }

            {
                // set edge
                IList<uint> eIds = cad.GetElementIds(CadElementType.Edge);
                for (int iEId = 0; iEId < eIds.Count; iEId++)
                {
                    uint eId = eIds[iEId];
                    double height = 0;
                    {
                        int layer = cad.GetLayer(CadElementType.Edge, eId);
                        height += (layer - minLayer + 0.01) * layerHeight;
                    }
                    int idp0 = 0;
                    for (; idp0 < oldDrawParts.Count; idp0++)
                    {
                        Cad2DDrawPart olddp = oldDrawParts[idp0];
                        if (olddp.Type == CadElementType.Edge && olddp.CadId == eId)
                        {
                            olddp.MeshId = 1;
                            olddp.Height = height;
                            DrawParts.Add(olddp);
                            break;
                        }
                    }
                    if (idp0 == oldDrawParts.Count)
                    {
                        Cad2DDrawPart dp = new Cad2DDrawPart();
                        dp.CadId = eId;
                        dp.Type = CadElementType.Edge;
                        dp.Height = height;
                        DrawParts.Add(dp);
                    }
                    {
                        Cad2DDrawPart dp = DrawParts[DrawParts.Count - 1];
                        Edge2D edge = cad.GetEdge(eId);
                        dp.CtrlPoints.Clear();
                        dp.CurveType = edge.CurveType;
                        if (edge.CurveType == CurveType.CurveArc)
                        {
                            OpenTK.Vector2d cPt;
                            double radius;
                            edge.GetCenterRadius(out cPt, out radius);
                            dp.CtrlPoints.Add(cPt);
                        }
                        else if (edge.CurveType == CurveType.CurveBezier)
                        {
                            IList<OpenTK.Vector2d> cos = edge.GetCurvePoint();
                            dp.CtrlPoints.Add(cos[0]);
                            dp.CtrlPoints.Add(cos[1]);
                        }
                    }
                }
            }

            { 
                // set vertex
                IList<uint> vIds = cad.GetElementIds(CadElementType.Vertex);
                for (int iVId = 0; iVId < vIds.Count; iVId++)
                {
                    uint vCadId = vIds[iVId];
                    int layer = cad.GetLayer(CadElementType.Vertex, vCadId);
                    double height = (layer - minLayer + 0.1) * layerHeight;
                    int idp0 = 0;
                    for (; idp0 < oldDrawParts.Count; idp0++)
                    {
                        Cad2DDrawPart olddp = oldDrawParts[idp0];
                        if (olddp.Type == CadElementType.Vertex && olddp.CadId == vCadId)
                        {
                            olddp.MeshId = 1;
                            olddp.Height = height;
                            DrawParts.Add(olddp);
                            break;
                        }
                    }
                    if (idp0 == oldDrawParts.Count)
                    {
                        Cad2DDrawPart dp = new Cad2DDrawPart();
                        dp.CadId = vCadId;
                        dp.Type = CadElementType.Vertex;
                        dp.Height = height;
                        DrawParts.Add(dp);
                    }
                }
            }

            oldDrawParts.Clear();

            UpdateCadGeometry(cad);
            return true;
        }

        public void UpdateCadGeometry(Cad2D cad)
        {
            Mesher2D mesh = new Mesher2D(cad);
            for (int idp = 0; idp < DrawParts.Count; idp++)
            {
                Cad2DDrawPart dp = DrawParts[idp];
                dp.Clear();
                uint cadId = dp.CadId;
                CadElementType cadType = dp.Type;
                if (!cad.IsElementId(cadType, cadId))
                {
                    continue;
                }
                uint meshId = mesh.GetIdFromCadId(cadId, cadType);
                if (meshId == 0)
                {
                    continue;
                }
                MeshType meshType;
                uint elemCnt;
                int loc;
                uint cadId0;
                mesh.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadId0);
                System.Diagnostics.Debug.Assert(cadId0 == cadId);
                if (meshType == MeshType.Tri)
                {
                    dp.SetTriArray(mesh.GetTriArrays()[loc]);
                    double[] color = cad.GetLoopColor(cadId0);
                    color.CopyTo(dp.Color, 0);
                }
                else if (meshType == MeshType.Bar)
                {
                    dp.SetBarArray(mesh.GetBarArrays()[loc]);
                    System.Diagnostics.Debug.Assert(cadType == CadElementType.Edge);
                    Edge2D edge = cad.GetEdge(cadId);
                    dp.CurveType = edge.CurveType;
                    dp.CtrlPoints.Clear();
                    // 2019-03-11 エッジの色 FIX
                    double[] color = edge.Color;
                    color.CopyTo(dp.Color, 0);
                    if (edge.CurveType == CurveType.CurveArc)
                    {
                        OpenTK.Vector2d cPt;
                        double radius;
                        edge.GetCenterRadius(out cPt, out radius);
                        dp.CtrlPoints.Add(cPt);
                    }
                    else if (edge.CurveType == CurveType.CurveBezier)
                    {
                        IList<OpenTK.Vector2d> cos = edge.GetCurvePoint();
                        dp.CtrlPoints.Add(cos[0]);
                        dp.CtrlPoints.Add(cos[1]);
                    }
                }
                else if (meshType == MeshType.Vertex)
                {
                    dp.SetVertex(mesh.GetVertexs()[loc]);
                }
            }

            {
                // 座標をセット
                IList<OpenTK.Vector2d> vecs = mesh.GetVectors();
                uint ptCnt = (uint)vecs.Count;
                uint ndim = 2;
                VertexArray.SetSize(ptCnt, ndim);
                for (int iPt = 0; iPt < ptCnt; iPt++)
                {
                    VertexArray.VertexCoords[iPt * ndim] = vecs[iPt].X;
                    VertexArray.VertexCoords[iPt * ndim + 1] = vecs[iPt].Y;
                }
                if (VertexArray.UVCoords != null)
                {
                    for (int iPt = 0; iPt < ptCnt; iPt++)
                    {
                        VertexArray.UVCoords[iPt * ndim] = vecs[iPt].X * TexScale;
                        VertexArray.UVCoords[iPt * ndim + 1] = vecs[iPt].Y * TexScale;
                    }
                }
            }
        }

        public BoundingBox3D GetBoundingBox(OpenTK.Matrix3d rot)
        {
            return VertexArray.GetBoundingBox(rot);
        }

        public void Draw()
        {
            GL.Enable(EnableCap.DepthTest);
            GL.Disable(EnableCap.CullFace);
            bool isLighting = GL.IsEnabled(EnableCap.Lighting);
            bool isTexture = GL.IsEnabled(EnableCap.Texture2D);
            bool isBlend = GL.IsEnabled(EnableCap.Blend);
            GL.Disable(EnableCap.Lighting);

            uint ndim = VertexArray.Dimension;

            ////////////////////////////////////////////////////////////////
            // モデルの描画

            /////////////
            // vertex arrayを登録する
            GL.EnableClientState(ArrayCap.VertexArray);
            GL.VertexPointer((int)ndim, VertexPointerType.Double, 0, VertexArray.VertexCoords);
            if (isTexture && VertexArray.UVCoords != null)
            {
                GL.EnableClientState(ArrayCap.TextureCoordArray);
                GL.TexCoordPointer(2, TexCoordPointerType.Double, 0, VertexArray.UVCoords);
                GL.MatrixMode(MatrixMode.Texture);
                GL.LoadIdentity();
                GL.Translate(-TexCentX, -TexCentY, 0.0);
            }
            GL.PointSize(PointSize);
            GL.LineWidth(LineWidth);
            for (int idp = 0; idp < DrawParts.Count; idp++)
            {
                Cad2DDrawPart dp = DrawParts[idp];
                if (dp.ShowMode == -1)
                {
                    continue;
                }
                double height = dp.Height;
                double dispX = dp.DispX;
                double dispY = dp.DispY;
                if (dp.Type == CadElementType.Vertex)
                {
                    GL.Disable(EnableCap.Texture2D);
                    if (dp.ShowMode == 1)
                    {
                        GL.Color3(1.0, 1.0, 0.0);
                    }
                    else if (dp.ShowMode == 2)
                    {
                        GL.Color3(SelectedColor);
                    }
                    else
                    {
                        GL.Color3(0.0, 0.0, 0.0);
                    }

                    GL.Translate(0.0, 0.0, height);
                    dp.DrawElements();

                    GL.Translate(0.0, 0.0, -height);
                    if (isTexture)
                    {
                        GL.Enable(EnableCap.Texture2D);
                    }
                }
                if (dp.Type == CadElementType.Edge)
                {
                    GL.Disable(EnableCap.Texture2D);
                    GL.LineWidth(LineWidth);
                    if (IsAntiAliasing)
                    {
                        GL.Enable(EnableCap.LineSmooth);
                        GL.Enable(EnableCap.Blend);
                        GL.BlendFunc(BlendingFactor.SrcAlpha, BlendingFactor.OneMinusSrcAlpha);
                        GL.Hint(HintTarget.LineSmoothHint, HintMode.DontCare);
                    }
                    if (dp.ShowMode == 1)
                    {
                        GL.Color3(1.0, 1.0, 0.0);
                    }
                    else if (dp.ShowMode == 2)
                    {
                        GL.Color3(SelectedColor);
                    }
                    else
                    {
                        //2019-03-11 エッジの色 FIX
                        //GL.Color3(0, 0, 0);
                        GL.Color3(dp.Color);
                    }

                    GL.Translate(0.0, 0.0, height);
                    dp.DrawElements();
                    if (dp.ShowMode > 0)
                    {
                        // draw ctrl point        
                        GL.Begin(PrimitiveType.Points);
                        for (int icp = 0; icp < dp.CtrlPoints.Count; icp++)
                        {
                            OpenTK.Vector2d cp = dp.CtrlPoints[icp];
                            GL.Vertex3(cp.X, cp.Y, 0.0);
                        }
                        GL.End();
                        // draw line between ctrl point and point
                        GL.Enable(EnableCap.LineStipple);
                        GL.LineStipple(1, 0xF0F0);
                        GL.PolygonStipple(Mask);
                        GL.LineWidth(1);
                        GL.Begin(PrimitiveType.Lines);
                        uint sVI = dp.Indexs[0];
                        uint eVI = dp.Indexs[dp.ElemCount * dp.ElemPtCount - 1];
                        double[] va = VertexArray.VertexCoords;
                        OpenTK.Vector2d sPt = new OpenTK.Vector2d(va[sVI * 2 + 0], va[sVI * 2 + 1]);
                        OpenTK.Vector2d ePt = new OpenTK.Vector2d(va[eVI * 2 + 0], va[eVI * 2 + 1]);
                        if (dp.CurveType == CurveType.CurveArc)
                        {
                            OpenGLUtils.GLVertex2(sPt);
                            OpenGLUtils.GLVertex2(dp.CtrlPoints[0]);
                            OpenGLUtils.GLVertex2(ePt);
                            OpenGLUtils.GLVertex2(dp.CtrlPoints[0]);
                        }
                        if (dp.CurveType == CurveType.CurveBezier)
                        {
                            OpenGLUtils.GLVertex2(sPt);
                            OpenGLUtils.GLVertex2(dp.CtrlPoints[0]);
                            OpenGLUtils.GLVertex2(ePt);
                            OpenGLUtils.GLVertex2(dp.CtrlPoints[1]);
                        }
                        GL.End();
                        GL.Disable(EnableCap.LineStipple);
                    }
                    GL.Translate(0.0, 0.0, -height);
                    GL.Disable(EnableCap.LineSmooth);
                    GL.Disable(EnableCap.Blend);
                    if (isTexture)
                    {
                        GL.Enable(EnableCap.Texture2D);
                    }
                }
                else if (dp.Type == CadElementType.Loop)
                {
                    GL.Disable(EnableCap.Blend);
                    if (dp.ShowMode > 0)
                    {
                        GL.Enable(EnableCap.PolygonStipple);
                        GL.PolygonStipple(Mask);
                        if (dp.ShowMode == 1)
                        {
                            GL.Color3(1.0, 1.0, 0.0);
                        }
                        else if (dp.ShowMode == 2)
                        {
                            GL.Color3(SelectedColor);
                        }
                        GL.Translate(0.0, 0.0, +height + 0.001);
                        dp.DrawElements();
                        GL.Translate(0.0, 0.0, -height - 0.001);
                        GL.Disable(EnableCap.PolygonStipple);
                    }
                    if (dp.ShowMode != 0)
                    {
                        continue;
                    }
                    float[] color1 = { (float)dp.Color[0], (float)dp.Color[1], (float)dp.Color[2] };
                    GL.Material(MaterialFace.FrontAndBack, MaterialParameter.Diffuse, color1);
                    GL.Color3(dp.Color);
                    GL.Translate(+dispX, +dispY, +height);
                    dp.DrawElements();
                    GL.Translate(-dispX, -dispY, -height);
                    if (isBlend)
                    {
                        GL.Enable(EnableCap.Blend);
                    }
                }
            }
            GL.DisableClientState(ArrayCap.VertexArray);
            GL.DisableClientState(ArrayCap.TextureCoordArray);

            if (isLighting)
            {
                GL.Enable(EnableCap.Lighting);
            }
            else
            {
                GL.Disable(EnableCap.Lighting);
            }
            if (isBlend)
            {
                GL.Enable(EnableCap.Blend);
            }
            else
            {
                GL.Disable(EnableCap.Blend);
            }
            if (isTexture)
            {
                GL.Enable(EnableCap.Texture2D);
            }
            else
            {
                GL.Disable(EnableCap.Texture2D);
            }
        }

        public void DrawSelection(uint iDraw)
        {
            bool isBlend = GL.IsEnabled(EnableCap.Blend);
            bool isLineSmooth = GL.IsEnabled(EnableCap.LineSmooth);
            bool isTexture = GL.IsEnabled(EnableCap.Texture2D);
            GL.Disable(EnableCap.Blend);
            GL.Disable(EnableCap.LineSmooth);
            GL.Disable(EnableCap.Texture2D);

            uint ndim = VertexArray.Dimension;
            GL.PushName(iDraw);
            // モデルの描画
            GL.EnableClientState(ArrayCap.VertexArray);
            GL.VertexPointer((int)ndim, VertexPointerType.Double, 0, VertexArray.VertexCoords);
            for (int idp = 0; idp < DrawParts.Count; idp++)
            {
                Cad2DDrawPart dp = DrawParts[idp];
                double height = dp.Height;
                GL.PushName(idp);
                GL.Translate(0.0, 0.0, +height);
                dp.DrawElements();
                if (dp.Type == CadElementType.Edge && dp.ShowMode == 2)
                {
                    for (int icp = 0; icp < dp.CtrlPoints.Count; icp++)
                    {
                        OpenTK.Vector2d cp = dp.CtrlPoints[icp];
                        GL.PushName(icp);
                        GL.Begin(PrimitiveType.Points);
                        OpenGLUtils.GLVertex2(cp);
                        GL.End();
                        GL.PopName();
                    }
                }
                GL.Translate(0.0, 0.0, -height);
                GL.PopName();
            }
            GL.DisableClientState(ArrayCap.VertexArray);
            GL.PopName();

            if (isBlend)
            {
                GL.Enable(EnableCap.Blend);
            }
            else
            {
                GL.Disable(EnableCap.Blend);
            }
            if (isTexture)
            {
                GL.Enable(EnableCap.Texture2D);
            }
            else
            {
                GL.Disable(EnableCap.Texture2D);
            }
        }

        public void GetPartCadId(int[] selectFlag, 
            out CadElementType partType, out uint partId, out int ctrlIndex)
        {
            uint idp = (uint)selectFlag[1];
            if (idp < DrawParts.Count)
            {
                Cad2DDrawPart dp = DrawParts[(int)idp];
                partType = dp.Type;
                partId = dp.CadId;
                ctrlIndex = selectFlag[2];
                return;
            }
            partType = CadElementType.NotSet;
            partId = 0;
            ctrlIndex = 0;
        }

        public void SetShowMode(CadElementType type, uint id, int showMode)
        {
            bool flag = false;
            for (int idp = 0; idp < DrawParts.Count; idp++)
            {
                if (DrawParts[idp].Type == type && DrawParts[idp].CadId == id)
                {
                    DrawParts[idp].ShowMode = showMode;
                    System.Diagnostics.Debug.Assert(flag == false);
                    flag = true;
                }
            }
        }

        public int GetShowMode(CadElementType type, uint id)
        {
            for (int idp = 0; idp < DrawParts.Count; idp++)
            {
                if (DrawParts[idp].Type == type && DrawParts[idp].CadId == id)
                {
                    return DrawParts[idp].ShowMode;
                }
            }
            return 0;
        }

        public void ClearShowMode(int curShowMode)
        {
            for (int idp = 0; idp < DrawParts.Count; idp++)
            {
                if (DrawParts[idp].ShowMode == curShowMode)
                {
                    DrawParts[idp].ShowMode = 0;
                }
            }
        }

        public void SetIsShow(bool isShow, CadElementType cadPartType, uint cadPartId)
        {
            int mode = (isShow) ? 0 : -1;
            for (int idp = 0; idp < DrawParts.Count; idp++)
            {
                if (DrawParts[idp].CadId == cadPartId &&
                    DrawParts[idp].Type == cadPartType)
                {
                    DrawParts[idp].ShowMode = mode;
                }
            }
        }

        public void SetIsShow(bool isShow, CadElementType cadPartType, IList<uint> partIds)
        {
            for (int i = 0; i < partIds.Count; i++)
            {
                uint id = partIds[i];
                SetIsShow(isShow, cadPartType, id);
            }
        }

        public void AddSelected(int[] selectFlag)
        {
            CadElementType partType;
            uint partId;
            int ctrlIndex;
            GetPartCadId(selectFlag, out partType, out partId, out ctrlIndex);
            if (partType == CadElementType.NotSet)
            {
                return;
            }
            int showMode = 2;
            SetShowMode(partType, partId, showMode);
        }

        public void ClearSelected()
        {
            int curShowMode = 2;
            ClearShowMode(curShowMode);
        }
    }
}
