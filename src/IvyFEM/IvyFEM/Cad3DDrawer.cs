using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    public class Cad3DDrawer : IDrawer
    {
        public RotMode SutableRotMode { get; private set; } = RotMode.RotMode3D;
        public bool IsAntiAliasing { get; set; } = false;

        private double[] SelectedColor = { 1.0, 0.5, 1.0 }; //{ 1.0, 1.0, 0.0 };
        public bool IsMask { get; set; } = false;
        private byte[] Mask = new byte[128];
        private IList<Cad3DDrawPart> DrawParts = new List<Cad3DDrawPart>();
        private VertexArray VertexArray = new VertexArray();
        public uint LineWidth { get; set; } = 3;
        private uint PointSize = 5;

        public Cad3DDrawer()
        {
            SetupMask();
        }

        public Cad3DDrawer(Cad3D cad)
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

        public bool UpdateCadTopologyGeometry(Cad3D cad)
        {
            SutableRotMode = RotMode.RotMode3D;
            IList<Cad3DDrawPart> oldDrawParts = new List<Cad3DDrawPart>();
            for (int i = 0; i < DrawParts.Count; i++)
            {
                oldDrawParts.Add(new Cad3DDrawPart(DrawParts[i]));
            }

            for (int idp = 0; idp < oldDrawParts.Count; idp++)
            {
                oldDrawParts[idp].MeshId = 0;
                oldDrawParts[idp].ShowMode = 0;
            }
            DrawParts.Clear();

            {
                // 面をセット
                IList<uint> lIds = cad.GetElementIds(CadElementType.Loop);
                for (int iLId = 0; iLId < lIds.Count; iLId++)
                {
                    uint lId = lIds[iLId];
                    int idp0 = 0;
                    for (; idp0 < oldDrawParts.Count; idp0++)
                    {
                        Cad3DDrawPart olddp = oldDrawParts[idp0];
                        if (olddp.Type == CadElementType.Loop && olddp.CadId == lId)
                        {
                            olddp.MeshId = 1;
                            DrawParts.Add(oldDrawParts[idp0]);
                            break;
                        }
                    }
                    if (idp0 == oldDrawParts.Count)
                    {
                        Cad3DDrawPart dp = new Cad3DDrawPart();
                        dp.CadId = lId;
                        dp.Type = CadElementType.Loop;
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
                    int idp0 = 0;
                    for (; idp0 < oldDrawParts.Count; idp0++)
                    {
                        Cad3DDrawPart olddp = oldDrawParts[idp0];
                        if (olddp.Type == CadElementType.Edge && olddp.CadId == eId)
                        {
                            olddp.MeshId = 1;
                            DrawParts.Add(olddp);
                            break;
                        }
                    }
                    if (idp0 == oldDrawParts.Count)
                    {
                        Cad3DDrawPart dp = new Cad3DDrawPart();
                        dp.CadId = eId;
                        dp.Type = CadElementType.Edge;
                        DrawParts.Add(dp);
                    }
                }
            }

            { 
                // set vertex
                IList<uint> vIds = cad.GetElementIds(CadElementType.Vertex);
                for (int iVId = 0; iVId < vIds.Count; iVId++)
                {
                    uint vCadId = vIds[iVId];
                    int idp0 = 0;
                    for (; idp0 < oldDrawParts.Count; idp0++)
                    {
                        Cad3DDrawPart olddp = oldDrawParts[idp0];
                        if (olddp.Type == CadElementType.Vertex && olddp.CadId == vCadId)
                        {
                            olddp.MeshId = 1;
                            DrawParts.Add(olddp);
                            break;
                        }
                    }
                    if (idp0 == oldDrawParts.Count)
                    {
                        Cad3DDrawPart dp = new Cad3DDrawPart();
                        dp.CadId = vCadId;
                        dp.Type = CadElementType.Vertex;
                        DrawParts.Add(dp);
                    }
                }
            }

            oldDrawParts.Clear();

            UpdateCadGeometry(cad);
            return true;
        }

        public void UpdateCadGeometry(Cad3D cad)
        {
            Mesher3D mesh = new Mesher3D(cad);
            for (int idp = 0; idp < DrawParts.Count; idp++)
            {
                Cad3DDrawPart dp = DrawParts[idp];
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
                    dp.Color = new double[color.Length];
                    color.CopyTo(dp.Color, 0);
                }
                else if (meshType == MeshType.Bar)
                {
                    dp.SetBarArray(mesh.GetBarArrays()[loc]);

                    System.Diagnostics.Debug.Assert(cadType == CadElementType.Edge);
                    Edge3D edge = cad.GetEdge(cadId);
                    double[] color = edge.Color;
                    color.CopyTo(dp.Color, 0);
                }
                else if (meshType == MeshType.Vertex)
                {
                    dp.SetVertex(mesh.GetVertexs()[loc]);

                    System.Diagnostics.Debug.Assert(cadType == CadElementType.Vertex);
                    Vertex3D v = cad.GetVertex(cadId);
                    double[] color = v.Color;
                    color.CopyTo(dp.Color, 0);
                }
            }

            {
                // 座標をセット
                IList<OpenTK.Vector3d> vecs = mesh.GetVectors();
                uint ptCnt = (uint)vecs.Count;
                uint ndim = 3;
                VertexArray.SetSize(ptCnt, ndim);
                for (int iPt = 0; iPt < ptCnt; iPt++)
                {
                    VertexArray.VertexCoords[iPt * ndim] = vecs[iPt].X;
                    VertexArray.VertexCoords[iPt * ndim + 1] = vecs[iPt].Y;
                    VertexArray.VertexCoords[iPt * ndim + 2] = vecs[iPt].Z;
                }
                if (VertexArray.UVCoords != null)
                {
                    // Not implemented
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

            //// alpha blend
            //GL.Enable(EnableCap.Blend);
            //GL.BlendFunc(BlendingFactor.SrcAlpha, BlendingFactor.OneMinusSrcAlpha);

            uint ndim = VertexArray.Dimension;

            ////////////////////////////////////////////////////////////////
            // モデルの描画

            /////////////
            // vertex arrayを登録する
            GL.EnableClientState(ArrayCap.VertexArray);
            GL.VertexPointer((int)ndim, VertexPointerType.Double, 0, VertexArray.VertexCoords);
            GL.PointSize(PointSize);
            GL.LineWidth(LineWidth);
            for (int idp = 0; idp < DrawParts.Count; idp++)
            {
                Cad3DDrawPart dp = DrawParts[idp];
                if (dp.ShowMode == -1)
                {
                    continue;
                }
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
                        //GL.Color3(0.0, 0.0, 0.0);
                        GL.Color3(dp.Color);
                    }

                    dp.DrawElements();
                    if (isTexture)
                    {
                        GL.Enable(EnableCap.Texture2D);
                    }
                }
                if (dp.Type == CadElementType.Edge)
                {
                    GL.Disable(EnableCap.Lighting);
                    GL.Disable(EnableCap.Texture2D);

                    GL.LineWidth(LineWidth);
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
                        GL.Color3(dp.Color);
                    }
                    dp.DrawElements();

                    if (isLighting)
                    {
                        GL.Enable(EnableCap.Lighting);
                    }
                    if (isTexture)
                    {
                        GL.Enable(EnableCap.Texture2D);
                    }
                }
                else if (dp.Type == CadElementType.Loop)
                {
                    if (GL.IsEnabled(EnableCap.Lighting))
                    {
                        GL.Normal3(dp.Normal[0], dp.Normal[1], dp.Normal[2]);
                    }
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

                        dp.DrawElements();
                        GL.Disable(EnableCap.PolygonStipple);
                    }
                    if (dp.ShowMode != 0)
                    {
                        continue;
                    }
                    if (IsMask)
                    {
                        GL.Enable(EnableCap.PolygonStipple);
                        GL.PolygonStipple(Mask);
                    }
                    _GLColor(dp.Color);
                    dp.DrawElements();
                    if (IsMask)
                    {
                        GL.Disable(EnableCap.PolygonStipple);
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

        private void _GLColor(double[] color)
        {
            if (color.Length == 3)
            {
                GL.Color3(color);
            }
            else if (color.Length == 4)
            {
                GL.Color4(color);
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
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
                Cad3DDrawPart dp = DrawParts[idp];
                GL.PushName(idp);
                dp.DrawElements();
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
                Cad3DDrawPart dp = DrawParts[(int)idp];
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
