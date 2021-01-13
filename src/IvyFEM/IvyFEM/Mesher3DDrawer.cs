using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    public class Mesher3DDrawer : IDrawer
    {
        public RotMode SutableRotMode { get; private set; } = RotMode.RotMode3D;
        public bool IsAntiAliasing { get; set; } = false;
        public bool IsMask { get; set; } = false;
        private byte[] Mask = new byte[128];
        private IList<Mesher3DDrawPart> DrawParts = new List<Mesher3DDrawPart>();
        private VertexArray VertexArray = new VertexArray();

        public Mesher3DDrawer()
        {
            SetupMask();
        }

        public Mesher3DDrawer(Mesher3D mesher)
        {
            SetupMask();
            Set(mesher);
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

        private bool Set(Mesher3D mesher)
        {
            SutableRotMode = RotMode.RotMode3D;

            Cad3D cad = mesher.Cad;
            {
                // 四面体要素をセット
                IList<MeshTetArray> tetArrays = mesher.GetTetArrays();

                IList<uint> sIds = cad.GetElementIds(CadElementType.Solid);
                for (int iSId = 0; iSId < sIds.Count; iSId++)
                {
                    uint sId = sIds[iSId];

                    uint cadId = sId;
                    CadElementType cadType = CadElementType.Solid;
                    if (!cad.IsElementId(cadType, cadId))
                    {
                        continue;
                    }
                    uint meshId = mesher.GetIdFromCadId(cadId, cadType);
                    if (meshId == 0)
                    {
                        continue;
                    }
                    MeshType meshType;
                    uint elemCnt;
                    int loc;
                    uint cadId0;
                    mesher.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadId0);
                    System.Diagnostics.Debug.Assert(cadId0 == cadId);
                    if (meshType != MeshType.Tet)
                    {
                        continue;
                    }

                    Mesher3DDrawPart dp = new Mesher3DDrawPart(tetArrays[loc]);
                    DrawParts.Add(dp);
                }
            }
            {
                // 三角形要素をセット
                IList<MeshTriArray3D> triArrays = mesher.GetTriArrays();

                IList<uint> lIds = cad.GetElementIds(CadElementType.Loop);
                for (int iLId = 0; iLId < lIds.Count; iLId++)
                {
                    uint lId = lIds[iLId];

                    uint cadId = lId;
                    CadElementType cadType = CadElementType.Loop;
                    if (!cad.IsElementId(cadType, cadId))
                    {
                        continue;
                    }
                    uint meshId = mesher.GetIdFromCadId(cadId, cadType);
                    if (meshId == 0)
                    {
                        continue;
                    }
                    MeshType meshType;
                    uint elemCnt;
                    int loc;
                    uint cadId0;
                    mesher.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadId0);
                    System.Diagnostics.Debug.Assert(cadId0 == cadId);
                    if (meshType != MeshType.Tri)
                    {
                        continue;
                    }

                    Mesher3DDrawPart dp = new Mesher3DDrawPart(triArrays[loc]);
                    DrawParts.Add(dp);
                }
            }

            {
                // 線要素をセット
                IList<MeshBarArray3D> barArrays = mesher.GetBarArrays();

                IList<uint> eIds = cad.GetElementIds(CadElementType.Edge);
                for (int iEId = 0; iEId < eIds.Count; iEId++)
                {
                    uint eId = eIds[iEId];

                    uint cadId = eId;
                    CadElementType cadType = CadElementType.Edge;
                    if (!cad.IsElementId(cadType, cadId))
                    {
                        continue;
                    }
                    uint meshId = mesher.GetIdFromCadId(cadId, cadType);
                    if (meshId == 0)
                    {
                        continue;
                    }
                    MeshType meshType;
                    uint elemCnt;
                    int loc;
                    uint cadId0;
                    mesher.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadId0);
                    System.Diagnostics.Debug.Assert(cadId0 == cadId);
                    if (meshType != MeshType.Bar)
                    {
                        continue;
                    }

                    Mesher3DDrawPart dp = new Mesher3DDrawPart(barArrays[loc]);
                    DrawParts.Add(dp);
                }
            }

            {
                // 頂点をセット
                IList<MeshVertex3D> vertexs = mesher.GetVertexs();

                IList<uint> vIds = cad.GetElementIds(CadElementType.Vertex);
                for (int iVId = 0; iVId < vIds.Count; iVId++)
                {
                    uint vCadId = vIds[iVId];

                    uint cadId = vCadId;
                    CadElementType cadType = CadElementType.Vertex;
                    if (!cad.IsElementId(cadType, cadId))
                    {
                        continue;
                    }
                    uint meshId = mesher.GetIdFromCadId(cadId, cadType);
                    if (meshId == 0)
                    {
                        continue;
                    }
                    MeshType meshType;
                    uint elemCnt;
                    int loc;
                    uint cadId0;
                    mesher.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadId0);
                    System.Diagnostics.Debug.Assert(cadId0 == cadId);
                    if (meshType != MeshType.Vertex)
                    {
                        continue;
                    }

                    Mesher3DDrawPart dp = new Mesher3DDrawPart(vertexs[loc]);
                    DrawParts.Add(dp);
                }
            }

            {
                // 座標をセット
                IList<OpenTK.Vector3d> vecs = mesher.GetVectors();
                uint nDim = 3;
                uint nVec = (uint)vecs.Count;
                VertexArray.SetSize(nVec, nDim);
                for (int ivec = 0; ivec < nVec; ivec++)
                {
                    VertexArray.VertexCoords[ivec * nDim] = vecs[ivec].X;
                    VertexArray.VertexCoords[ivec * nDim + 1] = vecs[ivec].Y;
                    VertexArray.VertexCoords[ivec * nDim + 2] = vecs[ivec].Z;
                }
            }
            return true;
        }

        public void UpdateCoord(Mesher3D mesher)
        {
            {
                // 座標をセット
                IList<OpenTK.Vector3d> vecs = mesher.GetVectors();
                uint nDim = 3;
                uint nVec = (uint)vecs.Count;
                if (VertexArray.Dimension != nDim)
                {
                    return;
                }
                if (VertexArray.PointCount != nVec)
                {
                    return;
                }
                for (int ivec = 0; ivec < nVec; ivec++)
                {
                    VertexArray.VertexCoords[ivec * nDim] = vecs[ivec].X;
                    VertexArray.VertexCoords[ivec * nDim + 1] = vecs[ivec].Y;
                    VertexArray.VertexCoords[ivec * nDim + 2] = vecs[ivec].Z;
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
            // ライティングの指定
            GL.Disable(EnableCap.Lighting);
            // 色の指定
            //GL.Color3(0.8, 0.8, 0.8);

            uint nDim = VertexArray.Dimension;

            // 頂点配列の設定
            GL.EnableClientState(ArrayCap.VertexArray);
            GL.VertexPointer((int)nDim, VertexPointerType.Double, 0, VertexArray.VertexCoords);

            GL.Disable(EnableCap.CullFace);

            for (int idp = 0; idp < DrawParts.Count; idp++)
            {
                Mesher3DDrawPart dp = DrawParts[idp];
                dp.Mask = IsMask ? Mask : null;
                dp.DrawElements();
            }

            GL.DisableClientState(ArrayCap.VertexArray);
        }

        public void DrawSelection(uint iDraw)
        {
            uint nDim = VertexArray.Dimension;

            // モデルの描画
            GL.EnableClientState(ArrayCap.VertexArray);

            GL.VertexPointer((int)nDim, VertexPointerType.Double, 0, VertexArray.VertexCoords);

            GL.PushName(iDraw);
            for (int idp = 0; idp < DrawParts.Count; idp++)
            {
                Mesher3DDrawPart dp = DrawParts[idp];

                GL.PushName(idp);
                dp.DrawElementsSelection();
                GL.PopName();
            }
            GL.PopName();

            GL.DisableClientState(ArrayCap.VertexArray);
        }

        public void AddSelected(int[] selectFlag)
        {
            int idp0 = selectFlag[1];
            int ielem0 = selectFlag[2];
            System.Diagnostics.Debug.WriteLine("Mesher3DDrawer.AddSelected idp0 = " +
                idp0 + " ielem0 = " + ielem0);
            if (idp0 >= 0 && idp0 < DrawParts.Count)
            {
                DrawParts[idp0].IsSelected = true;
                IList<uint> selectedElems = DrawParts[idp0].SelectedElems;
                if (selectedElems.IndexOf((uint)ielem0) < 0)
                {
                    selectedElems.Add((uint)ielem0);
                }
            }
        }

        public void ClearSelected()
        {
            for (int idp = 0; idp < DrawParts.Count; idp++)
            {
                Mesher3DDrawPart dp = DrawParts[idp];
                dp.IsSelected = false;
                dp.SelectedElems.Clear();
            }
        }
    }
}
