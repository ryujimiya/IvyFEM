using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    public class Mesher2DDrawer : IDrawer
    {
        public RotMode SutableRotMode { get; private set; } = RotMode.RotMode2D;
        public bool IsAntiAliasing { get; set; } = false;

        private IList<Mesher2DDrawPart> DrawParts = new List<Mesher2DDrawPart>();
        private VertexArray VertexArray = new VertexArray();
        private bool IsDrawFace;

        public Mesher2DDrawer()
        {

        }

        public Mesher2DDrawer(Mesher2D mesher, bool isDrawFace = true)
        {
            IsDrawFace = isDrawFace;
            Set(mesher);
        }

        private bool Set(Mesher2D mesher)
        {
            SutableRotMode = RotMode.RotMode2D;

            Cad2D cad = mesher.Cad;
            int layerMin = 0;
            int layerMax = 0;
            {
                bool isInited = false;
                IList<MeshTriArray2D> triArrays = mesher.GetTriArrays();
                for (int itri = 0; itri < triArrays.Count; itri++)
                {
                    int layer = triArrays[itri].Layer;
                    if (isInited)
                    {
                        layerMin = (layer < layerMin) ? layer : layerMin;
                        layerMax = (layer > layerMax) ? layer : layerMax;
                    }
                    else
                    {
                        layerMin = layer;
                        layerMax = layer;
                        isInited = true;
                    }
                }
            }
            double layerHeight = 1.0 / (layerMax - layerMin + 1);

            {
                // 三角形要素をセット
                IList<MeshTriArray2D> triArrays = mesher.GetTriArrays();

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

                    Mesher2DDrawPart dp = new Mesher2DDrawPart(triArrays[loc]);
                    int layer = triArrays[loc].Layer;
                    dp.Height = (layer - layerMin) * layerHeight;
                    DrawParts.Add(dp);
                }
            }

            {
                // 線要素をセット
                IList<MeshBarArray2D> barArrays = mesher.GetBarArrays();

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

                    Mesher2DDrawPart dp = new Mesher2DDrawPart(barArrays[loc]);
                    double height = 0;
                    {
                        int layer = barArrays[loc].Layer;
                        height += (layer - layerMin + 0.01) * layerHeight;
                    }
                    dp.Height = height;
                    DrawParts.Add(dp);
                }
            }

            {
                // 頂点をセット
                IList<MeshVertex2D> vertexs = mesher.GetVertexs();
                for (int iver = 0; iver < vertexs.Count; iver++)
                {
                    double height = 0;
                    /*
                    {
                        int layer = vertexs[iver].Layer;
                        height += (layer - layerMin + 0.1) * layerHeight;
                    }
                    */
                    height = 0.2;
                    Mesher2DDrawPart dp = new Mesher2DDrawPart(vertexs[iver]);
                    dp.Height = height;

                    DrawParts.Add(dp);
                }
            }

            {
                // 座標をセット
                IList<OpenTK.Vector2d> vecs = mesher.GetVectors();
                uint nDim = 2;
                uint nVec = (uint)vecs.Count;
                VertexArray.SetSize(nVec, nDim);
                for (int ivec = 0; ivec < nVec; ivec++)
                {
                    VertexArray.VertexCoords[ivec * nDim] = vecs[ivec].X;
                    VertexArray.VertexCoords[ivec * nDim + 1] = vecs[ivec].Y;
                }
            }
            return true;
        }

        public void UpdateCoord(Mesher2D mesher)
        {
            {
                // 座標をセット
                IList<OpenTK.Vector2d> vecs = mesher.GetVectors();
                uint nDim = 2;
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
                }
            }
        }

        public BoundingBox3D GetBoundingBox(OpenTK.Matrix3d rot)
        {
            return VertexArray.GetBoundingBox(rot);
        }

        public void Draw()
        {
            // ライティングの指定
            GL.Disable(EnableCap.Lighting);
            // 色の指定
            GL.Color3(0.8, 0.8, 0.8);

            // 片面かどうかの指定
            GL.Enable(EnableCap.CullFace);

            GL.CullFace(CullFaceMode.Back);
            //GL.Disable(EnableCap.CullFace);
            GL.Enable(EnableCap.DepthTest);
            //GL.Disable(EnableCap.DepthTest);

            uint nDim = VertexArray.Dimension;

            // 頂点配列の設定
            GL.EnableClientState(ArrayCap.VertexArray);
            GL.VertexPointer((int)nDim, VertexPointerType.Double, 0, VertexArray.VertexCoords);
            
            for (int idp = 0; idp < DrawParts.Count; idp++)
            {
                Mesher2DDrawPart dp = DrawParts[idp];
                if (!IsDrawFace && (dp.Dimension == 2))
                {
                    continue;
                }
                double height = dp.Height;

                GL.Translate(0, 0, +height);
                dp.DrawElements();
                GL.Translate(0, 0, -height);
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
                Mesher2DDrawPart dp = DrawParts[idp];
                double height = dp.Height;

                GL.PushName(idp);
                GL.Translate(0, 0, +height);
                dp.DrawElementsSelection();
                GL.Translate(0, 0, -height);
                GL.PopName();
            }
            GL.PopName();

            GL.DisableClientState(ArrayCap.VertexArray);
        }

        public void AddSelected(int[] selectFlag)
        {
            int idp0 = selectFlag[1];
            int ielem0 = selectFlag[2];
            System.Diagnostics.Debug.WriteLine("Mesher2DDrawer.AddSelected idp0 = " +
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
                Mesher2DDrawPart dp = DrawParts[idp];
                dp.IsSelected = false;
                dp.SelectedElems.Clear();
            }
        }
    }
}
