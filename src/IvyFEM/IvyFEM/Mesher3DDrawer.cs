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

        private IList<Mesher3DDrawPart> DrawParts = new List<Mesher3DDrawPart>();
        private VertexArray VertexArray = new VertexArray();

        public Mesher3DDrawer()
        {

        }

        public Mesher3DDrawer(Mesher3D mesher)
        {
            Set(mesher);
        }

        private bool Set(Mesher3D mesher)
        {
            SutableRotMode = RotMode.RotMode3D;

            {
                // 三角形要素をセット
                IList<MeshTriArray3D> triArrays = mesher.GetTriArrays();
                for (int itri = 0; itri < triArrays.Count; itri++)
                {
                    Mesher3DDrawPart dp = new Mesher3DDrawPart(triArrays[itri]);
                    DrawParts.Add(dp);
                }
            }

            {
                // 線要素をセット
                IList<MeshBarArray3D> barArrays = mesher.GetBarArrays();
                for (int ibar = 0; ibar < barArrays.Count; ibar++)
                {
                    Mesher3DDrawPart dp = new Mesher3DDrawPart(barArrays[ibar]);
                    DrawParts.Add(dp);
                }
            }

            {
                // 頂点をセット
                IList<MeshVertex3D> vertexs = mesher.GetVertexs();
                for (int iver = 0; iver < vertexs.Count; iver++)
                {
                    Mesher3DDrawPart dp = new Mesher3DDrawPart(vertexs[iver]);
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
