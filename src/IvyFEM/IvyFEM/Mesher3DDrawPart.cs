using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    public class Mesher3DDrawPart
    {
        public byte[] Mask { get; set; } = null;
        public bool IsSelected { get; set; } = false;
        public bool IsShown { get; set; } = true;
        public IList<uint> SelectedElems { get; set; } = new List<uint>();
        public uint MeshId { get; set; } = 0;
        public uint CadId { get; set; } = 0;
        public double[] Color { get; set;  } = new double[3] { 0.8, 0.8, 0.8 };
        public uint LineWidth { get; set; } = 1;
        public uint ElemCount { get; set; } = 0;
        public int[] ElemIndexs { get; set; } = null;
        public uint EdgeCount { get; set; } = 0;
        public int[] EdgeIndexs { get; set; } = null;
        public uint Dimension
        {
            get
            {
                if (Type == ElementType.Point)
                {
                    return 1;
                }
                else if (Type == ElementType.Line)
                {
                    return 1;
                }
                else if (Type == ElementType.Tri)
                {
                    return 2;
                }
                else if (Type == ElementType.Tet)
                {
                    return 3;
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
                return 0;
            }
        }

        public double Height { get; set; } = 0.0;

        private ElementType Type;

        public Mesher3DDrawPart()
        {

        }

        public Mesher3DDrawPart(Mesher3DDrawPart src)
        {
            IsSelected = src.IsSelected;
            IsShown = src.IsShown;
            SelectedElems = new List<uint>(src.SelectedElems);
            MeshId = src.MeshId;
            CadId = src.CadId;
            src.Color.CopyTo(Color, 0);
            LineWidth = src.LineWidth;
            ElemCount = src.ElemCount;
            ElemIndexs = null;
            if (src.ElemIndexs != null)
            {
                ElemIndexs = new int[src.EdgeIndexs.Length];
                src.ElemIndexs.CopyTo(ElemIndexs, 0);
            }
            EdgeCount = src.EdgeCount;
            EdgeIndexs = null;
            if (src.EdgeIndexs != null)
            {
                EdgeIndexs = new int[src.EdgeIndexs.Length];
                src.EdgeIndexs.CopyTo(EdgeIndexs, 0);
            }
            Type = src.Type;
        }

        public Mesher3DDrawPart(MeshTetArray tetArray)
        {
            IsSelected = false;
            IsShown = true;
            Color = new double[3] { 0.8, 0.8, 0.8 };
            LineWidth = 1;
            ElemIndexs = null;
            EdgeIndexs = null;

            CadId = tetArray.SCadId;

            MeshId = tetArray.Id;
            System.Diagnostics.Debug.Assert(MeshId != 0);
            Type = ElementType.Tet;

            Height = 0;

            ////////////////////////////////
            // 面のセット
            ElemCount = 0;
            for (uint iTet = 0; iTet < tetArray.Tets.Count; iTet++)
            {
                for (uint ielemtet = 0; ielemtet < 4; ielemtet++)
                {
                    if (tetArray.Tets[(int)iTet].G[ielemtet] == -2)
                    {
                        continue;
                    }
                    ElemCount++;
                }
            }
            ElemIndexs = new int[ElemCount * 3];
            {
                uint icnt = 0;
                for (uint iTet = 0; iTet < tetArray.Tets.Count; iTet++)
                {
                    for (uint iTetFace = 0; iTetFace < 4; iTetFace++)
                    {
                        if (tetArray.Tets[(int)iTet].G[iTetFace] == -2)
                        {
                            continue;
                        }
                        ElemIndexs[icnt * 3 + 0] =
                            (int)tetArray.Tets[(int)iTet].V[(int)MeshUtils3D.ElNoTetFace[iTetFace][0]];
                        ElemIndexs[icnt * 3 + 1] =
                            (int)tetArray.Tets[(int)iTet].V[(int)MeshUtils3D.ElNoTetFace[iTetFace][1]];
                        ElemIndexs[icnt * 3 + 2] =
                            (int)tetArray.Tets[(int)iTet].V[(int)MeshUtils3D.ElNoTetFace[iTetFace][2]];
                        icnt++;
                    }
                }
            }

            ////////////////////////////////
            // 辺のセット
            EdgeCount = ElemCount * 3;
            EdgeIndexs = new int[EdgeCount * 2];
            for (uint ielem = 0; ielem < ElemCount; ielem++)
            {
                EdgeIndexs[(ielem * 3) * 2 + 0] = ElemIndexs[ielem * 3];
                EdgeIndexs[(ielem * 3) * 2 + 1] = ElemIndexs[ielem * 3 + 1];

                EdgeIndexs[(ielem * 3 + 1) * 2 + 0] = ElemIndexs[ielem * 3 + 1];
                EdgeIndexs[(ielem * 3 + 1) * 2 + 1] = ElemIndexs[ielem * 3 + 2];

                EdgeIndexs[(ielem * 3 + 2) * 2 + 0] = ElemIndexs[ielem * 3 + 2];
                EdgeIndexs[(ielem * 3 + 2) * 2 + 1] = ElemIndexs[ielem * 3];
            }
        }

        public Mesher3DDrawPart(MeshTriArray3D triArray)
        {
            IsSelected = false;
            IsShown = true;
            Color = new double[3] { 0.8, 0.8, 0.8 };
            LineWidth = 1;
            ElemIndexs = null;
            EdgeIndexs = null;

            CadId = triArray.LCadId;

            MeshId = triArray.Id;
            System.Diagnostics.Debug.Assert(MeshId != 0);
            Type = ElementType.Tri;

            ElemCount = (uint)triArray.Tris.Count;

            {
                // 面のセット
                ElemIndexs = new int[ElemCount * 3];
                for (int iTri = 0; iTri < ElemCount; iTri++)
                {
                    ElemIndexs[iTri * 3 + 0] = (int)triArray.Tris[iTri].V[0];
                    ElemIndexs[iTri * 3 + 1] = (int)triArray.Tris[iTri].V[1];
                    ElemIndexs[iTri * 3 + 2] = (int)triArray.Tris[iTri].V[2];
                }
            }

            {
                // 辺のセット
                EdgeCount = ElemCount * 3;
                EdgeIndexs = new int[EdgeCount * 2];
                for (int iTri = 0; iTri < ElemCount; iTri++)
                {
                    EdgeIndexs[(iTri * 3) * 2 + 0] = (int)triArray.Tris[iTri].V[0];
                    EdgeIndexs[(iTri * 3) * 2 + 1] = (int)triArray.Tris[iTri].V[1];

                    EdgeIndexs[(iTri * 3 + 1) * 2 + 0] = (int)triArray.Tris[iTri].V[1];
                    EdgeIndexs[(iTri * 3 + 1) * 2 + 1] = (int)triArray.Tris[iTri].V[2];

                    EdgeIndexs[(iTri * 3 + 2) * 2 + 0] = (int)triArray.Tris[iTri].V[2];
                    EdgeIndexs[(iTri * 3 + 2) * 2 + 1] = (int)triArray.Tris[iTri].V[0];
                }
            }
        }

        public Mesher3DDrawPart(MeshBarArray3D barArray)
        {
            IsSelected = false;
            IsShown = true;
            Color = new double[3] { 0.8, 0.8, 0.8 };
            LineWidth = 1;
            ElemIndexs = null;
            EdgeIndexs = null;

            CadId = barArray.ECadId;

            MeshId = barArray.Id;
            System.Diagnostics.Debug.Assert(MeshId != 0);
            Type = ElementType.Line;

            ElemCount = (uint)barArray.Bars.Count;
            ElemIndexs = new int[ElemCount * 2];
            for (int ibar = 0; ibar < ElemCount; ibar++)
            {
                ElemIndexs[ibar * 2 + 0] = (int)barArray.Bars[ibar].V[0];
                ElemIndexs[ibar * 2 + 1] = (int)barArray.Bars[ibar].V[1];
            }
        }

        public Mesher3DDrawPart(MeshVertex3D vtx)
        {
            IsSelected = false;
            IsShown = true;
            Color = new double[3] { 0.8, 0.8, 0.8 };
            LineWidth = 1;
            ElemIndexs = null;
            EdgeIndexs = null;

            CadId = vtx.VCadId;

            MeshId = vtx.Id;
            System.Diagnostics.Debug.Assert(MeshId != 0);
            Type = ElementType.Point;

            ElemCount = 1;
            ElemIndexs = new int[ElemCount];
            ElemIndexs[0] = (int)vtx.V;
        }

        public void DrawElements()
        {
            if (Type == ElementType.Point)
            {
                DrawElementsVertex();
            }
            else if (Type == ElementType.Line)
            {
                DrawElementsBar();
            }
            else if (Type == ElementType.Tri)
            {
                DrawElementsTri();
            }
            else if (Type == ElementType.Tet)
            {
                DrawElementsTet();
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
        }

        public void DrawElementsSelection()
        {
            if (Type == ElementType.Point)
            {
                DrawElementsSelectionVertex();
            }
            else if (Type == ElementType.Line)
            {
                DrawElementsSelectionBar();
            }
            else if (Type == ElementType.Tri)
            {
                DrawElementsSelectionTri();
            }
            else if (Type == ElementType.Tet)
            {
                DrawElementsSelectionTet();
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
        }

        private void DrawElementsTet()
        {
            if (!IsShown)
            {
                return;
            }
            // 辺を描画
            GL.LineWidth(1);
            if (IsSelected)
            {
                GL.LineWidth(2);
                GL.Color3(1.0, 1.0, 0.0);
            }
            else
            {
                GL.LineWidth(1);
                GL.Color3(0.0, 0.0, 0.0);
            }

            GL.DrawElements(PrimitiveType.Lines, (int)EdgeCount * 2, DrawElementsType.UnsignedInt, EdgeIndexs);

            // ピックされた面を描画
            GL.Color3(1.0, 0.0, 0.0);
            for (uint iiElem = 0; iiElem < SelectedElems.Count; iiElem++)
            {
                uint iElem0 = SelectedElems[(int)iiElem];
                GL.Begin(PrimitiveType.Triangles);
                GL.ArrayElement(ElemIndexs[iElem0 * 3]);
                GL.ArrayElement(ElemIndexs[iElem0 * 3 + 1]);
                GL.ArrayElement(ElemIndexs[iElem0 * 3 + 2]);
                GL.End();
            }

            if (Mask != null)
            {
                //----------------------
                GL.Enable(EnableCap.PolygonStipple);
                GL.PolygonStipple(Mask);
                //----------------------
            }
            // 面を描画
            GL.Color3(Color);
            GL.DrawElements(PrimitiveType.Triangles, (int)ElemCount * 3, DrawElementsType.UnsignedInt, ElemIndexs);
            if (Mask != null)
            {
                //----------------------
                GL.Disable(EnableCap.PolygonStipple);
                //----------------------
            }
        }

        private void DrawElementsSelectionTet()
        {
            if (!IsShown)
            {
                return;
            }
            for (uint iElem = 0; iElem < ElemCount; iElem++)
            {
                GL.PushName(iElem);
                GL.Begin(PrimitiveType.Triangles);

                GL.ArrayElement(ElemIndexs[iElem * 3]);
                GL.ArrayElement(ElemIndexs[iElem * 3 + 1]);
                GL.ArrayElement(ElemIndexs[iElem * 3 + 2]);

                GL.End();
                GL.PopName();
            }
        }

        private void DrawElementsTri()
        {
            if (!IsShown)
            {
                return;
            }

            // 辺を描画
            /*
            if (IsSelected)
            {
                GL.LineWidth(LineWidth + 1);
                GL.Color3(1.0, 1.0, 0.0);
            }
            else
            {
                GL.LineWidth(LineWidth);
                GL.Color3(0.0, 0.0, 0.0);
            }
            */
            GL.LineWidth(LineWidth);
            GL.Color3(0.0, 0.0, 0.0);
            GL.DrawElements(PrimitiveType.Lines, (int)EdgeCount * 2, DrawElementsType.UnsignedInt, EdgeIndexs);

            GL.Color3(1.0, 0.0, 0.0);
            GL.Begin(PrimitiveType.Triangles);
            for (int iiElem = 0; iiElem < SelectedElems.Count; iiElem++)
            {
                uint iElem0 = SelectedElems[iiElem];

                GL.ArrayElement(ElemIndexs[iElem0 * 3]);
                GL.ArrayElement(ElemIndexs[iElem0 * 3 + 1]);
                GL.ArrayElement(ElemIndexs[iElem0 * 3 + 2]);
            }
            GL.End();

            if (Mask != null)
            {
                //----------------------
                GL.Enable(EnableCap.PolygonStipple);
                GL.PolygonStipple(Mask);
                //----------------------
            }
            // 面を描画
            GL.Color3(0.8, 0.8, 0.8);
            GL.DrawElements(PrimitiveType.Triangles, (int)ElemCount * 3, DrawElementsType.UnsignedInt, ElemIndexs);
            if (Mask != null)
            {
                //----------------------
                GL.Disable(EnableCap.PolygonStipple);
                //----------------------
            }
        }

        private void DrawElementsSelectionTri()
        {
            if (!IsShown)
            {
                return;
            }

            // 面を描画
            for (int iTri = 0; iTri < ElemCount; iTri++)
            {
                GL.PushName(iTri);
                GL.Begin(PrimitiveType.Triangles);

                GL.ArrayElement(ElemIndexs[iTri * 3]);
                GL.ArrayElement(ElemIndexs[iTri * 3 + 1]);
                GL.ArrayElement(ElemIndexs[iTri * 3 + 2]);

                GL.End();
                GL.PopName();
            }
        }

        private void DrawElementsBar()
        {
            if (!IsShown)
            {
                return;
            }

            GL.LineWidth(2);
            GL.Color3(1.0, 0.0, 0.0);
            GL.Begin(PrimitiveType.Lines);
            for (int iiElem = 0; iiElem < SelectedElems.Count; iiElem++)
            {
                uint iElem0 = SelectedElems[iiElem];

                GL.ArrayElement(ElemIndexs[iElem0 * 2]);
                GL.ArrayElement(ElemIndexs[iElem0 * 2 + 1]);
            }

            GL.End();

            /*
            if (IsSelected)
            {
                GL.Color3(1.0, 1.0, 0.0);
            }
            else
            {
                GL.Color3(0.0, 0.0, 0.0);
            }
            */
            GL.Color3(0.0, 0.0, 0.0);
            GL.DrawElements(PrimitiveType.Lines, (int)ElemCount * 2, DrawElementsType.UnsignedInt, ElemIndexs);
        }

        private void DrawElementsSelectionBar()
        {
            if (!IsShown)
            {
                return;
            }

            for (int ibar = 0; ibar < ElemCount; ibar++)
            {
                GL.PushName(ibar);
                GL.Begin(PrimitiveType.Lines);

                GL.ArrayElement(ElemIndexs[ibar * 2]);
                GL.ArrayElement(ElemIndexs[ibar * 2 + 1]);

                GL.End();
                GL.PopName();
            }
        }

        private void DrawElementsVertex()
        {
            if (!IsShown)
            {
                return;
            }

            GL.PointSize(5);
            GL.Color3(1.0, 0.0, 0.0);
            GL.Begin(PrimitiveType.Points);
            for (int iiElem = 0; iiElem < SelectedElems.Count; iiElem++)
            {
                uint iElem0 = SelectedElems[iiElem];

                GL.ArrayElement(ElemIndexs[iElem0]);
            }

            GL.End();

            /*
            if (IsSelected)
            {
                GL.Color3(1.0, 1.0, 0.0);
            }
            else
            {
                GL.Color3(0.0, 0.0, 0.0);
            }
            */
            GL.Color3(0.0, 0.0, 0.0);
            GL.DrawElements(PrimitiveType.Points, (int)ElemCount, DrawElementsType.UnsignedInt, ElemIndexs);
        }

        private void DrawElementsSelectionVertex()
        {
            if (!IsShown)
            {
                return;
            }

            GL.PushName(0);
            GL.Begin(PrimitiveType.Points);

            GL.ArrayElement(ElemIndexs[0]);
            GL.End();
            GL.PopName();
        }

    }
}
