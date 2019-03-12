using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    public class CadObject2DDrawPart
    {
        public int ShowMode { get; set; } = 0;  // -1:not_show 0:default 1:hilight 2:selected
        public float[] Color { get; } = new float[3] { 0, 0, 0 };

        public CadElementType Type { get; set; } = CadElementType.Vertex;
        public uint CadId { get; set; } = 0;
        public uint MeshId { get; set; } = 0;

        public uint ElemCount { get; set; } = 0;
        public uint ElemPtCount { get; set; } = 0;
        public uint[] Indexs { get; set; } = null;

        public double Height { get; set; } = 0;
        public double DispX { get; set; } = 0;
        public double DispY { get; set; } = 0;

        public CurveType CurveType;
        public IList<OpenTK.Vector2d> CtrlPoints = new List<OpenTK.Vector2d>();

        public CadObject2DDrawPart()
        {

        }

        public CadObject2DDrawPart(CadObject2DDrawPart src)
        {
            ShowMode = src.ShowMode;
            for (int i = 0; i < 3; i++)
            {
                Color[i] = src.Color[i];
            }
            Type = src.Type;
            CadId = src.CadId;
            MeshId = src.MeshId;
            ElemCount = src.ElemCount;
            ElemPtCount = src.ElemPtCount;
            Indexs = null;
            if (src.Indexs != null)
            {
                Indexs = new uint[src.Indexs.Length];
                src.Indexs.CopyTo(Indexs, 0);
            }
            Height = src.Height;
            DispX = src.DispX;
            DispY = src.DispY;
            CurveType = src.CurveType;
            CtrlPoints.Clear();
            for (int i = 0; i < src.CtrlPoints.Count; i++)
            {
                CtrlPoints.Add(src.CtrlPoints[i]);
            }
        }

        public void Clear()
        {
            Indexs = null;
            ElemCount = 0;
            ElemPtCount = 0;
        }

        public bool SetTriArray(MeshTriArray2D triArray)
        {
            CadId = triArray.LCadId;
            System.Diagnostics.Debug.Assert(CadId != 0);
            MeshId = triArray.Id;
            System.Diagnostics.Debug.Assert(MeshId != 0);
            Type = CadElementType.Loop;

            ElemPtCount = 3;
            ElemCount = (uint)triArray.Tris.Count;
            Indexs = new uint[ElemCount * ElemPtCount];
            for (uint ielem = 0; ielem < ElemCount; ielem++)
            {
                for (uint ipoel = 0; ipoel < ElemPtCount; ipoel++)
                {
                    Indexs[ielem * ElemPtCount + ipoel] = triArray.Tris[(int)ielem].V[ipoel];
                }
            }

            /*
            Color[0] = 0.8f;
            Color[1] = 0.8f;
            Color[2] = 0.8f;
            */
            Color[0] = 0.2f;
            Color[1] = 0.2f;
            Color[2] = 0.2f;

            return true;
        }

        public bool SetBarArray(MeshBarArray BarArray)
        {
            CadId = BarArray.ECadId;
            System.Diagnostics.Debug.Assert(CadId != 0);
            MeshId = BarArray.Id;
            System.Diagnostics.Debug.Assert(MeshId != 0);
            Type = CadElementType.Edge;

            ElemPtCount = 2;
            ElemCount = (uint)BarArray.Bars.Count;
            Indexs = new uint[ElemCount * ElemPtCount];
            for (uint ielem = 0; ielem < ElemCount; ielem++)
            {
                for (uint ipoel = 0; ipoel < ElemPtCount; ipoel++)
                {
                    Indexs[ielem * ElemPtCount + ipoel] = BarArray.Bars[(int)ielem].V[ipoel];
                }
            }

            Color[0] = 0.0f;
            Color[1] = 0.0f;
            Color[2] = 0.0f;

            return true;
        }

        public bool SetVertex(MeshVertex vertex)
        {
            CadId = vertex.VCadId;
            System.Diagnostics.Debug.Assert(CadId != 0);
            MeshId = vertex.Id;
            System.Diagnostics.Debug.Assert(MeshId != 0);
            Type = CadElementType.Vertex;

            ElemPtCount = 1;
            ElemCount = 1;
            Indexs = new uint[ElemCount * ElemPtCount];
            Indexs[0] = vertex.V;

            Color[0] = 0.0f;
            Color[1] = 0.0f;
            Color[2] = 0.0f;
            return true;
        }

        public void DrawElements()
        {
            if (ElemPtCount == 1)
            {
                GL.DrawElements(PrimitiveType.Points, (int)(ElemCount * ElemPtCount), DrawElementsType.UnsignedInt, Indexs);
                return;
            }
            else if (ElemPtCount == 2)
            {
                GL.DrawElements(PrimitiveType.Lines, (int)(ElemCount * ElemPtCount), DrawElementsType.UnsignedInt, Indexs);
                return;
            }
            else if (ElemPtCount == 3)
            {
                GL.DrawElements(PrimitiveType.Triangles, (int)(ElemCount * ElemPtCount), DrawElementsType.UnsignedInt, Indexs);
                return;
            }
        }

    }
}
