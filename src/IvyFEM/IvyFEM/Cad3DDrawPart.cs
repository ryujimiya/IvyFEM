﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    public class Cad3DDrawPart
    {
        public int ShowMode { get; set; } = 0;  // -1:not_show 0:default 1:hilight 2:selected
        public double[] Color { get; set; } = new double[3] { 0.0, 0.0, 0.0 };

        public CadElementType Type { get; set; } = CadElementType.Vertex;
        public uint CadId { get; set; } = 0;
        public uint MeshId { get; set; } = 0;

        public uint ElemCount { get; set; } = 0;
        public uint ElemPtCount { get; set; } = 0;
        public uint[] Indexs { get; set; } = null;

        public double[] Normal { get; set; } = new double[3];

        public Cad3DDrawPart()
        {

        }

        public Cad3DDrawPart(Cad3DDrawPart src)
        {
            ShowMode = src.ShowMode;
            src.Color.CopyTo(Color, 0);
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
        }

        public void Clear()
        {
            Indexs = null;
            ElemCount = 0;
            ElemPtCount = 0;
        }

        public bool SetTriArray(MeshTriArray3D triArray)
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

            Color = new double[3] { 0.8, 0.8, 0.8 };

            return true;
        }

        public bool SetBarArray(MeshBarArray3D BarArray)
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

            Color = new double[3] { 0.0, 0.0, 0.0 };

            return true;
        }

        public bool SetVertex(MeshVertex3D vertex)
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

            Color = new double[3] { 0.0, 0.0, 0.0 };
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
