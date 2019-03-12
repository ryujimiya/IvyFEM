using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    public class FaceFieldDrawPart
    {
        public bool IsSelected { get; set; } = false;
        public uint MeshId { get; set; } = 0; 
        public ElementType Type { get; private set; } = ElementType.NotSet;
        public int Layer { get; set; } = 0;
        public double[] Color { get; } = new double[3] { 0.8, 0.8, 0.8 };
        public uint ElemCount { get; set; } = 0;
        public uint ElemPtCount { get; set; } = 0;
        public uint[] Indexs { get; set; } = null;
        public float[] Colors { get; set; } = null;

        public uint Dimension
        {
            get
            {
                if (Type == ElementType.Line) { return 1; }
                if (Type == ElementType.Tri || Type == ElementType.Quad) { return 2; }
                if (Type == ElementType.Tet || Type == ElementType.Hex) { return 3; }
                return 0;
            }
        }

        public FaceFieldDrawPart()
        {

        }

        public FaceFieldDrawPart(FaceFieldDrawPart src)
        {
            IsSelected = src.IsSelected;
            MeshId = src.MeshId;
            Type = src.Type;
            Layer = src.Layer;
            for (int i = 0; i < 3; i++)
            {
                Color[i] = src.Color[i];
            }
            ElemCount = src.ElemCount;
            ElemPtCount = src.ElemPtCount;
            Indexs = null;
            if (src.Indexs != null)
            {
                Indexs = new uint[src.Indexs.Length];
                src.Indexs.CopyTo(Indexs, 0);
            }
            Colors = null;
            if (src.Colors != null)
            {
                Colors = new float[src.Colors.Length];
                src.Colors.CopyTo(Colors, 0);
            }
        }

        public FaceFieldDrawPart(uint meshId, FEWorld world, uint valueId)
        {
            var mesh = world.Mesh;
            if (!mesh.IsId(meshId))
            {
                return;
            }
            MeshId = meshId;

            uint cadId;
            int layer;
            uint elemCount;
            MeshType meshType;
            int loc;
            mesh.GetInfo(MeshId, out cadId, out layer);
            mesh.GetMeshInfo(MeshId, out elemCount, out meshType, out loc, out cadId);
            Layer = layer;
            ElemCount = elemCount;

            if (meshType == MeshType.Vertex)
            {
                Type = ElementType.Point;
                Color[0] = 0;
                Color[1] = 0;
                Color[2] = 0;
            }
            else if (meshType == MeshType.Bar)
            {
                Type = ElementType.Line;
                Color[0] = 0;
                Color[1] = 0;
                Color[2] = 0;
                SetLine(world, valueId);
            }
            else if (meshType == MeshType.Tri)
            {
                Type = ElementType.Tri;
                SetTri(world, valueId);
            }
            else if (meshType == MeshType.Quad)
            {
                Type = ElementType.Quad;
                SetQuad(world);
            }
            else
            {
                throw new NotImplementedException();
            }
        }

        private void SetLine(FEWorld world, uint valueId)
        {
            System.Diagnostics.Debug.Assert(Type == ElementType.Line);
            if (Type != ElementType.Line)
            {
                return;
            }

            FieldValue fv = world.GetFieldValue(valueId);
            uint quantityId = fv.QuantityId;
            int feOrder;
            {
                uint feId = world.GetLineFEIdFromMesh(quantityId, MeshId, 0); // 先頭の要素
                System.Diagnostics.Debug.Assert(feId != 0);
                LineFE lineFE = world.GetLineFE(quantityId, feId);
                feOrder = lineFE.Order;
                ElemPtCount = lineFE.NodeCount;
            }

            Indexs = new uint[ElemPtCount * ElemCount];
            for (int iEdge = 0; iEdge < ElemCount; iEdge++)
            {
                uint feId = world.GetLineFEIdFromMesh(quantityId, MeshId, (uint)iEdge);
                System.Diagnostics.Debug.Assert(feId != 0);
                LineFE lineFE = world.GetLineFE(quantityId, feId);
                for (int iPt = 0; iPt < ElemPtCount; iPt++)
                {
                    Indexs[iEdge * ElemPtCount + iPt] = (uint)lineFE.NodeCoordIds[iPt];
                }
            }
        }

        private void SetTri(FEWorld world, uint valueId)
        {
            System.Diagnostics.Debug.Assert(Type == ElementType.Tri);
            if (Type != ElementType.Tri)
            {
                return;
            }

            FieldValue fv = world.GetFieldValue(valueId);
            uint quantityId = fv.QuantityId;
            int feOrder;
            {
                uint feId = world.GetTriangleFEIdFromMesh(quantityId, MeshId, 0); // 先頭の要素
                System.Diagnostics.Debug.Assert(feId != 0);
                TriangleFE triFE = world.GetTriangleFE(quantityId, feId);
                feOrder = triFE.Order;
                ElemPtCount = triFE.NodeCount;
            }

            Indexs = new uint[ElemPtCount * ElemCount];
            for (int iTri = 0; iTri < ElemCount; iTri++)
            {
                uint feId = world.GetTriangleFEIdFromMesh(quantityId, MeshId, (uint)iTri);
                System.Diagnostics.Debug.Assert(feId != 0);
                TriangleFE triFE = world.GetTriangleFE(quantityId, feId);
                for (int iPt = 0; iPt < ElemPtCount; iPt++)
                {
                    Indexs[iTri * ElemPtCount + iPt] = (uint)triFE.NodeCoordIds[iPt];
                }
            }
        }

        private void SetQuad(FEWorld world)
        {
            System.Diagnostics.Debug.Assert(Type == ElementType.Quad);
            if (Type != ElementType.Quad)
            {
                return;
            }
            throw new NotImplementedException();
        }

        public void SetColors(uint bubbleValueId, FieldDerivativeType dt, FEWorld world, IColorMap colorMap)
        {
            FieldValue fv = world.GetFieldValue(bubbleValueId);
            System.Diagnostics.Debug.Assert(fv.IsBubble == true);
            uint quantityId = fv.QuantityId;
            var mesh = world.Mesh;
            MeshType meshType;
            int[] vertexs;
            mesh.GetConnectivity(MeshId, out meshType, out vertexs);

            if (Type == ElementType.Tri)
            {
                Colors = new float[ElemCount * 3];
                for (int iTri = 0; iTri < ElemCount; iTri++)
                {
                    // Bubble
                    uint feId = world.GetTriangleFEIdFromMesh(quantityId, MeshId, (uint)iTri);
                    System.Diagnostics.Debug.Assert(feId != 0);
                    double value = fv.GetShowValue((int)(feId - 1), 0, dt);
                    var color = colorMap.GetColor(value);
                    for (int iColor = 0; iColor < 3; iColor++)
                    {
                        Colors[iTri * 3 + iColor] = (float)color[iColor];
                    }
                }
            }
            else if (Type == ElementType.Quad)
            {
                // TRIと同じでよいが要素IDを取得するメソッドが現状ない
                throw new NotImplementedException();
            }
        }

        public void ClearColors()
        {
            Colors = null;
        }

        public uint[] GetVertexs(uint iElem)
        {
            // 頂点だけを対象とする
            int vertexCnt = 0;
            if (Type == ElementType.Tri)
            {
                vertexCnt = 3;
            }
            uint[] vertexs = new uint[vertexCnt]; 
            for (int iPt = 0; iPt < vertexCnt; iPt++)
            {
                vertexs[iPt] = Indexs[iElem * ElemPtCount + iPt];
            }
            return vertexs;
        }

        public void DrawElements()
        {
            if (Colors == null)
            {
                //GL.Color3(Color[0], Color[1], Color[2]);
                if (Type == ElementType.Line)
                {
                    if (ElemPtCount == 2) // 1次要素
                    {
                        GL.DrawElements(PrimitiveType.Lines, (int)(ElemCount * ElemPtCount), DrawElementsType.UnsignedInt, Indexs);
                    }
                    else if (ElemPtCount == 3) // 2次要素
                    {
                        GL.Begin(PrimitiveType.Lines);
                        for (int iEdge = 0; iEdge < ElemCount; iEdge++)
                        {
                            GL.ArrayElement((int)Indexs[iEdge * ElemPtCount + 0]);
                            GL.ArrayElement((int)Indexs[iEdge * ElemPtCount + 2]);
                            GL.ArrayElement((int)Indexs[iEdge * ElemPtCount + 2]);
                            GL.ArrayElement((int)Indexs[iEdge * ElemPtCount + 1]);
                        }
                        GL.End();
                    }
                }
                else if (Type == ElementType.Tri || Type == ElementType.Tet) 
                {
                    if (ElemPtCount == 3) // 1次要素
                    {
                        GL.DrawElements(PrimitiveType.Triangles, (int)(ElemCount * ElemPtCount), DrawElementsType.UnsignedInt, Indexs);
                    }
                    else if (ElemPtCount == 6) // 2次要素
                    {
                        GL.Begin(PrimitiveType.Triangles);
                        for (int iTri = 0; iTri < ElemCount; iTri++)
                        {
                            // 4つの三角形
                            GL.ArrayElement((int)Indexs[iTri * ElemPtCount + 0]);
                            GL.ArrayElement((int)Indexs[iTri * ElemPtCount + 3]);
                            GL.ArrayElement((int)Indexs[iTri * ElemPtCount + 5]);
                            GL.ArrayElement((int)Indexs[iTri * ElemPtCount + 3]);
                            GL.ArrayElement((int)Indexs[iTri * ElemPtCount + 4]);
                            GL.ArrayElement((int)Indexs[iTri * ElemPtCount + 5]);
                            GL.ArrayElement((int)Indexs[iTri * ElemPtCount + 3]);
                            GL.ArrayElement((int)Indexs[iTri * ElemPtCount + 1]);
                            GL.ArrayElement((int)Indexs[iTri * ElemPtCount + 4]);
                            GL.ArrayElement((int)Indexs[iTri * ElemPtCount + 5]);
                            GL.ArrayElement((int)Indexs[iTri * ElemPtCount + 4]);
                            GL.ArrayElement((int)Indexs[iTri * ElemPtCount + 2]);
                        }
                        GL.End();
                    }
                }
                else if (Type == ElementType.Quad || Type == ElementType.Hex)
                {
                    throw new NotImplementedException();
                }
                return;
            }

            if (Type == ElementType.Tri || Type == ElementType.Tet)
            {
                // 要素全体を塗り潰すので要素の次数に関係なく1つの三角形を描画
                GL.Begin(PrimitiveType.Triangles);
                for (int iTri = 0; iTri < ElemCount; iTri++)
                {
                    GL.Color3(
                        Colors[iTri * 3], Colors[iTri * 3 + 1], Colors[iTri * 3 + 2]);
                    // 先頭の3点は頂点
                    for (int iPt = 0; iPt < 3; iPt++)
                    {
                        GL.ArrayElement((int)Indexs[iTri * ElemPtCount + iPt]);
                    }
                }
                GL.End();
            }
            else if (Type == ElementType.Quad || Type == ElementType.Hex)
            {
                throw new NotImplementedException();
            }
        }

    }
}
