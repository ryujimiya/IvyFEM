using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    public class FaceFieldDrawer : IFieldDrawer
    {
        private IList<FaceFieldDrawPart> DrawParts = new List<FaceFieldDrawPart>();
        private VertexArray VertexArray = new VertexArray();
        private uint ValueId = 0;
        private FieldDerivativeType ValueDt = FieldDerivativeType.Value;
        private uint ColorValueId = 0;
        private FieldDerivativeType ColorValueDt = FieldDerivativeType.Value;
        private bool IsntDisplacementValue = false;
        private bool IsNsvDraw = false;
        private float[] ColorArray = null;
        public IColorMap ColorMap { get; private set; } = new ColorMap();
        private bool IsColorLegendDraw = false;
        private double[] NormalArray = null;
        private double[] UVArray = null;
        public double TexCentX { get; set; } = 0;
        public double TexCentY { get; set; } = 0;
        public double TexScale { get; private set; } = 1;
        public RotMode SutableRotMode { get; private set; } = RotMode.RotModeNotSet;
        public bool IsAntiAliasing { get; set; } = false;

        public FaceFieldDrawer()
        {

        }

        public FaceFieldDrawer(uint valueId, FieldDerivativeType valueDt, bool isntDisplacementValue,
            FEWorld world, 
            uint colorValueId = 0, FieldDerivativeType colorValueDt = FieldDerivativeType.Value)
        {
            Set(valueId, valueDt, isntDisplacementValue, world, colorValueId, colorValueDt);
        }

        public FaceFieldDrawer(uint valueId, FieldDerivativeType valueDt, bool isntDisplacementValue,
            FEWorld world,
            uint colorValueId, FieldDerivativeType colorValueDt,
            double min, double max)
        {
            ColorMap = new ColorMap(min, max);
            Set(valueId, valueDt, isntDisplacementValue, world, colorValueId, colorValueDt);
        }

        private void Set(uint valueId, FieldDerivativeType valueDt, bool isntDisplacementValue,
            FEWorld world,
            uint colorValueId, FieldDerivativeType colorValueDt)
        {
            var mesh = world.Mesh;

            if (!world.IsFieldValueId(valueId))
            {
                throw new ArgumentException();
                //return;
            }
            if (world.IsFieldValueId(colorValueId))
            {
                IsColorLegendDraw = true;
            }
            else
            {
                IsColorLegendDraw = false;
            }

            ValueId = valueId;
            ValueDt = valueDt;
            ColorValueId = colorValueId;
            ColorValueDt = colorValueDt;
            IsntDisplacementValue = isntDisplacementValue;

            var fv = world.GetFieldValue(valueId);
            uint quantityId = fv.QuantityId;
            uint dim = world.Dimension;
            uint ptCnt = 0;
            if (IsNsvDraw)
            {
                ptCnt = fv.GetPointCount();
                // あとで実装する
                throw new NotImplementedException();
            }
            else
            {
                ptCnt = world.GetCoordCount(quantityId);
            }

            uint drawDim;
            if (!IsntDisplacementValue 
                && dim == 2 
                && (fv.Type == FieldValueType.Scalar ||fv.Type== FieldValueType.ZScalar))
            {
                drawDim = 3;
            }
            else
            {
                drawDim = dim;
            }
            VertexArray.SetSize(ptCnt, drawDim);

            { 
                bool isNormal = (NormalArray != null);
                if (NormalArray != null)
                {
                    NormalArray = null;
                }
                if (isNormal)
                {
                    NormalArray = new double[ptCnt * 3];
                }
            }

            {
                bool isUv = (UVArray != null);
                if (UVArray != null)
                {
                    UVArray = null;
                }
                if (isUv)
                {
                    UVArray = new double[ptCnt * 2];
                }
            }

            if (drawDim == 2) { SutableRotMode = RotMode.RotMode2D; }
            else if (dim == 3) { SutableRotMode =RotMode.RotMode3D; }
            else { SutableRotMode = RotMode.RotMode2DH; }

            {
                DrawParts.Clear();
                IList<uint> meshIds = mesh.GetIds();
                foreach(uint meshId in meshIds)
                {
                    if (IsNsvDraw)
                    {
                        // あとで実装する
                        throw new NotImplementedException();
                    }
                    else
                    {

                    }
                    FaceFieldDrawPart dp = new FaceFieldDrawPart(meshId, world, valueId);
                    DrawParts.Add(dp);
                }
            }

            Update(world);
        }

        public void Update(FEWorld world)
        {
            FieldValue fv = world.GetFieldValue(ValueId);
            uint quantityId = fv.QuantityId;
            uint dim = world.Dimension;

            uint ptCnt = 0;
            if (IsNsvDraw)
            {
                ptCnt = fv.GetPointCount();
                // あとで実装する
                throw new NotImplementedException();
            }
            else
            {
                ptCnt = world.GetCoordCount(quantityId);
            }
            System.Diagnostics.Debug.Assert(VertexArray.PointCount == ptCnt);

            if (!IsntDisplacementValue)
            {
                // 変位を伴う場合

                if (dim == 2 &&
                    (fv.Type == FieldValueType.Scalar || fv.Type == FieldValueType.ZScalar))
                {
                    // 垂直方向の変位として捉える

                    System.Diagnostics.Debug.Assert(VertexArray.Dimension == 3);
                    for (int coId = 0; coId < ptCnt; coId++)
                    {
                        double[] coord = world.GetCoord(quantityId, coId);
                        FieldDerivativeType dt = ValueDt;
                        double value = fv.GetShowValue(coId, 0, dt);
                        VertexArray.VertexCoords[coId * 3 + 0] = coord[0];
                        VertexArray.VertexCoords[coId * 3 + 1] = coord[1];
                        VertexArray.VertexCoords[coId * 3 + 2] = value;
                    }
                }
                else
                {
                    System.Diagnostics.Debug.Assert(VertexArray.Dimension == dim);
                    for (int coId = 0; coId < ptCnt; coId++)
                    {
                        double[] coord = world.GetCoord(quantityId, coId);
                        FieldDerivativeType dt = ValueDt;
                        for (int iDim = 0; iDim < dim; iDim++)
                        {
                            double value = fv.GetShowValue(coId, iDim, dt);
                            VertexArray.VertexCoords[coId * dim + iDim] = coord[iDim] + value;
                        }
                    }
                }
            }
            else
            {
                System.Diagnostics.Debug.Assert(VertexArray.Dimension == dim);
                for (int coId = 0; coId < ptCnt; coId++)
                {
                    double[] coord = world.GetCoord(quantityId, coId);
                    for (int iDim = 0; iDim < dim; iDim++)
                    {
                        VertexArray.VertexCoords[coId * dim + iDim] = coord[iDim];
                    }
                }
            }

            if (world.IsFieldValueId(ColorValueId))
            {
                FieldValue colorfv = world.GetFieldValue(ColorValueId);
                FieldDerivativeType dt = ColorValueDt;
                bool isBubble = colorfv.IsBubble;
                uint colorQuantityId = colorfv.QuantityId;
                uint colorPtCnt = world.GetCoordCount(colorQuantityId);
                if (!ColorMap.IsFixedMinMax)
                {
                    double min;
                    double max;

                    colorfv.GetMinMaxShowValue(out min, out max, 0, dt);
                    ColorMap.MinValue = min;
                    ColorMap.MaxValue = max;
                }

                if (!isBubble)
                {
                    // color is assigned to vertex
                    if (ColorArray == null)
                    {
                        ColorArray = new float[colorPtCnt * 4];
                    }
                    for (int coId = 0; coId < colorPtCnt; coId++)
                    {
                        double value = colorfv.GetShowValue(coId, 0, dt);
                        double[] color = ColorMap.GetColor(value);
                        ColorArray[coId * 4] = (float)color[0];
                        ColorArray[coId * 4 + 1] = (float)color[1];
                        ColorArray[coId * 4 + 2] = (float)color[2];
                        ColorArray[coId * 4 + 3] = 0.0f;
                    }
                }
                else
                {
                    // color is assigned to face
                    for (int idp = 0; idp < DrawParts.Count; idp++)
                    {
                        FaceFieldDrawPart dp = DrawParts[idp];
                        dp.SetColors(ColorValueId, dt, world, ColorMap);
                    }
                }
            }

            if (NormalArray != null)
            {
                MakeNormal();
            }


            if (UVArray != null)
            {
                for (int coId = 0; coId < ptCnt; coId++)
                {
                    double[] coord = world.GetCoord(quantityId, coId);
                    UVArray[coId * 2 + 0] = coord[0] * TexScale;
                    UVArray[coId * 2 + 1] = coord[1] * TexScale;
                }
            }
        }

        private void MakeNormal()
        {
            uint ptCnt = VertexArray.PointCount;
            if (NormalArray == null)
            {
                NormalArray = new double[ptCnt * 3];
            }
            for (int i = 0; i < 3 * ptCnt; i++)
            {
                NormalArray[i] = 0;
            }
            for (int idp = 0; idp < DrawParts.Count; idp++)
            {
                FaceFieldDrawPart dp = DrawParts[idp];
                uint elemCnt = dp.ElemCount;
                for (int iElem = 0; iElem < elemCnt; iElem++)
                {
                    uint[] vertexs = dp.GetVertexs((uint)iElem);
                    OpenTK.Vector3d c0 = new OpenTK.Vector3d(
                        VertexArray.VertexCoords[vertexs[0] * 3],
                        VertexArray.VertexCoords[vertexs[0] * 3 + 1],
                        VertexArray.VertexCoords[vertexs[0] * 3 + 2]);
                    OpenTK.Vector3d c1 = new OpenTK.Vector3d(
                        VertexArray.VertexCoords[vertexs[1] * 3],
                        VertexArray.VertexCoords[vertexs[1] * 3 + 1],
                        VertexArray.VertexCoords[vertexs[1] * 3 + 2]);
                    OpenTK.Vector3d c2 = new OpenTK.Vector3d(
                        VertexArray.VertexCoords[vertexs[2] * 3],
                        VertexArray.VertexCoords[vertexs[2] * 3 + 1],
                        VertexArray.VertexCoords[vertexs[2] * 3 + 2]);
                    double[] n;
                    double area;
                    CadUtils.UnitNormalAreaTri3D(out n, out area, c0, c1, c2);
                    NormalArray[vertexs[0] * 3 + 0] += n[0];
                    NormalArray[vertexs[0] * 3 + 1] += n[1];
                    NormalArray[vertexs[0] * 3 + 2] += n[2];
                    NormalArray[vertexs[1] * 3 + 0] += n[0];
                    NormalArray[vertexs[1] * 3 + 1] += n[1];
                    NormalArray[vertexs[1] * 3 + 2] += n[2];
                    NormalArray[vertexs[2] * 3 + 0] += n[0];
                    NormalArray[vertexs[2] * 3 + 1] += n[1];
                    NormalArray[vertexs[2] * 3 + 2] += n[2];
                }
            }
            for (int iPt = 0; iPt < ptCnt; iPt++)
            {
                double[] p = new double[3];
                p[0] = NormalArray[iPt * 3];
                p[1] = NormalArray[iPt * 3 + 1];
                p[2] = NormalArray[iPt * 3 + 2];
                double invLen = 1.0 / Math.Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
                NormalArray[0] *= invLen;
                NormalArray[1] *= invLen;
                NormalArray[2] *= invLen;
            }
        }

        public void EnableNormal(bool isLighting)
        {
            if ((NormalArray != null) == isLighting)
            {
                return;
            }
            if (isLighting)
            {
                uint ptCnt = VertexArray.PointCount;
                NormalArray = new double[ptCnt * 3];
                MakeNormal();
            }
            else
            {
                NormalArray = null;
            }
        }

        public void EnableUVMap(bool isUVMap, FEWorld world)
        {
            if ((UVArray != null) == isUVMap)
            {
                return;
            }
            if (isUVMap)
            {
                uint ptCnt = VertexArray.PointCount;
                UVArray = new double[ptCnt * 2];
                if (!world.IsFieldValueId(ValueId))
                {
                    return;
                }
                FieldValue fv = world.GetFieldValue(ValueId);
                uint quantityId = fv.QuantityId;
                uint coCnt = world.GetCoordCount(quantityId);
                for (int coId = 0; coId < coCnt; coId++)
                {
                    double[] c = world.GetCoord(quantityId, coId);
                    UVArray[coId * 2 + 0] = c[0] * TexScale;
                    UVArray[coId * 2 + 1] = c[1] * TexScale;
                }
            }
            else
            {
                UVArray = null;
            }
        }

        public void SetTexScale(double scale, FEWorld world)
        {
            TexScale = scale;

            if (UVArray != null)
            {
                if (!world.IsFieldValueId(ValueId))
                {
                    return;
                }
                FieldValue fv = world.GetFieldValue(ValueId);
                uint quantityId = fv.QuantityId;
                uint coCnt = world.GetCoordCount(quantityId);
                for (int coId = 0; coId < coCnt; coId++)
                {
                    double[] c = world.GetCoord(quantityId, coId);
                    UVArray[coId * 2 + 0] = c[0] * TexScale;
                    UVArray[coId * 2 + 1] = c[1] * TexScale;
                }
            }
        }

        public void AddSelected(int[] selectFlag)
        {
            throw new NotImplementedException();
        }

        public void ClearSelected()
        {
            throw new NotImplementedException();
        }

        public void Draw()
        {
            if (VertexArray.Dimension == 2)
            {
                // cannot see the opposite side
                GL.Enable(EnableCap.CullFace);
                GL.CullFace(CullFaceMode.Back);
            }
            else
            {
                GL.Disable(EnableCap.CullFace);
            }

            int minLayer;
            int maxLayer;
            {
                if (DrawParts.Count > 0)
                {
                    minLayer = DrawParts[0].Layer;
                    maxLayer = minLayer;
                }
                else
                {
                    minLayer = 0; maxLayer = 0;
                }
                for (int idp = 1; idp < DrawParts.Count; idp++)
                {
                    int layer = DrawParts[idp].Layer;
                    minLayer = (layer < minLayer) ? layer : minLayer;
                    maxLayer = (layer > maxLayer) ? layer : maxLayer;
                }
            }
            double layerHeight = 1.0 / (maxLayer - minLayer + 1);

            if (ColorArray == null)
            {
                // color is assigned to face
                GL.LineWidth(3);
                GL.EnableClientState(ArrayCap.VertexArray);
                GL.VertexPointer((int)VertexArray.Dimension, VertexPointerType.Double,
                    0, VertexArray.VertexCoords);
                if (NormalArray != null && GL.IsEnabled(EnableCap.Lighting))
                {
                    GL.EnableClientState(ArrayCap.NormalArray);
                    GL.NormalPointer(NormalPointerType.Double, 0, NormalArray);
                    float[] shine = { 0, 0, 0, 0 };
                    GL.Material(MaterialFace.FrontAndBack, MaterialParameter.Specular, shine);
                    GL.Material(MaterialFace.FrontAndBack, MaterialParameter.Shininess, 20.0f);
                }
                if (UVArray != null && GL.IsEnabled(EnableCap.Texture2D))
                {
                    GL.EnableClientState(ArrayCap.TextureCoordArray);
                    GL.TexCoordPointer(2, TexCoordPointerType.Double, 0, UVArray);
                    GL.MatrixMode(MatrixMode.Texture);
                    GL.LoadIdentity();
                    GL.Translate(-TexCentX, -TexCentY, 0.0);
                    GL.MatrixMode(MatrixMode.Modelview);
                }
                for (int idp = 0; idp < DrawParts.Count; idp++)
                {
                    FaceFieldDrawPart dp = DrawParts[idp];
                    if (dp.Dimension == 1)
                    {
                        // draw line
                        GL.Color3(0.0, 0.0, 0.0);
                        GL.LineWidth(3);
                    }
                    if (dp.Dimension == 2 || dp.Dimension == 3)
                    {  
                        // draw face
                        GL.Color3(dp.Color);
                    }
                    float[] color1 = { (float)dp.Color[0], (float)dp.Color[1], (float)dp.Color[2] };
                    GL.Material(MaterialFace.FrontAndBack, MaterialParameter.Diffuse, color1);
                    //GL.Color3(0.8, 0.8, 0.8);      
                    dp.DrawElements();
                }

                GL.DisableClientState(ArrayCap.VertexArray);
                GL.DisableClientState(ArrayCap.NormalArray);
                GL.DisableClientState(ArrayCap.TextureCoordArray);
            }
            else
            {
                // color is assigned to vertex
                GL.ShadeModel(ShadingModel.Smooth);
                GL.EnableClientState(ArrayCap.VertexArray);
                GL.VertexPointer((int)VertexArray.Dimension, VertexPointerType.Double,
                    0, VertexArray.VertexCoords);
                if (NormalArray != null && GL.IsEnabled(EnableCap.Lighting))
                {
                    GL.EnableClientState(ArrayCap.NormalArray);
                    GL.NormalPointer(NormalPointerType.Double, 0, NormalArray);
                    GL.Enable(EnableCap.ColorMaterial);
                    GL.ColorMaterial(MaterialFace.FrontAndBack, ColorMaterialParameter.AmbientAndDiffuse);
                    float[] shine = { 1, 1, 1, 1 };
                    GL.Material(MaterialFace.FrontAndBack, MaterialParameter.Specular, shine);
                    GL.Material(MaterialFace.FrontAndBack, MaterialParameter.Shininess, 100.0f);
                }
                if (UVArray != null && GL.IsEnabled(EnableCap.Texture2D))
                {
                    GL.EnableClientState(ArrayCap.TextureCoordArray);
                    GL.TexCoordPointer(2, TexCoordPointerType.Double, 0, UVArray);
                }
                GL.EnableClientState(ArrayCap.ColorArray);
                GL.ColorPointer(4, ColorPointerType.Float, 0, ColorArray);

                for (int idp = 0; idp < DrawParts.Count; idp++)
                {
                    FaceFieldDrawPart dp = DrawParts[idp];
                    int layer = dp.Layer;
                    double height = (layer - minLayer) * layerHeight;
                    GL.Translate(0, 0, +height);
                    dp.DrawElements();
                    GL.Translate(0, 0, -height);
                }

                GL.DisableClientState(ArrayCap.ColorArray);
                GL.DisableClientState(ArrayCap.VertexArray);
                GL.DisableClientState(ArrayCap.NormalArray);
                GL.DisableClientState(ArrayCap.TextureCoordArray);
            }

            if (IsColorLegendDraw)
            {
                // draw legend
                var legend = new ColorLegend(ColorMap);
                legend.Draw();
            }
        }

        public void DrawSelection(uint idraw)
        {
            throw new NotImplementedException();
        }

        public BoundingBox3D GetBoundingBox(OpenTK.Matrix3d rot)
        {
            return VertexArray.GetBoundingBox(rot);
        }

    }
}
