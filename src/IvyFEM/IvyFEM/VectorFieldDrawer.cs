using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    public class VectorFieldDrawer : IFieldDrawer
    {
        private IList<VectorFieldDrawPart> DrawParts = new List<VectorFieldDrawPart>();
        private uint ValueId = 0;
        private FieldDerivativeType ValueDt = FieldDerivativeType.Value;

        public RotMode SutableRotMode { get; private set; } = RotMode.RotModeNotSet;
        public bool IsAntiAliasing { get; set; } = false;

        public VectorFieldDrawerType Type { get; private set; } = VectorFieldDrawerType.NotSet;

        public VectorFieldDrawer() : base()
        {

        }

        public VectorFieldDrawer(uint valueId, FieldDerivativeType valueDt, FEWorld world)
        {
            Set(valueId, valueDt, world);
        }

        private void Set(uint valueId, FieldDerivativeType valueDt, FEWorld world)
        {
            System.Diagnostics.Debug.Assert(world.IsFieldValueId(valueId));
            ValueId = valueId;
            var mesh = world.Mesh;

            uint dim = world.Dimension;
            {
                if (dim == 2)
                {
                    SutableRotMode = RotMode.RotMode2D;
                }
                else if (dim == 3)
                {
                    SutableRotMode = RotMode.RotMode3D;
                }
            }
            FieldValue fv = world.GetFieldValue(valueId);
            if (fv.Type == FieldValueType.Vector2 || fv.Type == FieldValueType.Vector3)
            {
                Type = VectorFieldDrawerType.Vector;
            }
            else if (fv.Type == FieldValueType.SymmetricTensor2)
            {
                Type = VectorFieldDrawerType.SymmetricTensor2;
            }

            {
                DrawParts.Clear();
                IList<uint> meshIds = mesh.GetIds();
                foreach (uint meshId in meshIds)
                {
                    VectorFieldDrawPart dp = new VectorFieldDrawPart(meshId, world);
                    DrawParts.Add(dp);
                }
            }

            Update(world);
        }

        public void Update(FEWorld world)
        {
            for (int idp = 0; idp < DrawParts.Count; idp++)
            {
                VectorFieldDrawPart dp = DrawParts[idp];
                dp.Update(ValueId, ValueDt, Type, world);
            }
        }

        public void Draw()
        {
            bool isTexture = GL.IsEnabled(EnableCap.Texture2D);
            //GL.Enable(EnableCap.DepthTest);
            GL.Disable(EnableCap.Texture2D);

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

            GL.LineWidth(1);
            GL.Begin(PrimitiveType.Lines);
            for (int idp = 0; idp < DrawParts.Count; idp++)
            {
                VectorFieldDrawPart dp = DrawParts[idp];
                int layer = dp.Layer;
                double height = (layer - minLayer) * layerHeight;
                GL.Translate(0, 0, +height);
                dp.DrawElements();
                GL.Translate(0, 0, -height);
            }
            GL.End();

            if (isTexture)
            {
                GL.Enable(EnableCap.Texture2D);
            }
        }

        public void DrawSelection(uint idraw)
        {
            throw new NotImplementedException();
        }

        public void AddSelected(int[] selectFlag)
        {
            throw new NotImplementedException();
        }

        public void ClearSelected()
        {
            throw new NotImplementedException();
        }

        public BoundingBox3D GetBoundingBox(OpenTK.Matrix3d rot)
        {
            //throw new NotImplementedException();
            return new BoundingBox3D();
        }

    }
}
