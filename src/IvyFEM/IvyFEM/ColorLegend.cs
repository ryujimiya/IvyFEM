using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    public class ColorLegend
    {
        public IColorMap ColorMap { get; set; } = null;
        public double ScaleX { get; set; } = 0.03;
        public double ScaleY { get; set; } = 0.05;
        public double FontSize { get; set; } = 23.5;

        public ColorLegend()
        {

        }

        public ColorLegend(IColorMap colorMap)
        {
            ColorMap = colorMap;
        }

        public void Draw()
        {
            GL.ShadeModel(ShadingModel.Smooth);
            GL.Disable(EnableCap.CullFace);
            GL.LineWidth(1);
            GL.Color3(0, 0, 0);

            int[] viewport = new int[4];
            GL.GetInteger(GetPName.Viewport, viewport);
            double asp = (viewport[2] + 1.0) / (viewport[3] + 1.0);

            GL.MatrixMode(MatrixMode.Projection);
            GL.PushMatrix();
            GL.LoadIdentity();
            GL.Ortho(-asp, asp, -1, 1, -1, 1);

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();
            GL.LoadIdentity();

            double intervalN = 17.0;
            uint divC = 20;
            GL.Scale(ScaleX, ScaleY, 1.0);
            GL.Translate(((asp - 0.05) / ScaleX - 9), ((1 - 0.05) / ScaleY - 10.5 * 1.7), 1.0);
            GL.Begin(PrimitiveType.Quads);
            for (int i = 0; i < divC; i++)
            {
                double val0 = (ColorMap.MaxValue - ColorMap.MinValue) * (1.0 / divC) * (i + 0) + ColorMap.MinValue;
                double val1 = (ColorMap.MaxValue - ColorMap.MinValue) * (1.0 / divC) * (i + 1) + ColorMap.MinValue;
                double[] color0 = ColorMap.GetColor(val0);
                double[] color1 = ColorMap.GetColor(val1);
                GL.Color3(color0);
                GL.Vertex2(-3, intervalN * (1.0 / divC) * i);
                GL.Vertex2(-0, intervalN * (1.0 / divC) * i);
                GL.Color3(color1);
                GL.Vertex2(-0, intervalN * (1.0 / divC) * (i + 1));
                GL.Vertex2(-3, intervalN * (1.0 / divC) * (i + 1));
            }
            GL.End();

            GL.Color3(0, 0, 0);
            GL.Translate(0, -0.5, 0);
            uint divN = 10;
            for (int i = 0; i < divN + 1; i++)
            {
                double val = (ColorMap.MaxValue - ColorMap.MinValue) * i / divN + ColorMap.MinValue;
                string str = string.Format("{0,5:N1}", val);
                OpenGLUtils.DrawString(str, FontSize);
                GL.Translate(0, +intervalN / divN, 0);
            }
            GL.MatrixMode(MatrixMode.Projection);
            GL.PopMatrix();
            GL.MatrixMode(MatrixMode.Modelview);
            GL.PopMatrix();
        }


    }
}
