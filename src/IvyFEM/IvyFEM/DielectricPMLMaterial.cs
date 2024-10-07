using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class DielectricPMLMaterial : Material
    {
        public double Epxx { get => Values[0]; set => Values[0] = value; }
        public double Epyy { get => Values[1]; set => Values[1] = value; }
        public double Epzz { get => Values[2]; set => Values[2] = value; }
        public double Muxx { get => Values[3]; set => Values[3] = value; }
        public double Muyy { get => Values[4]; set => Values[4] = value; }
        public double Muzz { get => Values[5]; set => Values[5] = value; }
        public double Reflection0 { get => Values[6]; set => Values[6] = value; }
        public double DistanceOrder { get => Values[7]; set => Values[7] = value; }
        protected double OriginPointX { get => Values[8]; set => Values[8] = value; }
        protected double OriginPointY { get => Values[9]; set => Values[9] = value; }
        protected double OriginPointZ { get => Values[10]; set => Values[10] = value; } // 未使用
        public OpenTK.Vector2d OriginPoint
        {
            get
            {
                return new OpenTK.Vector2d(OriginPointX, OriginPointY);
            }
            set
            {
                OriginPointX = value.X;
                OriginPointY = value.Y;
            }
        }
        public OpenTK.Vector3d OriginPoint3D
        {
            get
            {
                return new OpenTK.Vector3d(OriginPointX, OriginPointY, OriginPointZ);
            }
            set
            {
                OriginPointX = value.X;
                OriginPointY = value.Y;
                OriginPointZ = value.Z;
            }
        }
        public double XThickness { get => Values[11]; set => Values[11] = value; }
        public double YThickness { get => Values[12]; set => Values[12] = value; }
        public double ZThickness { get => Values[13]; set => Values[13] = value; }
        // 回転
        protected double RotOriginPointX { get => Values[14]; set => Values[14] = value; }
        protected double RotOriginPointY { get => Values[15]; set => Values[15] = value; }
        protected double RotOriginPointZ { get => Values[16]; set => Values[16] = value; }
        public OpenTK.Vector2d RotOriginPoint
        {
            get
            {
                return new OpenTK.Vector2d(RotOriginPointX, RotOriginPointY);
            }
            set
            {
                RotOriginPointX = value.X;
                RotOriginPointY = value.Y;
            }
        }
        public OpenTK.Vector3d RotOriginPoint3D
        {
            get
            {
                return new OpenTK.Vector3d(RotOriginPointX, RotOriginPointY, RotOriginPointZ);
            }
            set
            {
                RotOriginPointX = value.X;
                RotOriginPointY = value.Y;
            }
        }
        public double RotAngle { get => Values[14]; set => Values[14] = value; }

        public DielectricPMLMaterial() : base()
        {
            int len = 17;
            Values = new double[len];
            for (int i = 0; i < len; i++)
            {
                Values[i] = 0.0;
            }
            Epxx = 1.0;
            Epyy = 1.0;
            Epzz = 1.0;
            Muxx = 1.0;
            Muyy = 1.0;
            Muzz = 1.0;
            Reflection0 = 1.0e-8;
            DistanceOrder = 2.0;

            int intLen = 1;
            IntValues = new int[intLen];
            for (int i = 0; i < intLen; i++)
            {
                IntValues[i] = 0;
            }
        }

        public DielectricPMLMaterial(DielectricMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }

        public bool IsXDirection()
        {
            return Math.Abs(XThickness) >= Constants.PrecisionLowerLimit;
        }

        public bool IsYDirection()
        {
            return Math.Abs(YThickness) >= Constants.PrecisionLowerLimit;
        }

        public bool IsZDirection()
        {
            return Math.Abs(ZThickness) >= Constants.PrecisionLowerLimit;
        }

        private double CalcSigma(
            double thickness,
            double ep, double mu,
            double pos, double origin)
        {
            double order = DistanceOrder;
            // σmax
            System.Diagnostics.Debug.Assert(thickness >= 0);
            double sigmaMax = ((1.0 + order) / (2.0 * thickness)) *
                Math.Sqrt(Constants.Ep0 / Constants.Mu0) *
                Math.Sqrt(ep / mu) *
                Math.Log(1.0 / Reflection0);

            // σ
            double s = Math.Abs(pos - origin) / thickness;
            System.Diagnostics.Debug.Assert(s <= 1.0);

            double sigma = sigmaMax * Math.Pow(s, order);

            return sigma;
        }

        public double CalcSigmaX(OpenTK.Vector2d pt)
        {
            double sigmaX = 0.0;
            if (!IsXDirection())
            {
                return sigmaX;
            }

            double[] pmlOriginCoord = { OriginPointX, OriginPointY };
            double[] rotOriginCoord = { RotOriginPointX, RotOriginPointY };
            pmlOriginCoord = CadUtils2D.GetRotCoord2D(pmlOriginCoord, RotAngle, rotOriginCoord);
            double thickness = XThickness;
            double ep = Epxx;
            double mu = Muxx;
            double pos = pt.X;
            double origin = pmlOriginCoord[0]; // X
            sigmaX = CalcSigma(
                thickness,
                ep, mu,
                pos, origin);

            return sigmaX;
        }

        public double CalcSigmaY(OpenTK.Vector2d pt)
        {
            double sigmaY = 0.0;
            if (!IsYDirection())
            {
                return sigmaY;
            }

            double[] pmlOriginCoord = { OriginPointX, OriginPointY };
            double[] rotOriginCoord = { RotOriginPointX, RotOriginPointY };
            pmlOriginCoord = CadUtils2D.GetRotCoord2D(pmlOriginCoord, RotAngle, rotOriginCoord);
            double thickness = YThickness;
            double ep = Epyy;
            double mu = Muyy;
            double pos = pt.Y;
            double origin = pmlOriginCoord[1]; // Y
            sigmaY = CalcSigma(
                thickness,
                ep, mu,
                pos, origin);

            return sigmaY;
        }

        // Time Domain
        private void CalcSigmaForTD(
            double thickness,
            double ep, double mu,
            double pos, double origin,
            double dt,
            out double sigma,
            out double c1P, out double c2P,
            out double c1V)
        {
            sigma = CalcSigma(thickness, ep, mu, pos, origin);

            // c1ψ,c2ψ
            c1P = 1.0 - Math.Exp(-1.0 * (sigma / (Constants.Ep0 * ep)) * dt);
            c2P = Math.Exp(-1.0 * (sigma / (Constants.Ep0 * ep)) * dt);

            // c1v
            c1V = (sigma / (Constants.Ep0 * ep)) * dt;
        }

        // Time Domain
        public void CalcSigmaXForTD(
            OpenTK.Vector2d pt,
            double dt,
            out double sigmaX,
            out double c1PX, out double c2PX,
            out double c1VX)
        {
            sigmaX = 0.0;
            c1PX = 0.0;
            c2PX = 0.0;
            c1VX = 0.0;
            if (!IsXDirection())
            {
                return;
            }

            double[] pmlOriginCoord = { OriginPointX, OriginPointY };
            double[] rotOriginCoord = { RotOriginPointX, RotOriginPointY };
            pmlOriginCoord = CadUtils2D.GetRotCoord2D(pmlOriginCoord, RotAngle, rotOriginCoord);
            double thickness = XThickness;
            double ep = Epxx;
            double mu = Muxx;
            double pos = pt.X;
            double origin = pmlOriginCoord[0]; // X
            CalcSigmaForTD(
                thickness,
                ep, mu,
                pos, origin,
                dt,
                out sigmaX,
                out c1PX, out c2PX,
                out c1VX);
        }

        // Time Domain
        public void CalcSigmaYForTD(
            OpenTK.Vector2d pt,
            double dt,
            out double sigmaY,
            out double c1PY, out double c2PY,
            out double c1VY)
        {
            sigmaY = 0.0;
            c1PY = 0.0;
            c2PY = 0.0;
            c1VY = 0.0;
            if (!IsYDirection())
            {
                return;
            }

            double[] pmlOriginCoord = { OriginPointX, OriginPointY };
            double[] rotOriginCoord = { RotOriginPointX, RotOriginPointY };
            pmlOriginCoord = CadUtils2D.GetRotCoord2D(pmlOriginCoord, RotAngle, rotOriginCoord);
            double thickness = YThickness;
            double ep = Epyy;
            double mu = Muyy;
            double pos = pt.Y;
            double origin = pmlOriginCoord[1]; // Y
            CalcSigmaForTD(
                thickness,
                ep, mu,
                pos, origin,
                dt,
                out sigmaY,
                out c1PY, out c2PY,
                out c1VY);
        }

        public double CalcSigmaXFor3D(OpenTK.Vector3d pt)
        {
            double sigmaX = 0.0;
            if (!IsXDirection())
            {
                return sigmaX;
            }

            double[] pmlOriginCoord = { OriginPointX, OriginPointY, OriginPointZ };
            //double[] rotOriginCoord = { RotOriginPointX, RotOriginPointY, RotOriginPointZ };
            //pmlOriginCoord = CadUtils2D.GetRotCoord2D(pmlOriginCoord, RotAngle, rotOriginCoord);
            double thickness = XThickness;
            double ep = Epxx;
            double mu = Muxx;
            double pos = pt.X;
            double origin = pmlOriginCoord[0]; // X
            sigmaX = CalcSigma(
                thickness,
                ep, mu,
                pos, origin);

            return sigmaX;
        }

        /*同じなので削除
        // 3D Time Domain
        // X方向PMLのみ対応
        private void CalcSigmaFor3DTD(
            double thickness,
            double ep, double mu,
            double pos, double origin,
            double dt,
            out double sigma,
            out double c1P, out double c2P,
            out double c1V)
        {
            sigma = CalcSigma(thickness, ep, mu, pos, origin);

            // c1ψ,c2ψ
            c1P = 1.0 - Math.Exp(-1.0 * (sigma / (Constants.Ep0 * ep)) * dt);
            c2P = Math.Exp(-1.0 * (sigma / (Constants.Ep0 * ep)) * dt);

            // c1v
            c1V = (sigma / (Constants.Ep0 * ep)) * dt;
        }
        */

        // 3D Time Domain
        public void CalcSigmaXFor3DTD(
            OpenTK.Vector3d pt,
            double dt,
            out double sigmaX,
            out double c1PX, out double c2PX,
            out double c1VX)
        {
            sigmaX = 0.0;
            c1PX = 0.0;
            c2PX = 0.0;
            c1VX = 0.0;
            if (!IsXDirection())
            {
                return;
            }

            double[] pmlOriginCoord = { OriginPointX, OriginPointY, OriginPointZ };
            //double[] rotOriginCoord = { RotOriginPointX, RotOriginPointY };
            //pmlOriginCoord = CadUtils2D.GetRotCoord2D(pmlOriginCoord, RotAngle, rotOriginCoord);
            double thickness = XThickness;
            double ep = Epxx;
            double mu = Muxx;
            double pos = pt.X;
            double origin = pmlOriginCoord[0]; // X
            /*同じなので削除
            CalcSigmaFor3DTD(
                thickness,
                ep, mu,
                pos, origin,
                dt,
                out sigmaX,
                out c1PX, out c2PX,
                out c1VX);
            */
            CalcSigmaForTD(
                thickness,
                ep, mu,
                pos, origin,
                dt,
                out sigmaX,
                out c1PX, out c2PX,
                out c1VX);
        }

    }
}
