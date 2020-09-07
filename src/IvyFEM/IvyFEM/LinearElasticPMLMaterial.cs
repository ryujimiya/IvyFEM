using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class LinearElasticPMLMaterial : ElasticBaseMaterial
    {
        // PML
        public double Reflection0 { get => Values[5]; set => Values[5] = value; }
        public double ScalingFactor0 { get => Values[6]; set => Values[6] = value; }
        public double DistanceOrder { get => Values[7]; set => Values[7] = value; }
        protected double OriginPointX { get => Values[8]; set => Values[8] = value; }
        protected double OriginPointY { get => Values[9]; set => Values[9] = value; }
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
        public double XThickness { get => Values[10]; set => Values[10] = value; }
        public double YThickness { get => Values[11]; set => Values[11] = value; }
        // 回転
        protected double RotOriginPointX { get => Values[12]; set => Values[12] = value; }
        protected double RotOriginPointY { get => Values[13]; set => Values[13] = value; }
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
        public double RotAngle { get => Values[14]; set => Values[14] = value; }

        public LinearElasticPMLMaterial() : base()
        {
            //-----------------------------
            // baseクラスの値
            int baseLen = Values.Length;
            System.Diagnostics.Debug.Assert(baseLen == 5);
            double[] baseValues = Values;
            int intBaseLen = IntValues.Length;
            System.Diagnostics.Debug.Assert(intBaseLen == 1);
            //int[] baseIntValues = IntValues;

            // このクラスの値の数に変更
            int len = 15;
            Values = new double[len];
            for (int i = 0; i < len; i++)
            {
                Values[i] = 0.0;
            }
            // baseクラスの値をコピー
            baseValues.CopyTo(Values, 0);
            //
            //intValuesは同じ
            //-----------------------------

            // PML
            Reflection0 = 1.0e-8;
            ScalingFactor0 = 1.0;
            DistanceOrder = 2.0;
            OriginPointX = 0.0;
            OriginPointY = 0.0;
            XThickness = 0.0;
            YThickness = 0.0;
            // 回転
            RotOriginPoint = new OpenTK.Vector2d(0.0, 0.0);
        }

        public bool IsXDirection()
        {
            return Math.Abs(XThickness) >= Constants.PrecisionLowerLimit;
        }

        public bool IsYDirection()
        {
            return Math.Abs(YThickness) >= Constants.PrecisionLowerLimit;
        }

        private void CalcScalingDampingFactor(
            double thickness,
            double rho, double lambda, double mu,
            double pos, double origin,
            out double A, out double D)
        {
            double order = DistanceOrder;
            //
            double s = Math.Abs(pos - origin) / thickness;
            System.Diagnostics.Debug.Assert(s <= 1.0);

            // damping factor
            //-----------------------------------------------
            // Vp
            double vp = Math.Sqrt((lambda + 2.0 * mu) / rho); 
            // D0
            System.Diagnostics.Debug.Assert(thickness >= 0);
            double D0 =
                (vp * (order + 1) / (2.0 * thickness)) *
                Math.Log(1.0 / Reflection0);

            // D
            D = D0 * Math.Pow(s, order);

            System.Diagnostics.Debug.Assert(D > 0.0); // 正の値
            //-----------------------------------------------

            // scaling factor
            //-----------------------------------------------
            double A0 = ScalingFactor0;
            A = 1.0 + (A0 - 1.0) * Math.Pow(s, order);
            //-----------------------------------------------
        }

        public void CalcScalingDampingFactorX(
            OpenTK.Vector2d pt, out double AX, out double DX)
        {
            AX = 1.0;
            DX = 0.0;
            if (!IsXDirection())
            {
                return;
            }

            double[] pmlOriginCoord = { OriginPointX, OriginPointY };
            double[] rotOriginCoord = { RotOriginPointX, RotOriginPointY };
            pmlOriginCoord = CadUtils2D.GetRotCoord2D(pmlOriginCoord, RotAngle, rotOriginCoord);
            double thickness = XThickness;
            double rho = MassDensity;
            double lambda = LameLambda;
            double mu = LameMu;
            double pos = pt.X;
            double origin = pmlOriginCoord[0]; // X
            CalcScalingDampingFactor(
                thickness,
                rho, lambda, mu,
                pos, origin,
                out AX, out DX);
        }

        public void CalcScalingDampingFactorY(
            OpenTK.Vector2d pt, out double AY, out double DY)
        {
            AY = 1.0;
            DY = 0.0;
            if (!IsYDirection())
            {
                return;
            }

            double[] pmlOriginCoord = { OriginPointX, OriginPointY };
            double[] rotOriginCoord = { RotOriginPointX, RotOriginPointY };
            pmlOriginCoord = CadUtils2D.GetRotCoord2D(pmlOriginCoord, RotAngle, rotOriginCoord);
            double thickness = YThickness;
            double rho = MassDensity;
            double lambda = LameLambda;
            double mu = LameMu;
            double pos = pt.Y;
            double origin = pmlOriginCoord[1]; // Y
            CalcScalingDampingFactor(
                thickness,
                rho, lambda, mu,
                pos, origin,
                out AY, out DY);
        }

        // Time Domain
        private void CalcScalingDampingFactorForTD(
            double thickness,
            double rho, double lambda, double mu,
            double pos, double origin,
            double dt,
            out double A, out double D,
            out double c1PXi, out double c2PXi,
            out double c1WXi)
        {
            CalcScalingDampingFactor(
                thickness,
                rho, lambda, mu,
                pos, origin,
                out A, out D);

            // c1ψξ,c2ψξ
            c1PXi = 1.0 - Math.Exp(-1.0 * D * dt);
            c2PXi = Math.Exp(-1.0 * D * dt);

            // c1wξ
            c1WXi = D * dt;
        }

        // Time Domain
        public void CalcScalingDampingFactorXForTD(
            OpenTK.Vector2d pt,
            double dt,
            out double AX, out double DX,
            out double c1PX, out double c2PX,
            out double c1WX)
        {
            AX = 1.0;
            DX = 0.0;
            c1PX = 0.0;
            c2PX = 0.0;
            c1WX = 0.0;
            if (!IsXDirection())
            {
                return;
            }

            double[] pmlOriginCoord = { OriginPointX, OriginPointY };
            double[] rotOriginCoord = { RotOriginPointX, RotOriginPointY };
            pmlOriginCoord = CadUtils2D.GetRotCoord2D(pmlOriginCoord, RotAngle, rotOriginCoord);
            double thickness = XThickness;
            double rho = MassDensity;
            double lambda = LameLambda;
            double mu = LameMu;
            double pos = pt.X;
            double origin = pmlOriginCoord[0]; // X
            CalcScalingDampingFactorForTD(
                thickness,
                rho, lambda, mu,
                pos, origin,
                dt,
                out AX, out DX,
                out c1PX, out c2PX,
                out c1WX);
        }

        // Time Domain
        public void CalcScalingDampingFactorYForTD(
            OpenTK.Vector2d pt,
            double dt,
            out double AY, out double DY,
            out double c1PY, out double c2PY,
            out double c1WY)
        {
            AY = 1.0;
            DY = 0.0;
            c1PY = 0.0;
            c2PY = 0.0;
            c1WY = 0.0;
            if (!IsYDirection())
            {
                return;
            }

            double[] pmlOriginCoord = { OriginPointX, OriginPointY };
            double[] rotOriginCoord = { RotOriginPointX, RotOriginPointY };
            pmlOriginCoord = CadUtils2D.GetRotCoord2D(pmlOriginCoord, RotAngle, rotOriginCoord);
            double thickness = YThickness;
            double rho = MassDensity;
            double lambda = LameLambda;
            double mu = LameMu;
            double pos = pt.Y;
            double origin = pmlOriginCoord[1]; // Y
            CalcScalingDampingFactorForTD(
                thickness,
                rho, lambda, mu,
                pos, origin,
                dt,
                out AY, out DY,
                out c1PY, out c2PY,
                out c1WY);
        }
    }
}
