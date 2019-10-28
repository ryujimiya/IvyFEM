﻿using System;
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
        protected double OriginPointX { get => Values[6]; set => Values[6] = value; }
        protected double OriginPointY { get => Values[7]; set => Values[7] = value; }
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
        public double XThickness { get => Values[8]; set => Values[8] = value; }
        public double YThickness { get => Values[9]; set => Values[9] = value; }
        public double Reflection0 { get => Values[10]; set => Values[10] = value; }

        public bool IsTMMode { get => IntValues[0] == 1; set => IntValues[0] = value ? 1 : 0; }

        public DielectricPMLMaterial() : base()
        {
            int len = 11;
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

            int intLen = 1;
            IntValues = new int[intLen];
            for (int i = 0; i < intLen; i++)
            {
                IntValues[i] = 0;
            }
            IsTMMode = false;
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

        private double CalcSigma(
            double thickness,
            double ep, double mu,
            double pos, double origin)
        {
            // σmax
            System.Diagnostics.Debug.Assert(thickness >= 0);
            double sigmaMax = (3.0 / (2.0 * thickness)) *
                Math.Sqrt(Constants.Ep0 / Constants.Mu0) *
                Math.Sqrt(ep / mu) *
                Math.Log(1.0 / Reflection0);

            // σ
            double s = Math.Abs(pos - origin) / thickness;
            System.Diagnostics.Debug.Assert(s <= 1.0);

            double sigma = sigmaMax * s * s;

            return sigma;
        }

        public double CalcSigmaX(OpenTK.Vector2d pt)
        {
            double sigmaX = 0.0;
            if (!IsXDirection())
            {
                return sigmaX;
            }

            double thickness = XThickness;
            double ep = 0.0;
            double mu = 0.0;
            if (IsTMMode)
            {
                // TMモードの場合、MuにEpが格納されている
                ep = Muxx;
                mu = Epxx;
            }
            else
            {
                ep = Epxx;
                mu = Muxx;
            }
            double pos = pt.X;
            double origin = OriginPointX;
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

            double thickness = YThickness;
            double ep = 0.0;
            double mu = 0.0;
            if (IsTMMode)
            {
                // TMモードの場合、MuにEpが格納されている
                ep = Muyy;
                mu = Epyy;
            }
            else
            {
                ep = Epyy;
                mu = Muyy;
            }
            double pos = pt.Y;
            double origin = OriginPointY;
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

            double thickness = XThickness;
            double ep = 0.0;
            double mu = 0.0;
            if (IsTMMode)
            {
                // TMモードの場合、MuにEpが格納されている
                ep = Muxx;
                mu = Epxx;
            }
            else
            {
                ep = Epxx;
                mu = Muxx;
            }
            double pos = pt.X;
            double origin = OriginPointX;
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

            double thickness = YThickness;
            double ep = 0.0;
            double mu = 0.0;
            if (IsTMMode)
            {
                // TMモードの場合、MuにEpが格納されている
                ep = Muyy;
                mu = Epyy;
            }
            else
            {
                ep = Epyy;
                mu = Muyy;
            }
            double pos = pt.Y;
            double origin = OriginPointY;
            CalcSigmaForTD(
                thickness,
                ep, mu,
                pos, origin,
                dt,
                out sigmaY,
                out c1PY, out c2PY,
                out c1VY);
        }
    }
}