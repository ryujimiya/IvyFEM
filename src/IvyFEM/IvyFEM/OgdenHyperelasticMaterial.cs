using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class OgdenHyperelasticMaterial : Material
    {
        public bool IsCompressible { get => IntValues[0] == 1; set => IntValues[0] = (value ? 1 : 0); }
        public int Order { get => IntValues[1]; private set => IntValues[1] = value; }
        private const int AlphaOffset = 4;
        public double MassDensity { get => Values[0]; set => Values[0] = value; }
        public double GravityX { get => Values[1]; set => Values[1] = value; }
        public double GravityY { get => Values[2]; set => Values[2] = value; }
        public double D1 { get => Values[3]; set => Values[3] = value; }
        public double[] Alphas
        {
            get
            {
                if (Order == 0)
                {
                    return null;
                }
                double[] alphas = new double[Order];
                for (int i = 0; i < Order; i++)
                {
                    alphas[i] = Values[AlphaOffset + i];
                }
                return alphas;
            }
        }
        public double[] Mus
        {
            get
            {
                if (Order == 0)
                {
                    return null;
                }
                double[] mus = new double[Order];
                for (int i = 0; i < Order; i++)
                {
                    mus[i] = Values[AlphaOffset + Order + i];
                }
                return mus;
            }
        }

        public OgdenHyperelasticMaterial()
        {
            int len = AlphaOffset;
            Values = new double[len];
            int intLen = 2;
            IntValues = new int[intLen];

            MassDensity = 1.0;
            GravityX = 0.0;
            GravityY = 0.0;
            IsCompressible = false;
            D1 = 1.0;
        }

        public OgdenHyperelasticMaterial(OgdenHyperelasticMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }

        public void SetAlphaMu(int order, double[] alphas, double[] mus)
        {
            Order = order;
            System.Diagnostics.Debug.Assert(alphas.Length == Order);
            System.Diagnostics.Debug.Assert(mus.Length == Order);
            double[] oldValues = new double[Values.Length];
            Values.CopyTo(oldValues, 0);
            Values = new double[AlphaOffset + Order * 2];
            for (int i = 0; i < AlphaOffset; i++)
            {
                Values[i] = oldValues[i];
            }
            for (int i = 0; i < Order; i++)
            {
                Values[AlphaOffset + i] = alphas[i];
                Values[AlphaOffset + Order + i] = mus[i];
            }
        }
    }
}
