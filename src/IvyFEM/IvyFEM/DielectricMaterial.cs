using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class DielectricMaterial : Material
    {
        public double Epxx { get => Values[0]; set => Values[0] = value; }
        public double Epyy { get => Values[1]; set => Values[1] = value; }
        public double Epzz { get => Values[2]; set => Values[2] = value; }
        public double Muxx { get => Values[3]; set => Values[3] = value; }
        public double Muyy { get => Values[4]; set => Values[4] = value; }
        public double Muzz { get => Values[5]; set => Values[5] = value; }

        public double ImagEpxx { get => Values[6]; set => Values[6] = value; }
        public double ImagEpyy { get => Values[7]; set => Values[7] = value; }
        public double ImagEpzz { get => Values[8]; set => Values[8] = value; }
        public double ImagMuxx { get => Values[9]; set => Values[9] = value; }
        public double ImagMuyy { get => Values[10]; set => Values[10] = value; }
        public double ImagMuzz { get => Values[11]; set => Values[11] = value; }

        // helper
        public System.Numerics.Complex ComplexEpxx
        {
            get => new System.Numerics.Complex(Epxx, ImagEpxx);
            set
            {
                Epxx = value.Real;
                ImagEpxx = value.Imaginary;
            }
        }
        public System.Numerics.Complex ComplexEpyy
        {
            get => new System.Numerics.Complex(Epyy, ImagEpyy);
            set
            {
                Epyy = value.Real;
                ImagEpyy = value.Imaginary;
            }
        }
        public System.Numerics.Complex ComplexEpzz
        {
            get => new System.Numerics.Complex(Epzz, ImagEpzz);
            set
            {
                Epzz = value.Real;
                ImagEpzz = value.Imaginary;
            }
        }
        public System.Numerics.Complex ComplexMuxx
        {
            get => new System.Numerics.Complex(Muxx, ImagMuxx);
            set
            {
                Muxx = value.Real;
                ImagMuxx = value.Imaginary;
            }
        }
        public System.Numerics.Complex ComplexMuyy
        {
            get => new System.Numerics.Complex(Muyy, ImagMuyy);
            set
            {
                Muyy = value.Real;
                ImagMuyy = value.Imaginary;
            }
        }
        public System.Numerics.Complex ComplexMuzz
        {
            get => new System.Numerics.Complex(Muzz, ImagMuzz);
            set
            {
                Muzz = value.Real;
                ImagMuzz = value.Imaginary;
            }
        }

        public DielectricMaterial() : base()
        {
            int len = 12;
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

            ImagEpxx = 0.0;
            ImagEpyy = 0.0;
            ImagEpzz = 0.0;
            ImagMuxx = 0.0;
            ImagMuyy = 0.0;
            ImagMuzz = 0.0;
        }

        public DielectricMaterial(DielectricMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }
    }
}
