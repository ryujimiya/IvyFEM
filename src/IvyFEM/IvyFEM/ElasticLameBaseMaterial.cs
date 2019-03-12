using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class ElasticLameBaseMaterial : Material
    {
        public double MassDensity { get => Values[0]; set => Values[0] = value; }
        public double GravityX { get => Values[1]; set => Values[1] = value; }
        public double GravityY { get => Values[2]; set => Values[2] = value; }
        public double LameLambda { get => Values[3]; set => Values[3] = value; }
        public double LameMu { get => Values[4]; set => Values[4] = value; }

        public ElasticLameBaseMaterial()
        {
            int len = 5;
            Values = new double[len];

            MassDensity = 1.0;
            GravityX = 0.0;
            GravityY = 0.0;
            LameLambda = 0.0;
            LameMu = 1.0;
        }

        public void SetYoungPoisson(double young, double poisson)
        {
            LameLambda = young * poisson / ((1.0 + poisson) * (1 - 2.0 * poisson));
            LameMu = young / (2.0 * (1.0 + poisson));
        }

        public void GetYoungPoisson(out double young, out double poisson)
        {
            poisson = LameLambda * 0.5 / (LameLambda + LameMu);
            young = 2 * LameMu * (1 + poisson);
        }
    }
}
