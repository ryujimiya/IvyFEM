using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class ElasticBaseMaterial : Material
    {
        public double MassDensity { get => Values[0]; set => Values[0] = value; }
        public double GravityX { get => Values[1]; set => Values[1] = value; }
        public double GravityY { get => Values[2]; set => Values[2] = value; }
        public double GravityZ { get => Values[3]; set => Values[3] = value; }
        public double Young { get => Values[4]; set => Values[4] = value; }
        public double Poisson { get => Values[5]; set => Values[5] = value; }
        public double LameLambda {
            get
            {
                double lambda = Young * Poisson / ((1.0 + Poisson) * (1 - 2.0 * Poisson));
                if (IsPlainStressLame)
                {
                    double lambdaPlaneStress = 2.0 * lambda * LameMu / (lambda + 2.0 * LameMu);
                    lambda = lambdaPlaneStress;
                }
                return lambda;
            }
        }
        public double LameMu => Young / (2.0 * (1.0 + Poisson));
        //
        public bool IsPlainStressLame { get => (IntValues[0] == 1); set => IntValues[0] = value? 1 : 0; } 

        public ElasticBaseMaterial()
        {
            int len = 6;
            Values = new double[len];

            int intLen = 1;
            IntValues = new int[intLen];

            MassDensity = 1.0;
            GravityX = 0.0;
            GravityY = 0.0;
            GravityZ = 0.0;
            Young = 0.0;
            Poisson = 0.0;
            IsPlainStressLame = false;
        }

        public ElasticBaseMaterial(ElasticBaseMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }
    }
}
