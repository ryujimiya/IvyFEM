using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DFEMUtils
    {
        //////////////////////////////////////////
        // stiffness matrix
        public static void CalcSimpleCorotationalFrameKl(
            double l0, double E, double Ae, double Iz,
            double barU, double barT1, double barT2,
            out double[] fl, out IvyFEM.Lapack.DoubleMatrix kl)
        {
            kl = new IvyFEM.Lapack.DoubleMatrix(3, 3);
            double n = (E * Ae / l0) * barU;
            double m1 = (2.0 * E * Iz / l0) * (2.0 * barT1 + barT2);
            double m2 = (2.0 * E * Iz / l0) * (barT1 + 2.0 * barT2);
            kl[0, 0] = E * Ae / l0;
            kl[0, 1] = 0.0;
            kl[0, 2] = 0.0;
            kl[1, 0] = kl[0, 1];
            kl[1, 1] = 4.0 * E * Iz / l0;
            kl[1, 2] = 2.0 * E * Iz / l0;
            kl[2, 0] = kl[0, 2];
            kl[2, 1] = kl[1, 2];
            kl[2, 2] = 4.0 * E * Iz / l0;

            fl = new double[3];
            fl[0] = n;
            fl[1] = m1;
            fl[2] = m2;
        }

        public static void CalcBernoulliCorotationalFrameKl(
            double l0, double E, double Ae, double Iz,
            double barU, double barT1, double barT2,
            out double[] fl, out IvyFEM.Lapack.DoubleMatrix kl)
        {
            kl = new IvyFEM.Lapack.DoubleMatrix(3, 3);
            double n = 0.0;
            double m1 = 0.0;
            double m2 = 0.0;
            IntegrationPoints ip = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);
            System.Diagnostics.Debug.Assert(ip.Ls.Length == 5);
            for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
            {
                double[] L = ip.Ls[ipPt];
                double lineLen = l0;
                double weight = ip.Weights[ipPt];
                double detJWeight = (lineLen / 2.0) * weight;

                double x = l0 * L[1];

                double sigma = E * Ae * barU / l0;
                double sigmaY = E * Iz * (
                    (4.0 / l0 - (6.0 / (l0 * l0)) * x) * barT1 +
                    (2.0 / l0 - (6.0 / (l0 * l0)) * x) * barT2
                    );

                n += detJWeight * (1.0 / l0) * sigma;
                m1 += detJWeight * (4.0 / l0 - (6.0 / (l0 * l0)) * x) * sigmaY;
                m2 += detJWeight * (2.0 / l0 - (6.0 / (l0 * l0)) * x) * sigmaY;

                kl[0, 0] += detJWeight * (E * Ae / (l0 * l0));
                kl[0, 1] += 0.0;
                kl[0, 2] += 0.0;
                kl[1, 1] += detJWeight * E * Iz *
                    (4.0 / l0 - (6.0 / (l0 * l0)) * x) * (4.0 / l0 - (6.0 / (l0 * l0)) * x);
                kl[1, 2] += detJWeight * E * Iz *
                    (2.0 / l0 - (6.0 / (l0 * l0)) * x) * (4.0 / l0 - (6.0 / (l0 * l0)) * x);
                kl[2, 2] += detJWeight * E * Iz *
                    (2.0 / l0 - (6.0 / (l0 * l0)) * x) * (2.0 / l0 - (6.0 / (l0 * l0)) * x);
            }
            kl[1, 0] = kl[0, 1];
            kl[2, 0] = kl[0, 2];
            kl[2, 1] = kl[1, 2];

            fl = new double[3];
            fl[0] = n;
            fl[1] = m1;
            fl[2] = m2;
        }

        public static void CalcShallowArchCorotationalFrameKl(
            double l0, double E, double Ae, double Iz,
            double barU, double barT1, double barT2,
            out double[] fl, out IvyFEM.Lapack.DoubleMatrix kl)
        {
            kl = new IvyFEM.Lapack.DoubleMatrix(3, 3);
            double n = 0.0;
            double m1 = 0.0;
            double m2 = 0.0;
            IntegrationPoints ip = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);
            System.Diagnostics.Debug.Assert(ip.Ls.Length == 5);
            for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
            {
                double[] L = ip.Ls[ipPt];
                double lineLen = l0;
                double weight = ip.Weights[ipPt];
                double detJWeight = (lineLen / 2.0) * weight;

                double x = l0 * L[1];

                double sigma = E * Ae * (
                    barU / l0 +
                    (1.0 / 15.0) * barT1 * barT1 -
                    (1.0 / 30.0) * barT1 * barT2 +
                    (1.0 / 15.0) * barT2 * barT2
                    );
                double sigmaY = E * Iz * (
                    (4.0 / l0 - (6.0 / (l0 * l0)) * x) * barT1 +
                    (2.0 / l0 - (6.0 / (l0 * l0)) * x) * barT2
                    );

                n += detJWeight * (1.0 / l0) * sigma;
                m1 += detJWeight * ((2.0 / 15.0) * barT1 - (1.0 / 30.0) * barT2) * sigma +
                    detJWeight * (4.0 / l0 - (6.0 / (l0 * l0)) * x) * sigmaY;
                m2 += detJWeight * ((2.0 / 15.0) * barT2 - (1.0 / 30.0) * barT1) * sigma +
                    detJWeight * (2.0 / l0 - (6.0 / (l0 * l0)) * x) * sigmaY;

                kl[0, 0] += detJWeight * (E * Ae / (l0 * l0));
                kl[0, 1] += detJWeight * (E * Ae / l0) * ((2.0 / 15.0) * barT1 - (1.0 / 30.0) * barT2);
                kl[0, 2] += detJWeight * (E * Ae / l0) * ((2.0 / 15.0) * barT2 - (1.0 / 30.0) * barT1);
                kl[1, 1] += detJWeight * E * Ae *
                    ((2.0 / 15.0) * barT1 - (1.0 / 30.0) * barT2) * ((2.0 / 15.0) * barT1 - (1.0 / 30.0) * barT2) +
                    detJWeight * (2.0 / 15.0) * sigma +
                    detJWeight * E * Iz *
                    (4.0 / l0 - (6.0 / (l0 * l0)) * x) * (4.0 / l0 - (6.0 / (l0 * l0)) * x);
                kl[1, 2] += detJWeight * E * Ae *
                    ((2.0 / 15.0) * barT1 - (1.0 / 30.0) * barT2) * ((2.0 / 15.0) * barT2 - (1.0 / 30.0) * barT1) -
                    detJWeight * (1.0 / 30.0) * sigma +
                    detJWeight * E * Iz *
                    (2.0 / l0 - (6.0 / (l0 * l0)) * x) * (4.0 / l0 - (6.0 / (l0 * l0)) * x);
                kl[2, 2] += detJWeight * E * Ae *
                    ((2.0 / 15.0) * barT2 - (1.0 / 30.0) * barT1) * ((2.0 / 15.0) * barT2 - (1.0 / 30.0) * barT1) +
                    detJWeight * (2.0 / 15.0) * sigma +
                    detJWeight * E * Iz *
                    (2.0 / l0 - (6.0 / (l0 * l0)) * x) * (2.0 / l0 - (6.0 / (l0 * l0)) * x);
            }
            kl[1, 0] = kl[0, 1];
            kl[2, 0] = kl[0, 2];
            kl[2, 1] = kl[1, 2];

            fl = new double[3];
            fl[0] = n;
            fl[1] = m1;
            fl[2] = m2;
        }

        //////////////////////////////////////////
        // mass matrix
        public static IvyFEM.Lapack.DoubleMatrix CalcCorotationalFrameMl(
            double l0, double rho, double Ae, double Iz)
        {
            var ml1 = new IvyFEM.Lapack.DoubleMatrix(6, 6);
            double c1 = rho * Ae * l0 / 420.0;
            ml1[0, 0] = c1 * 140.0;
            ml1[0, 1] = c1 * 0.0;
            ml1[0, 2] = c1 * 0.0;
            ml1[0, 3] = c1 * 70.0;
            ml1[0, 4] = c1 * 0.0;
            ml1[0, 5] = c1 * 0.0;
            ml1[1, 1] = c1 * 156.0;
            ml1[1, 2] = c1 * 22.0 * l0;
            ml1[1, 3] = c1 * 0.0;
            ml1[1, 4] = c1 * 54.0;
            ml1[1, 5] = c1 * (-13.0 * l0);
            ml1[2, 2] = c1 * 4.0 * l0 * l0;
            ml1[2, 3] = c1 * 0.0;
            ml1[2, 4] = c1 * 13.0 * l0;
            ml1[2, 5] = c1 * (-3.0 * l0 * l0);
            ml1[3, 3] = c1 * 140.0;
            ml1[3, 4] = c1 * 0.0;
            ml1[3, 5] = c1 * 0.0;
            ml1[4, 4] = c1 * 156.0;
            ml1[4, 5] = c1 * (-22.0 * l0);
            ml1[5, 5] = c1 * 4.0 * l0 * l0;
            for (int i = 0; i < 6; i++)
            {
                for (int j = i + 1; j < 6; j++)
                {
                    ml1[j, i] = ml1[i, j];
                }
            }
            var ml2 = new IvyFEM.Lapack.DoubleMatrix(6, 6);
            double c2 = rho * Iz / (30.0 * l0);
            ml2[0, 0] = c2 * 0.0;
            ml2[0, 1] = c2 * 0.0;
            ml2[0, 2] = c2 * 0.0;
            ml2[0, 3] = c2 * 0.0;
            ml2[0, 4] = c2 * 0.0;
            ml2[0, 5] = c2 * 0.0;
            ml2[1, 1] = c2 * 36.0;
            ml2[1, 2] = c2 * 3.0 * l0;
            ml2[1, 3] = c2 * 0.0;
            ml2[1, 4] = c2 * (-36.0);
            ml2[1, 5] = c2 * 3.0 * l0;
            ml2[2, 2] = c2 * 4.0 * l0 * l0;
            ml2[2, 3] = c2 * 0.0;
            ml2[2, 4] = c2 * (-3.0 * l0);
            ml2[2, 5] = c2 * (-l0 * l0);
            ml2[3, 3] = c2 * 0.0;
            ml2[3, 4] = c2 * 0.0;
            ml2[3, 5] = c2 * 0.0;
            ml2[4, 4] = c2 * 36.0;
            ml2[4, 5] = c2 * (-3.0 * l0);
            ml2[5, 5] = c2 * 4.0 * l0 * l0;
            for (int i = 0; i < 6; i++)
            {
                for (int j = i + 1; j < 6; j++)
                {
                    ml2[j, i] = ml2[i, j];
                }
            }

            var ml = new IvyFEM.Lapack.DoubleMatrix(6, 6);
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    ml[i, j] = ml1[i, j] + ml2[i, j];
                }
            }
            return ml;
        }
    }
}
