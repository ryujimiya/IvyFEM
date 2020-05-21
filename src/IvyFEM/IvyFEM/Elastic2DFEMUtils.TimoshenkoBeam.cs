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
        public static IvyFEM.Lapack.DoubleMatrix CalcTimoshenkoBeamKl(
            uint dElemNodeCnt, uint rElemNodeCnt, int dDof, int rDof,
            LineFE dLineFE, LineFE rLineFE,
            double E, double G, double kappa, double Ae, double Iz)
        {
            // Timoshenko Beam
            int localBeamMaxNodeCnt = (int)Math.Max(dElemNodeCnt, rElemNodeCnt);
            var localKeBeam = new IvyFEM.Lapack.DoubleMatrix(
                localBeamMaxNodeCnt * (dDof + rDof), localBeamMaxNodeCnt * (dDof + rDof));
            IntegrationPoints ip;
            if (dLineFE.Order == 1 && rLineFE.Order == 1)
            {
                // 低減積分
                ip = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point1);
                System.Diagnostics.Debug.Assert(ip.Ls.Length == 1);
            }
            else
            {
                ip = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);
                System.Diagnostics.Debug.Assert(ip.Ls.Length == 5);
            }
            int beamDof = 2;
            int localBeamROffset = 1;
            System.Diagnostics.Debug.Assert(beamDof == (dDof + rDof));
            System.Diagnostics.Debug.Assert(localBeamROffset == dDof);
            for (int ipPt = 0; ipPt < ip.PointCount; ipPt++)
            {
                double[] L = ip.Ls[ipPt];
                double[] d2N = dLineFE.CalcN(L);
                double[][] d2Nu = dLineFE.CalcNu(L);
                double[] d2Nx = d2Nu[0];
                double[] rN = rLineFE.CalcN(L);
                double[][] rNu = rLineFE.CalcNu(L);
                double[] rNx = rNu[0];
                double lineLen = dLineFE.GetLineLength();
                double weight = ip.Weights[ipPt];
                double detJWeight = (lineLen / 2.0) * weight;

                // displacement
                for (int row = 0; row < dElemNodeCnt; row++)
                {
                    // displacement
                    for (int col = 0; col < dElemNodeCnt; col++)
                    {
                        double kValue = detJWeight * kappa * G * Ae * d2Nx[row] * d2Nx[col];
                        localKeBeam[row * beamDof, col * beamDof] += kValue;
                    }
                    // rotation
                    for (int col = 0; col < rElemNodeCnt; col++)
                    {
                        double kValue = -1.0 * detJWeight * kappa * G * Ae * d2Nx[row] * rN[col];
                        localKeBeam[row * beamDof, col * beamDof + localBeamROffset] += kValue;
                    }
                }
                // rotation
                for (int row = 0; row < rElemNodeCnt; row++)
                {
                    // displacement
                    for (int col = 0; col < dElemNodeCnt; col++)
                    {
                        double kValue = -1.0 * detJWeight * kappa * G * Ae * rN[row] * d2Nx[col];
                        localKeBeam[row * beamDof + localBeamROffset, col * beamDof] += kValue;
                    }
                    // rotation
                    for (int col = 0; col < rElemNodeCnt; col++)
                    {
                        double kValue1 = detJWeight * kappa * G * Ae * rN[row] * rN[col];
                        double kValue2 = detJWeight * E * Iz * rNx[row] * rNx[col];
                        localKeBeam[row * beamDof + localBeamROffset, col * beamDof + localBeamROffset] +=
                            kValue1 + kValue2;
                    }
                }
            }
            return localKeBeam;
        }

        //////////////////////////////////////////
        // mass matrix
        public static IvyFEM.Lapack.DoubleMatrix CalcTimoshenkoBeamMl(
            uint dElemNodeCnt, uint rElemNodeCnt, int dDof, int rDof,
            LineFE dLineFE, LineFE rLineFE,
            double rho, double Ae, double Iz)
        {
            int localBeamMaxNodeCnt = (int)Math.Max(dElemNodeCnt, rElemNodeCnt);
            var localKeBeam = new IvyFEM.Lapack.DoubleMatrix(
                localBeamMaxNodeCnt * (dDof + rDof), localBeamMaxNodeCnt * (dDof + rDof));
            var localMeBeam = new IvyFEM.Lapack.DoubleMatrix(
                localBeamMaxNodeCnt * (dDof + rDof), localBeamMaxNodeCnt * (dDof + rDof));
            int beamDof = 2;
            int localBeamROffset = 1;
            System.Diagnostics.Debug.Assert(beamDof == (dDof + rDof));
            IntegrationPoints ipM = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);
            System.Diagnostics.Debug.Assert(ipM.Ls.Length == 5);
            for (int ipPt = 0; ipPt < ipM.PointCount; ipPt++)
            {
                double[] L = ipM.Ls[ipPt];
                double[] d2N = dLineFE.CalcN(L);
                double[][] d2Nu = dLineFE.CalcNu(L);
                double[] d2Nx = d2Nu[0];
                double[] rN = rLineFE.CalcN(L);
                double[][] rNu = rLineFE.CalcNu(L);
                double[] rNx = rNu[0];
                double lineLen = dLineFE.GetLineLength();
                double weight = ipM.Weights[ipPt];
                double detJWeight = (lineLen / 2.0) * weight;

                // displacement
                for (int row = 0; row < dElemNodeCnt; row++)
                {
                    // displacement
                    for (int col = 0; col < dElemNodeCnt; col++)
                    {
                        double mValue = detJWeight * rho * Ae * d2N[row] * d2N[col];
                        localMeBeam[row * beamDof, col * beamDof] += mValue;
                    }
                    // rotation
                    for (int col = 0; col < rElemNodeCnt; col++)
                    {
                        double mValue = 0.0;
                        localMeBeam[row * beamDof, col * beamDof + localBeamROffset] += mValue;
                    }
                }
                // rotation
                for (int row = 0; row < rElemNodeCnt; row++)
                {
                    // displacement
                    for (int col = 0; col < dElemNodeCnt; col++)
                    {
                        double mValue = 0.0;
                        localMeBeam[row * beamDof + localBeamROffset, col * beamDof] += mValue;
                    }
                    // rotation
                    for (int col = 0; col < rElemNodeCnt; col++)
                    {
                        double mValue = detJWeight * rho * Iz * rN[row] * rN[col];
                        localMeBeam[row * beamDof + localBeamROffset, col * beamDof + localBeamROffset] +=
                            mValue;
                    }
                }
            }
            return localMeBeam;
        }
    }
}
