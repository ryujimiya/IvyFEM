using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DFEMUtils
    {
        // Linear / Saint Venant
        public static void SetStressValue(
            FEWorld world,
            uint displacementValueId, uint stressValueId, uint equivStressValueId)
        {
            System.Diagnostics.Debug.Assert(world.IsFieldValueId(displacementValueId));
            FieldValue uFV = world.GetFieldValue(displacementValueId);
            uint uQuantityId = uFV.QuantityId;

            FieldValue sigmaFV = null;
            if (stressValueId != 0)
            {
                System.Diagnostics.Debug.Assert(world.IsFieldValueId(stressValueId));
                sigmaFV = world.GetFieldValue(stressValueId);
                System.Diagnostics.Debug.Assert(sigmaFV.Type == FieldValueType.SymmetricTensor2);
                System.Diagnostics.Debug.Assert(sigmaFV.Dof == 3);
            }
            FieldValue eqSigmaFV = null;
            if (equivStressValueId != 0)
            {
                System.Diagnostics.Debug.Assert(world.IsFieldValueId(equivStressValueId));
                eqSigmaFV = world.GetFieldValue(equivStressValueId);
                System.Diagnostics.Debug.Assert(eqSigmaFV.Type == FieldValueType.Scalar);
                System.Diagnostics.Debug.Assert(eqSigmaFV.Dof == 1);
            }

            IList<uint> feIds = world.GetTriangleFEIds(uQuantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = world.GetTriangleFE(uQuantityId, feId);
                int[] coIds = triFE.NodeCoordIds;
                Material ma = world.GetMaterial(triFE.MaterialId);
                double lambda = 0;
                double mu = 0;
                if (ma is LinearElasticMaterial)
                {
                    var ma1 = ma as LinearElasticMaterial;
                    lambda = ma1.LameLambda;
                    mu = ma1.LameMu;
                }
                else if (ma is StVenantHyperelasticMaterial)
                {
                    var ma1 = ma as StVenantHyperelasticMaterial;
                    lambda = ma1.LameLambda;
                    mu = ma1.LameMu;
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                    throw new NotImplementedException();
                }

                var ip = TriangleFE.GetIntegrationPoints(TriangleIntegrationPointCount.Point1);
                System.Diagnostics.Debug.Assert(ip.PointCount == 1);
                double[] L = ip.Ls[0];
                double[][] Nu = triFE.CalcNu(L);
                double[] Nx = Nu[0];
                double[] Ny = Nu[1];
                double[,] uu = new double[2, 2];
                for (int iNode = 0; iNode < coIds.Length; iNode++)
                {
                    int coId = coIds[iNode];
                    double[] u = uFV.GetDoubleValue(coId, FieldDerivativeType.Value);
                    uu[0, 0] += u[0] * Nx[iNode];
                    uu[0, 1] += u[0] * Ny[iNode];
                    uu[1, 0] += u[1] * Nx[iNode];
                    uu[1, 1] += u[1] * Ny[iNode];
                }

                //ε strain
                double[,] eps = new double[2, 2];
                if (ma is LinearElasticMaterial)
                {
                    eps[0, 0] = (1.0 / 2.0) * (uu[0, 0] + uu[0, 0]);
                    eps[0, 1] = (1.0 / 2.0) * (uu[0, 1] + uu[1, 0]);
                    eps[1, 0] = (1.0 / 2.0) * (uu[1, 0] + uu[0, 1]);
                    eps[1, 1] = (1.0 / 2.0) * (uu[1, 1] + uu[1, 1]);
                }
                else if (ma is StVenantHyperelasticMaterial)
                {
                    eps[0, 0] = (1.0 / 2.0) * (uu[0, 0] + uu[0, 0] + uu[0, 0] * uu[0, 0] + uu[1, 0] * uu[1, 0]);
                    eps[0, 1] = (1.0 / 2.0) * (uu[0, 1] + uu[1, 0] + uu[0, 0] * uu[0, 1] + uu[1, 1] * uu[1, 0]);
                    eps[1, 0] = (1.0 / 2.0) * (uu[1, 0] + uu[0, 1] + uu[0, 1] * uu[0, 0] + uu[1, 0] * uu[1, 1]);
                    eps[1, 1] = (1.0 / 2.0) * (uu[1, 1] + uu[1, 1] + uu[0, 1] * uu[0, 1] + uu[1, 1] * uu[1, 1]);
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }

                // σ stress
                double[,] sigma = new double[2, 2];
                {
                    sigma[0, 0] = mu * eps[0, 0];
                    sigma[0, 1] = mu * eps[0, 1];
                    sigma[1, 0] = mu * eps[1, 0];
                    sigma[1, 1] = mu * eps[1, 1];
                    double tmp = lambda * (eps[0, 0] + eps[1, 1]);
                    sigma[0, 0] += tmp;
                    sigma[1, 1] += tmp;
                }

                double misesStress = Math.Sqrt(
                    (1.0 / 2.0) * (
                    sigma[0, 0] * sigma[0, 0] + sigma[1, 1] * sigma[1, 1] +
                    (sigma[1, 1] - sigma[0, 0]) * (sigma[1, 1] - sigma[0, 0])
                    ) +
                    3.0 * sigma[0, 1] * sigma[0, 1]);
                if (stressValueId != 0)
                {
                    double[] Sigma = sigmaFV.GetDoubleValues(FieldDerivativeType.Value);
                    uint dof = sigmaFV.Dof;
                    Sigma[(feId - 1) * dof + 0] = sigma[0, 0]; // σxx
                    Sigma[(feId - 1) * dof + 1] = sigma[1, 1]; // σyy
                    Sigma[(feId - 1) * dof + 2] = sigma[0, 1]; // τxy
                }
                if (equivStressValueId != 0)
                {
                    double[] EqSigma = eqSigmaFV.GetDoubleValues(FieldDerivativeType.Value);
                    EqSigma[feId - 1] = misesStress;
                }
            }
        }
    }
}
