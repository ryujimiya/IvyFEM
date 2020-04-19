using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DEigenFEM
    {
        protected void CalcTimoshenkoBeamElementKMForLine(
            uint feId, IvyFEM.Lapack.DoubleMatrix K, IvyFEM.Lapack.DoubleMatrix M)
        {
            {
                // Note: feIdは線要素
                uint quantityId0 = 0;
                LineFE workLineFE = World.GetLineFE(quantityId0, feId);
                uint workMaId = workLineFE.MaterialId;
                if (!World.IsMaterialId(workMaId))
                {
                    return;
                }
                Material workMa0 = World.GetMaterial(workMaId);
                if (!(workMa0 is TimoshenkoBeamMaterial))
                {
                    return;
                }
            }

            uint dQuantityId = 0; // displacement
            uint rQuantityId = 1; // rotation
            System.Diagnostics.Debug.Assert(World.GetDof(dQuantityId) == 1);
            System.Diagnostics.Debug.Assert(World.GetDof(rQuantityId) == 1);
            int dDof = 1; // w
            int rDof = 1; //θ
            int dNodeCnt = (int)World.GetNodeCount(dQuantityId);
            int rNodeCnt = (int)World.GetNodeCount(rQuantityId);
            int offset = dNodeCnt * dDof;

            // Note: feIdは線要素
            LineFE dLineFE = World.GetLineFE(dQuantityId, feId);
            LineFE rLineFE = World.GetLineFE(rQuantityId, feId);
            uint maId = dLineFE.MaterialId;
            if (!World.IsMaterialId(maId))
            {
                return;
            }
            Material ma0 = World.GetMaterial(maId);
            System.Diagnostics.Debug.Assert(ma0 is TimoshenkoBeamMaterial);

            int[] dCoIds = dLineFE.NodeCoordIds;
            uint dElemNodeCnt = dLineFE.NodeCount;
            int[] dNodes = new int[dElemNodeCnt];
            for (int iNode = 0; iNode < dElemNodeCnt; iNode++)
            {
                int coId = dCoIds[iNode];
                int nodeId = World.Coord2Node(dQuantityId, coId);
                dNodes[iNode] = nodeId;
            }
            int[] rCoIds = rLineFE.NodeCoordIds;
            uint rElemNodeCnt = rLineFE.NodeCount;
            int[] rNodes = new int[rElemNodeCnt];
            for (int iNode = 0; iNode < rElemNodeCnt; iNode++)
            {
                int coId = rCoIds[iNode];
                int nodeId = World.Coord2Node(rQuantityId, coId);
                rNodes[iNode] = nodeId;
            }
            int[] vCoIds = dLineFE.VertexCoordIds;
            uint elemVertexCnt = dLineFE.VertexCount;
            double[][] vCoords = new double[elemVertexCnt][];
            for (int iVertex = 0; iVertex < elemVertexCnt; iVertex++)
            {
                int coId = vCoIds[iVertex];
                double[] coord = World.GetCoord(dQuantityId, coId);
                vCoords[iVertex] = coord;
            }

            var ma = ma0 as TimoshenkoBeamMaterial;
            double Ae = ma.Area;
            double Iz = ma.SecondMomentOfArea;
            double Ix = ma.PolarSecondMomentOfArea;
            double rho = ma.MassDensity;
            double E = ma.Young;
            double G = ma.ShearCoefficient;
            double kappa = ma.TimoshenkoShearCoefficient;

            double[] pt1 = vCoords[0];
            double[] pt2 = vCoords[1];
            // y座標は0固定
            System.Diagnostics.Debug.Assert(Math.Abs(pt1[1]) < IvyFEM.Constants.PrecisionLowerLimit);
            System.Diagnostics.Debug.Assert(Math.Abs(pt2[1]) < IvyFEM.Constants.PrecisionLowerLimit);
            double le = Math.Sqrt(
                (pt2[0] - pt1[0]) * (pt2[0] - pt1[0]) +
                (pt2[1] - pt1[1]) * (pt2[1] - pt1[1]));

            // local dof
            int localDof = 2;
            System.Diagnostics.Debug.Assert(localDof == (dDof + rDof));
            int localOffset = dDof;

            IntegrationPoints ipK;
            if (dLineFE.Order == 1 && rLineFE.Order == 1)
            {
                // 低減積分
                ipK = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point1);
                System.Diagnostics.Debug.Assert(ipK.Ls.Length == 1);
            }
            else
            {
                ipK = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);
                System.Diagnostics.Debug.Assert(ipK.Ls.Length == 5);
            }
            for (int ipPt = 0; ipPt < ipK.PointCount; ipPt++)
            {
                double[] L = ipK.Ls[ipPt];
                double[] dN = dLineFE.CalcN(L);
                double[][] dNu = dLineFE.CalcNu(L);
                double[] dNx = dNu[0];
                double[] rN = rLineFE.CalcN(L);
                double[][] rNu = rLineFE.CalcNu(L);
                double[] rNx = rNu[0];
                double lineLen = dLineFE.GetLineLength();
                double weight = ipK.Weights[ipPt];
                double detJWeight = (lineLen / 2.0) * weight;

                // displacement
                for (int row = 0; row < dElemNodeCnt; row++)
                {
                    int rowNodeId = dNodes[row];
                    if (rowNodeId == -1)
                    {
                        continue;
                    }
                    // displacement
                    for (int col = 0; col < dElemNodeCnt; col++)
                    {
                        int colNodeId = dNodes[col];
                        if (colNodeId == -1)
                        {
                            continue;
                        }
                        double kValue = detJWeight * kappa * G * Ae * dNx[row] * dNx[col];
                        K[rowNodeId, colNodeId] += kValue;
                    }
                    // rotation
                    for (int col = 0; col < rElemNodeCnt; col++)
                    {
                        int colNodeId = rNodes[col];
                        if (colNodeId == -1)
                        {
                            continue;
                        }
                        double kValue = -1.0 * detJWeight * kappa * G * Ae * dNx[row] * rN[col];
                        K[rowNodeId, colNodeId + offset] += kValue;
                    }
                }
                // rotation
                for (int row = 0; row < rElemNodeCnt; row++)
                {
                    int rowNodeId = rNodes[row];
                    if (rowNodeId == -1)
                    {
                        continue;
                    }
                    // displacement
                    for (int col = 0; col < dElemNodeCnt; col++)
                    {
                        int colNodeId = dNodes[col];
                        if (colNodeId == -1)
                        {
                            continue;
                        }
                        double kValue = -1.0 * detJWeight * kappa * G * Ae * rN[row] * dNx[col];
                        K[rowNodeId + offset, colNodeId] += kValue;
                    }
                    // rotation
                    for (int col = 0; col < rElemNodeCnt; col++)
                    {
                        int colNodeId = rNodes[col];
                        if (colNodeId == -1)
                        {
                            continue;
                        }
                        double kValue1 = detJWeight * kappa * G * Ae * rN[row] * rN[col];
                        double kValue2 = detJWeight * E * Iz * rNx[row] * rNx[col];
                        K[rowNodeId + offset, colNodeId + offset] +=
                            kValue1 + kValue2;
                    }
                }
            }
            IntegrationPoints ipM = LineFE.GetIntegrationPoints(LineIntegrationPointCount.Point5);
            System.Diagnostics.Debug.Assert(ipM.Ls.Length == 5);
            for (int ipPt = 0; ipPt < ipM.PointCount; ipPt++)
            {
                double[] L = ipM.Ls[ipPt];
                double[] dN = dLineFE.CalcN(L);
                double[][] dNu = dLineFE.CalcNu(L);
                double[] dNx = dNu[0];
                double[] rN = rLineFE.CalcN(L);
                double[][] rNu = rLineFE.CalcNu(L);
                double[] rNx = rNu[0];
                double lineLen = dLineFE.GetLineLength();
                double weight = ipM.Weights[ipPt];
                double detJWeight = (lineLen / 2.0) * weight;

                // displacement
                for (int row = 0; row < dElemNodeCnt; row++)
                {
                    int rowNodeId = dNodes[row];
                    if (rowNodeId == -1)
                    {
                        continue;
                    }
                    // displacement
                    for (int col = 0; col < dElemNodeCnt; col++)
                    {
                        int colNodeId = dNodes[col];
                        if (colNodeId == -1)
                        {
                            continue;
                        }
                        double mValue = detJWeight * rho * Ae * dN[row] * dN[col];
                        M[rowNodeId, colNodeId] += mValue;
                    }
                    // rotation
                    for (int col = 0; col < rElemNodeCnt; col++)
                    {
                        int colNodeId = rNodes[col];
                        if (colNodeId == -1)
                        {
                            continue;
                        }
                        double mValue = 0.0;
                        M[rowNodeId, colNodeId + offset] += mValue;
                    }
                }
                // rotation
                for (int row = 0; row < rElemNodeCnt; row++)
                {
                    int rowNodeId = rNodes[row];
                    if (rowNodeId == -1)
                    {
                        continue;
                    }
                    // displacement
                    for (int col = 0; col < dElemNodeCnt; col++)
                    {
                        int colNodeId = dNodes[col];
                        if (colNodeId == -1)
                        {
                            continue;
                        }
                        double mValue = 0.0;
                        M[rowNodeId + offset, colNodeId] += mValue;
                    }
                    // rotation
                    for (int col = 0; col < rElemNodeCnt; col++)
                    {
                        int colNodeId = rNodes[col];
                        if (colNodeId == -1)
                        {
                            continue;
                        }
                        double mValue = detJWeight * rho * Ix * rN[row] * rN[col];
                        M[rowNodeId + offset, colNodeId + offset] += mValue;
                    }
                }
            }
        }
    }
}
