using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DEigenFEM
    {
        protected IvyFEM.Lapack.DoubleMatrix CalcBeamLocalKe(double le, double E, double I)
        {
            return Elastic2DFEMUtils.CalcBeamLocalKe(le, E, I);
        }

        protected IvyFEM.Lapack.DoubleMatrix CalcBeamLocalMe(double le, double rho, double Ae)
        {
            return Elastic2DFEMUtils.CalcBeamLocalMe(le, rho, Ae);
        }

        protected void CalcBeamElementKMForLine(
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
                if (!(workMa0 is BeamMaterial))
                {
                    return;
                }
            }

            uint dQuantityId = 0; // displacement
            uint rQuantityId = 1; // rotation
            int dNodeCnt = (int)World.GetNodeCount(dQuantityId);
            int rNodeCnt = (int)World.GetNodeCount(rQuantityId);
            int offset = dNodeCnt;
            System.Diagnostics.Debug.Assert(World.GetDof(dQuantityId) == 1);
            System.Diagnostics.Debug.Assert(World.GetDof(rQuantityId) == 1);
            int dDof = 1; // w
            int rDof = 1; //θ

            // Note: feIdは線要素
            LineFE dLineFE = World.GetLineFE(dQuantityId, feId);
            LineFE rLineFE = World.GetLineFE(rQuantityId, feId);
            uint maId = dLineFE.MaterialId;
            if (!World.IsMaterialId(maId))
            {
                return;
            }
            Material ma0 = World.GetMaterial(maId);
            System.Diagnostics.Debug.Assert(ma0 is BeamMaterial);

            int[] dCoIds = dLineFE.NodeCoordIds;
            uint dElemNodeCnt = dLineFE.NodeCount;
            // FIXME: エルミート要素にすべき
            // 暫定で要素マトリクスをこのメソッド内に直接記述
            System.Diagnostics.Debug.Assert(dElemNodeCnt == 2); // 1次要素
            int[] dNodes = new int[dElemNodeCnt];
            for (int iNode = 0; iNode < dElemNodeCnt; iNode++)
            {
                int coId = dCoIds[iNode];
                int nodeId = World.Coord2Node(dQuantityId, coId);
                dNodes[iNode] = nodeId;
            }
            int[] rCoIds = rLineFE.NodeCoordIds;
            uint rElemNodeCnt = rLineFE.NodeCount;
            // FIXME: エルミート要素にすべき
            // 暫定で要素マトリクスをこのメソッド内に直接記述
            System.Diagnostics.Debug.Assert(rElemNodeCnt == 2); // 1次要素
            int[] rNodes = new int[rElemNodeCnt];
            for (int iNode = 0; iNode < dElemNodeCnt; iNode++)
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

            var ma = ma0 as BeamMaterial;
            double Ae = ma.Area;
            double I = ma.SecondMomentOfArea;
            double rho = ma.MassDensity;
            double E = ma.Young;

            double[] pt1 = vCoords[0];
            double[] pt2 = vCoords[1];
            // y座標は0固定
            System.Diagnostics.Debug.Assert(Math.Abs(pt1[1]) < IvyFEM.Constants.PrecisionLowerLimit);
            System.Diagnostics.Debug.Assert(Math.Abs(pt2[1]) < IvyFEM.Constants.PrecisionLowerLimit);
            double le = Math.Sqrt(
                (pt2[0] - pt1[0]) * (pt2[0] - pt1[0]) +
                (pt2[1] - pt1[1]) * (pt2[1] - pt1[1]));

            var ke = CalcBeamLocalKe(le, E, I);
            var me = CalcBeamLocalMe(le, rho, Ae);

            // local dof
            int localDof = 2;
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

                    double kValue = ke[row * localDof, col * localDof];
                    double mValue = me[row * localDof, col * localDof];
                    K[rowNodeId, colNodeId] += kValue;
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

                    double kValue = ke[row * localDof, col * localDof + 1];
                    double mValue = me[row * localDof, col * localDof + 1];
                    K[rowNodeId, colNodeId + offset] += kValue;
                    M[rowNodeId, colNodeId + offset] += mValue;
                }
            }
            // rotation
            for (int row = 0; row < dElemNodeCnt; row++)
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

                    double kValue = ke[row * localDof + 1, col * localDof];
                    double mValue = me[row * localDof + 1, col * localDof];
                    K[rowNodeId + offset, colNodeId] += kValue;
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

                    double kValue = ke[row * localDof + 1, col * localDof + 1];
                    double mValue = me[row * localDof + 1, col * localDof + 1];
                    K[rowNodeId + offset, colNodeId + offset] += kValue;
                    M[rowNodeId + offset, colNodeId + offset] += mValue;
                }
            }
        }
    }
}
