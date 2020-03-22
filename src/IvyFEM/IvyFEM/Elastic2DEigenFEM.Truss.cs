using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DEigenFEM
    {
        protected IvyFEM.Lapack.DoubleMatrix CalcTrussLocalKe(double le, double E, double Ae)
        {
            return Elastic2DFEMUtils.CalcTrussLocalKe(le, E, Ae);
        }

        protected IvyFEM.Lapack.DoubleMatrix CalcTrussLocalMe(double le, double rho, double Ae)
        {
            return Elastic2DFEMUtils.CalcTrussLocalMe(le, rho, Ae);
        }

        protected void CalcTrussElementKMForLine(
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
                if (!(workMa0 is TrussMaterial))
                {
                    return;
                }
            }

            uint quantityId = 0;
            System.Diagnostics.Debug.Assert(World.GetDof(quantityId) == 2);
            int dof = 2;
            int nodeCnt = (int)World.GetNodeCount(quantityId);

            // Note: feIdは線要素
            LineFE lineFE = World.GetLineFE(quantityId, feId);
            uint maId = lineFE.MaterialId;
            if (!World.IsMaterialId(maId))
            {
                return;
            }
            Material ma0 = World.GetMaterial(maId);
            if (!(ma0 is TrussMaterial))
            {
                return;
            }
            int[] coIds = lineFE.NodeCoordIds;
            uint elemNodeCnt = lineFE.NodeCount;
            // FIXME:
            // 暫定で要素マトリクスをこのメソッド内に直接記述
            System.Diagnostics.Debug.Assert(elemNodeCnt == 2); // 1次要素
            int[] nodes = new int[elemNodeCnt];
            for (int iNode = 0; iNode < elemNodeCnt; iNode++)
            {
                int coId = coIds[iNode];
                int nodeId = World.Coord2Node(quantityId, coId);
                nodes[iNode] = nodeId;
            }
            int[] vCoIds = lineFE.VertexCoordIds;
            uint elemVertexCnt = lineFE.VertexCount;
            double[][] vCoords = new double[elemVertexCnt][];
            for (int iVertex = 0; iVertex < elemNodeCnt; iVertex++)
            {
                int coId = vCoIds[iVertex];
                double[] coord = World.GetCoord(quantityId, coId);
                vCoords[iVertex] = coord;
            }

            var ma = ma0 as TrussMaterial;
            double Ae = ma.Area;
            double rho = ma.MassDensity;
            double E = ma.Young;

            double[] pt1 = vCoords[0];
            double[] pt2 = vCoords[1];
            double le = Math.Sqrt(
                (pt2[0] - pt1[0]) * (pt2[0] - pt1[0]) +
                (pt2[1] - pt1[1]) * (pt2[1] - pt1[1]));
            double cosXX = (pt2[0] - pt1[0]) / le;
            double cosXY = (pt2[1] - pt1[1]) / le;

            var T = new IvyFEM.Lapack.DoubleMatrix(2, 4);
            T[0, 0] = cosXX;
            T[0, 1] = cosXY;
            T[1, 2] = cosXX;
            T[1, 3] = cosXY;
            var transT = IvyFEM.Lapack.DoubleMatrix.Transpose(T);

            var localKe = CalcTrussLocalKe(le, E, Ae);
            var localMe = CalcTrussLocalMe(le, rho, Ae);
            var Ke = (transT * localKe) * T;
            var Me = (transT * localMe) * T;

            for (int row = 0; row < elemNodeCnt; row++)
            {
                int rowNodeId = nodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                for (int col = 0; col < elemNodeCnt; col++)
                {
                    int colNodeId = nodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }

                    for (int rowDof = 0; rowDof < dof; rowDof++)
                    {
                        for (int colDof = 0; colDof < dof; colDof++)
                        {
                            double kValue = Ke[row * dof + rowDof, col * dof + colDof];
                            double mValue = Me[row * dof + rowDof, col * dof + colDof];
                            K[rowNodeId * dof + rowDof, colNodeId * dof + colDof] += kValue;
                            M[rowNodeId * dof + rowDof, colNodeId * dof + colDof] += mValue;
                        }
                    }
                }
            }
        }
    }
}
