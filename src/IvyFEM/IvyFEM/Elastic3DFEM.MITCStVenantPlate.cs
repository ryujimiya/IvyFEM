using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Animation;
using System.Windows.Media.TextFormatting;

namespace IvyFEM
{
    partial class Elastic3DFEM
    {
        protected void _InitMITCNonlinearPlateNodeValues(int coId, FieldValue uFV)
        {
            Elastic3DFEMUtils.InitMITCNonlinearPlateNodeValues(coId, uFV);
        }

        protected void _InitMITCNonlinearPlateElementNodeValues(
            uint feId, OpenTK.Vector3d[] xPt0s, FieldValue vnFV, FieldValue v1FV, FieldValue v2FV)
        {
            Elastic3DFEMUtils.InitMITCNonlinearPlateElementNodeValues(feId, xPt0s, vnFV, v1FV, v2FV);
        }

        protected void _UpdateMITCNonlinearPlateNodeValues(int coId, double[] curNodeUg, FieldValue uFV)
        {
            Elastic3DFEMUtils.UpdateMITCNonlinearPlateNodeValues(coId, curNodeUg, uFV);
        }

        protected void _UpdateMITCNonlinearPlateElementNodeValues(
            uint feId, int nodeCnt, double[] curUg, FieldValue vnFV, FieldValue v1FV, FieldValue v2FV)
        {
            Elastic3DFEMUtils.UpdateMITCNonlinearPlateElementNodeValues(
                feId, nodeCnt, curUg, vnFV, v1FV, v2FV);
        }

        protected void InitMITCStVenantPlateNodeValues(uint quantityId, int coId)
        {
            if (quantityId != 0)
            {
                return;
            }
            TriangleFE hitTrFE = null;
            IList<uint> feIds = World.GetTriangleFEIdsFromCoord(quantityId, coId);
            foreach (uint feId in feIds)
            {
                TriangleFE workTriFE = World.GetTriangleFE(quantityId, feId);
                uint workMaId = workTriFE.MaterialId;
                if (!World.IsMaterialId(workMaId))
                {
                    continue;
                }
                Material workMa = World.GetMaterial(workMaId);
                if (workMa is MITCStVenantPlateMaterial)
                {
                    hitTrFE = workTriFE;
                    break;
                }
            }
            if (hitTrFE == null)
            {
                return;
            }

            System.Diagnostics.Debug.Assert(TimeIndexForInit >= 0);
            if (TimeIndexForInit != 0)
            {
                // 最初の初期処理以外
                return;
            }

            System.Diagnostics.Debug.Assert(AdditionalValueIds.Count == 4);
            uint uValueId = AdditionalValueIds[0];
            uint vnValueId = AdditionalValueIds[1];
            uint v1ValueId = AdditionalValueIds[2];
            uint v2ValueId = AdditionalValueIds[3];
            FieldValue uFV = World.GetFieldValue(uValueId);
            FieldValue vnFV = World.GetFieldValue(vnValueId);
            FieldValue v1FV = World.GetFieldValue(v1ValueId);
            FieldValue v2FV = World.GetFieldValue(v2ValueId);

            _InitMITCNonlinearPlateNodeValues(coId, uFV);
        }

        protected void InitMITCStVenantPlateElementValues(uint feId)
        {
            {
                uint quantityId0 = 0;

                TriangleFE workTriFE = World.GetTriangleFE(quantityId0, feId);
                uint workMaId = workTriFE.MaterialId;
                if (!World.IsMaterialId(workMaId))
                {
                    return;
                }
                Material workMa = World.GetMaterial(workMaId);
                if (!(workMa is MITCStVenantPlateMaterial))
                {
                    return;
                }
            }

            System.Diagnostics.Debug.Assert(TimeIndexForInit >= 0);
            if (TimeIndexForInit != 0)
            {
                // 最初の初期処理以外
                return;
            }

            System.Diagnostics.Debug.Assert(AdditionalValueIds.Count == 4);
            uint uValueId = AdditionalValueIds[0];
            uint vnValueId = AdditionalValueIds[1];
            uint v1ValueId = AdditionalValueIds[2];
            uint v2ValueId = AdditionalValueIds[3];
            FieldValue uFV = World.GetFieldValue(uValueId);
            FieldValue vnFV = World.GetFieldValue(vnValueId);
            FieldValue v1FV = World.GetFieldValue(v1ValueId);
            FieldValue v2FV = World.GetFieldValue(v2ValueId);

            uint dQuantityId = 0; // displacement
            uint rQuantityId = 1; // rotation
            System.Diagnostics.Debug.Assert(World.GetDof(dQuantityId) == 3);
            System.Diagnostics.Debug.Assert(World.GetDof(rQuantityId) == 3);
            int dDof = 3; // u v w
            int rDof = 3; //θx θy θz
            int dNodeCnt = (int)World.GetNodeCount(dQuantityId);
            int rNodeCnt = (int)World.GetNodeCount(rQuantityId);
            int rOffset = dNodeCnt * dDof;

            TriangleFE dTriFE = World.GetTriangleFE(dQuantityId, feId);
            TriangleFE rTriFE = World.GetTriangleFE(rQuantityId, feId);

            int[] dCoIds = dTriFE.NodeCoordIds;
            uint dElemNodeCnt = dTriFE.NodeCount;
            OpenTK.Vector3d[] dPts = new OpenTK.Vector3d[dElemNodeCnt];
            for (int iNode = 0; iNode < dElemNodeCnt; iNode++)
            {
                int coId = dCoIds[iNode];
                double[] coord = World.GetCoord(dQuantityId, coId);
                dPts[iNode] = new OpenTK.Vector3d(coord[0], coord[1], coord[2]);
            }

            _InitMITCNonlinearPlateElementNodeValues(feId, dPts, vnFV, v1FV, v2FV);
        }

        protected void UpdateMITCStVenantPlateNodeValues(uint quantityId, int coId)
        {
            if (quantityId != 0)
            {
                return;
            }
            TriangleFE hitTrFE = null;
            IList<uint> feIds = World.GetTriangleFEIdsFromCoord(quantityId, coId);
            foreach (uint feId in feIds)
            {
                TriangleFE workTriFE = World.GetTriangleFE(quantityId, feId);
                uint workMaId = workTriFE.MaterialId;
                if (!World.IsMaterialId(workMaId))
                {
                    continue;
                }
                Material workMa = World.GetMaterial(workMaId);
                if (workMa is MITCStVenantPlateMaterial)
                {
                    hitTrFE = workTriFE;
                    break;
                }
            }
            if (hitTrFE == null)
            {
                return;
            }

            System.Diagnostics.Debug.Assert(AdditionalValueIds.Count == 4);
            uint uValueId = AdditionalValueIds[0];
            uint vnValueId = AdditionalValueIds[1];
            uint v1ValueId = AdditionalValueIds[2];
            uint v2ValueId = AdditionalValueIds[3];
            FieldValue uFV = World.GetFieldValue(uValueId);
            FieldValue vnFV = World.GetFieldValue(vnValueId);
            FieldValue v1FV = World.GetFieldValue(v1ValueId);
            FieldValue v2FV = World.GetFieldValue(v2ValueId);

            uint dQuantityId = 0; // displacement
            uint rQuantityId = 1; // rotation
            System.Diagnostics.Debug.Assert(World.GetDof(dQuantityId) == 3);
            System.Diagnostics.Debug.Assert(World.GetDof(rQuantityId) == 3);
            int dDof = 3; // u v w
            int rDof = 3; //θx θyθz
            int localROffset = dDof;
            int localDof = dDof + rDof;
            int dNodeCnt = (int)World.GetNodeCount(dQuantityId);
            int rNodeCnt = (int)World.GetNodeCount(rQuantityId);
            int rOffset = dNodeCnt * dDof;

            int nodeId = World.Coord2Node(quantityId, coId);
            double[] nodeUg = null;
            if (nodeId == -1)
            {
                nodeUg = new double[localDof];
            }
            else
            {
                nodeUg = new double[localDof];
                for (int iDof = 0; iDof < dDof; iDof++)
                {
                    nodeUg[iDof] = U[nodeId * dDof + iDof];
                }
                for (int iDof = 0; iDof < rDof; iDof++)
                {
                    nodeUg[iDof + localROffset] = U[nodeId * rDof + rOffset];
                }
            }
            _UpdateMITCNonlinearPlateNodeValues(coId, nodeUg, uFV);
        }

        protected void UpdateMITCStVenantPlateElementValues(uint feId)
        {
            {
                uint quantityId0 = 0;

                TriangleFE workTriFE = World.GetTriangleFE(quantityId0, feId);
                uint workMaId = workTriFE.MaterialId;
                if (!World.IsMaterialId(workMaId))
                {
                    return;
                }
                Material workMa = World.GetMaterial(workMaId);
                if (!(workMa is MITCStVenantPlateMaterial))
                {
                    return;
                }
            }

            System.Diagnostics.Debug.Assert(AdditionalValueIds.Count == 4);
            uint uValueId = AdditionalValueIds[0];
            uint vnValueId = AdditionalValueIds[1];
            uint v1ValueId = AdditionalValueIds[2];
            uint v2ValueId = AdditionalValueIds[3];
            FieldValue uFV = World.GetFieldValue(uValueId);
            FieldValue vnFV = World.GetFieldValue(vnValueId);
            FieldValue v1FV = World.GetFieldValue(v1ValueId);
            FieldValue v2FV = World.GetFieldValue(v2ValueId);

            uint dQuantityId = 0; // displacement
            uint rQuantityId = 1; // rotation
            System.Diagnostics.Debug.Assert(World.GetDof(dQuantityId) == 3);
            System.Diagnostics.Debug.Assert(World.GetDof(rQuantityId) == 3);
            int dDof = 3; // u v w
            int rDof = 3; //θx θyθz
            int localROffset = dDof;
            int localDof = dDof + rDof;
            int dNodeCnt = (int)World.GetNodeCount(dQuantityId);
            int rNodeCnt = (int)World.GetNodeCount(rQuantityId);
            int rOffset = dNodeCnt * dDof;

            TriangleFE dTriFE = World.GetTriangleFE(dQuantityId, feId);
            TriangleFE rTriFE = World.GetTriangleFE(rQuantityId, feId);

            int[] dCoIds = dTriFE.NodeCoordIds;
            uint dElemNodeCnt = dTriFE.NodeCount;
            int[] rCoIds = rTriFE.NodeCoordIds;
            uint rElemNodeCnt = rTriFE.NodeCount;
            System.Diagnostics.Debug.Assert(dElemNodeCnt == 3);
            System.Diagnostics.Debug.Assert(rElemNodeCnt == 3);
            double[] curUg = new double[dElemNodeCnt * localDof];
            for (int iNode = 0; iNode < dElemNodeCnt; iNode++)
            {
                int coId = dCoIds[iNode];
                int nodeId = World.Coord2Node(dQuantityId, coId);
                if (nodeId == -1)
                {
                    continue;
                }
                for (int iDof = 0; iDof < dDof; iDof++)
                {
                    curUg[iNode * localDof + iDof] = U[nodeId * dDof + iDof];
                }
            }
            for (int iNode = 0; iNode < rElemNodeCnt; iNode++)
            {
                int coId = rCoIds[iNode];
                int nodeId = World.Coord2Node(rQuantityId, coId);
                if (nodeId == -1)
                {
                    continue;
                }
                for (int iDof = 0; iDof < rDof; iDof++)
                {
                    curUg[iNode * localDof + iDof + localROffset] = U[nodeId * rDof + iDof + rOffset];
                }
            }

            _UpdateMITCNonlinearPlateElementNodeValues(
                feId, (int)rElemNodeCnt, curUg, vnFV, v1FV, v2FV);
        }


        protected void CalcMITCStVenantPateKe(
            double lambda, double mu,
            uint feId, TriangleFE dTriFE, IList<int> dCoIds, double h, OpenTK.Vector3d[] xPt0s,
            FieldValue uFV, FieldValue vnFV, FieldValue v1FV, FieldValue v2FV,
            out double[] Qe, out IvyFEM.Lapack.DoubleMatrix Ke)
        {
            Elastic3DFEMUtils.CalcMITCStVenantPateKe(
                lambda, mu,
                feId, dTriFE, dCoIds, h, xPt0s,
                uFV, vnFV, v1FV, v2FV,
                out Qe, out Ke);
        }

        protected void CalcMITCStVenantPlateElementAB(
            uint feId, IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            {
                uint quantityId0 = 0;
                TriangleFE workTriFE = World.GetTriangleFE(quantityId0, feId);
                uint workMaId = workTriFE.MaterialId;
                if (!World.IsMaterialId(workMaId))
                {
                    return;
                }
                Material workMa0 = World.GetMaterial(workMaId);
                if (!(workMa0 is MITCStVenantPlateMaterial))
                {
                    return;
                }
            }

            System.Diagnostics.Debug.Assert(DisplacementQuantityIds.Count == 1);
            System.Diagnostics.Debug.Assert(DisplacementQuantityIds[0] == 0);
            uint dQuantityId = 0; // displacement
            uint rQuantityId = 1; // rotation
            System.Diagnostics.Debug.Assert(World.GetDof(dQuantityId) == 3);
            System.Diagnostics.Debug.Assert(World.GetDof(rQuantityId) == 3);
            int dDof = 3; // u v w
            int rDof = 3; //θx θy θz
            int localROffset = dDof;
            int localDof = dDof + rDof;
            int dNodeCnt = (int)World.GetNodeCount(dQuantityId);
            int rNodeCnt = (int)World.GetNodeCount(rQuantityId);
            int rOffset = dNodeCnt * dDof;

            System.Diagnostics.Debug.Assert(AdditionalValueIds.Count == 4);
            uint uValueId = AdditionalValueIds[0];
            uint vnValueId = AdditionalValueIds[1];
            uint v1ValueId = AdditionalValueIds[2];
            uint v2ValueId = AdditionalValueIds[3];
            FieldValue uFV = World.GetFieldValue(uValueId);
            FieldValue vnFV = World.GetFieldValue(vnValueId);
            FieldValue v1FV = World.GetFieldValue(v1ValueId);
            FieldValue v2FV = World.GetFieldValue(v2ValueId);

            TriangleFE dTriFE = World.GetTriangleFE(dQuantityId, feId);
            TriangleFE rTriFE = World.GetTriangleFE(rQuantityId, feId);
            uint maId = dTriFE.MaterialId;
            if (!World.IsMaterialId(maId))
            {
                return;
            }
            Material ma0 = World.GetMaterial(maId);
            System.Diagnostics.Debug.Assert(ma0 is MITCStVenantPlateMaterial);

            int[] dCoIds = dTriFE.NodeCoordIds;
            uint dElemNodeCnt = dTriFE.NodeCount;
            // 暫定で要素マトリクスをこのメソッド内に直接記述
            System.Diagnostics.Debug.Assert(dElemNodeCnt == 3); // 1次要素
            int[] dNodes = new int[dElemNodeCnt];
            for (int iNode = 0; iNode < dElemNodeCnt; iNode++)
            {
                int coId = dCoIds[iNode];
                int nodeId = World.Coord2Node(dQuantityId, coId);
                dNodes[iNode] = nodeId;
            }
            int[] rCoIds = rTriFE.NodeCoordIds;
            uint rElemNodeCnt = rTriFE.NodeCount;
            // 暫定で要素マトリクスをこのメソッド内に直接記述
            System.Diagnostics.Debug.Assert(rElemNodeCnt == 3); // 1次要素
            int[] rNodes = new int[rElemNodeCnt];
            for (int iNode = 0; iNode < rElemNodeCnt; iNode++)
            {
                int coId = rCoIds[iNode];
                int nodeId = World.Coord2Node(rQuantityId, coId);
                rNodes[iNode] = nodeId;
            }
            OpenTK.Vector3d[] dPts = new OpenTK.Vector3d[dElemNodeCnt];
            for (int iNode = 0; iNode < dElemNodeCnt; iNode++)
            {
                int coId = dCoIds[iNode];
                double[] coord = World.GetCoord(dQuantityId, coId);
                dPts[iNode] = new OpenTK.Vector3d(coord[0], coord[1], coord[2]);
            }

            var ma = ma0 as MITCStVenantPlateMaterial;
            double h = ma.Thickness;
            double rho = ma.MassDensity;
            double lambda = ma.LameLambda;
            double mu = ma.LameMu;

            //------------------------------------------
            // local
            double[] Qe;
            IvyFEM.Lapack.DoubleMatrix Ke;
            CalcMITCStVenantPateKe(
                lambda, mu,
                feId, dTriFE, dCoIds, h, dPts,
                uFV, vnFV, v1FV, v2FV,
                out Qe, out Ke);
            // Note: グローバル座標へ変換済み

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
                    for (int rowDof = 0; rowDof < dDof; rowDof++)
                    {
                        for (int colDof = 0; colDof < dDof; colDof++)
                        {
                            double kValue = Ke[row * localDof + rowDof, col * localDof + colDof];
                            A[rowNodeId * dDof + rowDof, colNodeId * dDof + colDof] += kValue;
                        }
                    }
                }
                // rotation
                for (int col = 0; col < rElemNodeCnt; col++)
                {
                    int colNodeId = rNodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    for (int rowDof = 0; rowDof < dDof; rowDof++)
                    {
                        for (int colDof = 0; colDof < rDof; colDof++)
                        {
                            double kValue = Ke[row * localDof + rowDof, col * localDof + colDof + localROffset];
                            A[rowNodeId * dDof + rowDof, colNodeId * rDof + colDof + rOffset] += kValue;
                        }
                    }
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
                    for (int rowDof = 0; rowDof < rDof; rowDof++)
                    {
                        for (int colDof = 0; colDof < dDof; colDof++)
                        {
                            double kValue = Ke[row * localDof + rowDof + localROffset, col * localDof + colDof];
                            A[rowNodeId * rDof + rowDof + rOffset, colNodeId * dDof + colDof] += kValue;
                        }
                    }
                }
                // rotation
                for (int col = 0; col < rElemNodeCnt; col++)
                {
                    int colNodeId = rNodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    for (int rowDof = 0; rowDof < rDof; rowDof++)
                    {
                        for (int colDof = 0; colDof < rDof; colDof++)
                        {
                            double kValue =
                                Ke[row * localDof + rowDof + localROffset, col * localDof + colDof + localROffset];
                            A[rowNodeId * rDof + rowDof + rOffset, colNodeId * rDof + colDof + rOffset] += kValue;
                        }
                    }
                }
            }

            // internal force
            for (int row = 0; row < dElemNodeCnt; row++)
            {
                int rowNodeId = dNodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                for (int rowDof = 0; rowDof < dDof; rowDof++)
                {
                    B[rowNodeId * dDof + rowDof] += -Qe[row * localDof + rowDof];
                }
            }
            for (int row = 0; row < rElemNodeCnt; row++)
            {
                int rowNodeId = rNodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                for (int rowDof = 0; rowDof < rDof; rowDof++)
                {
                    B[rowNodeId * rDof + rowDof + rOffset] += -Qe[row * localDof + rowDof + localROffset];
                }
            }
        }
    }
}
