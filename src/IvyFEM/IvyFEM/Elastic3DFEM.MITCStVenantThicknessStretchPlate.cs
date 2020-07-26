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
        protected void _InitMITCThicknessStretchPlateNodeValues(int coId, FieldValue uFV)
        {
            Elastic3DFEMUtils.InitMITCThicknessStretchPlateNodeValues(coId, uFV);
        }

        protected void _InitMITCThicknessStretchPlateElementNodeValues(
            uint feId, OpenTK.Vector3d[] xPt0s,
            FieldValue vnFV, FieldValue v1FV, FieldValue v2FV, FieldValue lambdaFV)
        {
            Elastic3DFEMUtils.InitMITCThicknessStretchPlateElementNodeValues(
                feId, xPt0s, vnFV, v1FV, v2FV, lambdaFV);
        }

        protected void _UpdateMITCThicknessStretchPlateNodeValues(int coId, double[] curNodeUg, FieldValue uFV)
        {
            Elastic3DFEMUtils.UpdateMITCThicknessStretchPlateNodeValues(coId, curNodeUg, uFV);
        }

        protected void _UpdateMITCThicknessStretchPlateElementNodeValues(
            uint feId, int nodeCnt, double[] curUg,
            FieldValue vnFV, FieldValue v1FV, FieldValue v2FV, FieldValue lambdaFV)
        {
            Elastic3DFEMUtils.UpdateMITCThicknessStretchPlateElementNodeValues(
                feId, nodeCnt, curUg, vnFV, v1FV, v2FV, lambdaFV);
        }

        protected void InitMITCStVenantThicknessStretchPlateNodeValues(uint quantityId, int coId)
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
                if (workMa is MITCStVenantThicknessStretchPlateMaterial)
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

            System.Diagnostics.Debug.Assert(AdditionalValueIds.Count == 5);
            uint uValueId = AdditionalValueIds[0];
            uint vnValueId = AdditionalValueIds[1];
            uint v1ValueId = AdditionalValueIds[2];
            uint v2ValueId = AdditionalValueIds[3];
            uint lValueId = AdditionalValueIds[4];
            FieldValue uFV = World.GetFieldValue(uValueId);
            FieldValue vnFV = World.GetFieldValue(vnValueId);
            FieldValue v1FV = World.GetFieldValue(v1ValueId);
            FieldValue v2FV = World.GetFieldValue(v2ValueId);
            FieldValue lambdaFV = World.GetFieldValue(lValueId);

            _InitMITCThicknessStretchPlateNodeValues(coId, uFV);
        }

        protected void InitMITCStVenantThicknessStretchPlateElementValues(uint feId)
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
                if (!(workMa is MITCStVenantThicknessStretchPlateMaterial))
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

            System.Diagnostics.Debug.Assert(AdditionalValueIds.Count == 5);
            uint uValueId = AdditionalValueIds[0];
            uint vnValueId = AdditionalValueIds[1];
            uint v1ValueId = AdditionalValueIds[2];
            uint v2ValueId = AdditionalValueIds[3];
            uint lValueId = AdditionalValueIds[4];
            FieldValue uFV = World.GetFieldValue(uValueId);
            FieldValue vnFV = World.GetFieldValue(vnValueId);
            FieldValue v1FV = World.GetFieldValue(v1ValueId);
            FieldValue v2FV = World.GetFieldValue(v2ValueId);
            FieldValue lambdaFV = World.GetFieldValue(lValueId);

            uint dQuantityId = 0; // displacement
            uint rQuantityId = 1; // rotation
            uint lQuantityId = 2; // stretch
            System.Diagnostics.Debug.Assert(World.GetDof(dQuantityId) == 3);
            System.Diagnostics.Debug.Assert(World.GetDof(rQuantityId) == 3);
            System.Diagnostics.Debug.Assert(World.GetDof(lQuantityId) == 2);
            uint elemNodeCnt = 3;
            int dDof = 3; // u v w
            int rDof = 3; //θx θy θz
            int lDof = 2;  // deltaλ0、deltaλ1
            int additionalDof = lDof;
            int localROffset = dDof;
            int localDof = dDof + rDof;
            int localLOffset = (int)(elemNodeCnt * localDof);
            int dNodeCnt = (int)World.GetNodeCount(dQuantityId);
            int rNodeCnt = (int)World.GetNodeCount(rQuantityId);
            int lNodeCnt = (int)World.GetNodeCount(lQuantityId);
            int rOffset = dNodeCnt * dDof;
            int lOffset = rOffset + rNodeCnt * rDof;

            TriangleFE dTriFE = World.GetTriangleFE(dQuantityId, feId);
            TriangleFE rTriFE = World.GetTriangleFE(rQuantityId, feId);
            TriangleFE lTriFE = World.GetTriangleFE(lQuantityId, feId);

            int[] dCoIds = dTriFE.NodeCoordIds;
            uint dElemNodeCnt = dTriFE.NodeCount;
            OpenTK.Vector3d[] dPts = new OpenTK.Vector3d[dElemNodeCnt];
            for (int iNode = 0; iNode < dElemNodeCnt; iNode++)
            {
                int coId = dCoIds[iNode];
                double[] coord = World.GetCoord(dQuantityId, coId);
                dPts[iNode] = new OpenTK.Vector3d(coord[0], coord[1], coord[2]);
            }

            _InitMITCThicknessStretchPlateElementNodeValues(feId, dPts, vnFV, v1FV, v2FV, lambdaFV);
        }

        protected void UpdateMITCStVenantThicknessStretchPlateNodeValues(uint quantityId, int coId)
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
                if (workMa is MITCStVenantThicknessStretchPlateMaterial)
                {
                    hitTrFE = workTriFE;
                    break;
                }
            }
            if (hitTrFE == null)
            {
                return;
            }

            System.Diagnostics.Debug.Assert(AdditionalValueIds.Count == 5);
            uint uValueId = AdditionalValueIds[0];
            uint vnValueId = AdditionalValueIds[1];
            uint v1ValueId = AdditionalValueIds[2];
            uint v2ValueId = AdditionalValueIds[3];
            uint lValueId = AdditionalValueIds[4];
            FieldValue uFV = World.GetFieldValue(uValueId);
            FieldValue vnFV = World.GetFieldValue(vnValueId);
            FieldValue v1FV = World.GetFieldValue(v1ValueId);
            FieldValue v2FV = World.GetFieldValue(v2ValueId);
            FieldValue lambdaFV = World.GetFieldValue(lValueId);

            uint dQuantityId = 0; // displacement
            uint rQuantityId = 1; // rotation
            uint lQuantityId = 2; // stretch
            System.Diagnostics.Debug.Assert(World.GetDof(dQuantityId) == 3);
            System.Diagnostics.Debug.Assert(World.GetDof(rQuantityId) == 3);
            System.Diagnostics.Debug.Assert(World.GetDof(lQuantityId) == 2);
            uint elemNodeCnt = 3;
            int dDof = 3; // u v w
            int rDof = 3; //θx θy θz
            int lDof = 2;  // deltaλ0、deltaλ1
            int additionalDof = lDof;
            int localROffset = dDof;
            int localDof = dDof + rDof;
            int localLOffset = (int)(elemNodeCnt * localDof);
            int dNodeCnt = (int)World.GetNodeCount(dQuantityId);
            int rNodeCnt = (int)World.GetNodeCount(rQuantityId);
            int lNodeCnt = (int)World.GetNodeCount(lQuantityId);
            int rOffset = dNodeCnt * dDof;
            int lOffset = rOffset + rNodeCnt * rDof;

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
                // lambdaはない
            }
            _UpdateMITCThicknessStretchPlateNodeValues(coId, nodeUg, uFV);
        }

        protected void UpdateMITCStVenantThicknessStretchPlateElementValues(uint feId)
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
                if (!(workMa is MITCStVenantThicknessStretchPlateMaterial))
                {
                    return;
                }
            }

            System.Diagnostics.Debug.Assert(AdditionalValueIds.Count == 5);
            uint uValueId = AdditionalValueIds[0];
            uint vnValueId = AdditionalValueIds[1];
            uint v1ValueId = AdditionalValueIds[2];
            uint v2ValueId = AdditionalValueIds[3];
            uint lValueId = AdditionalValueIds[4];
            FieldValue uFV = World.GetFieldValue(uValueId);
            FieldValue vnFV = World.GetFieldValue(vnValueId);
            FieldValue v1FV = World.GetFieldValue(v1ValueId);
            FieldValue v2FV = World.GetFieldValue(v2ValueId);
            FieldValue lambdaFV = World.GetFieldValue(lValueId);

            uint dQuantityId = 0; // displacement
            uint rQuantityId = 1; // rotation
            uint lQuantityId = 2; // stretch
            System.Diagnostics.Debug.Assert(World.GetDof(dQuantityId) == 3);
            System.Diagnostics.Debug.Assert(World.GetDof(rQuantityId) == 3);
            System.Diagnostics.Debug.Assert(World.GetDof(lQuantityId) == 2);
            uint elemNodeCnt = 3;
            int dDof = 3; // u v w
            int rDof = 3; //θx θy θz
            int lDof = 2;  // deltaλ0、deltaλ1
            int additionalDof = lDof;
            int localROffset = dDof;
            int localDof = dDof + rDof;
            int localLOffset = (int)(elemNodeCnt * localDof);
            int dNodeCnt = (int)World.GetNodeCount(dQuantityId);
            int rNodeCnt = (int)World.GetNodeCount(rQuantityId);
            int lNodeCnt = (int)World.GetNodeCount(lQuantityId);
            int rOffset = dNodeCnt * dDof;
            int lOffset = rOffset + rNodeCnt * rDof;

            TriangleFE dTriFE = World.GetTriangleFE(dQuantityId, feId);
            TriangleFE rTriFE = World.GetTriangleFE(rQuantityId, feId);
            TriangleFE lTriFE = World.GetTriangleFE(lQuantityId, feId);

            int[] dCoIds = dTriFE.NodeCoordIds;
            uint dElemNodeCnt = dTriFE.NodeCount;
            int[] rCoIds = rTriFE.NodeCoordIds;
            uint rElemNodeCnt = rTriFE.NodeCount;
            int[] lCoIds = lTriFE.NodeCoordIds;
            uint lElemNodeCnt = lTriFE.NodeCount;
            System.Diagnostics.Debug.Assert(dElemNodeCnt == 3);
            System.Diagnostics.Debug.Assert(rElemNodeCnt == 3);
            System.Diagnostics.Debug.Assert(lElemNodeCnt == 1);
            double[] curUg = new double[dElemNodeCnt * localDof + additionalDof];
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
            for (int iNode = 0; iNode < lElemNodeCnt; iNode++)
            {
                int coId = lCoIds[iNode];
                int nodeId = World.Coord2Node(lQuantityId, coId);
                if (nodeId == -1)
                {
                    continue;
                }
                for (int iDof = 0; iDof < lDof; iDof++)
                {
                    curUg[iDof + localLOffset] = U[nodeId * lDof + iDof + lOffset];
                }
            }

            _UpdateMITCThicknessStretchPlateElementNodeValues(
                feId, (int)rElemNodeCnt, curUg, vnFV, v1FV, v2FV, lambdaFV);
        }


        protected void CalcMITCStVenantThicknessStretchPateKe(
            double c1, double c2,
            uint feId, TriangleFE dTriFE, IList<int> dCoIds, double h, OpenTK.Vector3d[] xPt0s,
            double[] curUg,
            FieldValue uFV, FieldValue vnFV, FieldValue v1FV, FieldValue v2FV, FieldValue lambdaFV,
            out double[] Qe0, out double[] Qe2, out IvyFEM.Lapack.DoubleMatrix Ke1, out IvyFEM.Lapack.DoubleMatrix Ke2)
        {
            Elastic3DFEMUtils.CalcMITCStVenantThicknessStretchPateKe(
                c1, c2,
                feId, dTriFE, dCoIds, h, xPt0s,
                uFV, vnFV, v1FV, v2FV, lambdaFV,
                curUg,
                out Qe0, out Qe2, out Ke1, out Ke2);
        }

        protected void CalcMITCStVenantThicknessStretchPlateElementAB(
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
                if (!(workMa0 is MITCStVenantThicknessStretchPlateMaterial))
                {
                    return;
                }
            }

            uint dQuantityId = 0; // displacement
            uint rQuantityId = 1; // rotation
            uint lQuantityId = 2; // stretch
            System.Diagnostics.Debug.Assert(World.GetDof(dQuantityId) == 3);
            System.Diagnostics.Debug.Assert(World.GetDof(rQuantityId) == 3);
            System.Diagnostics.Debug.Assert(World.GetDof(lQuantityId) == 2);
            uint elemNodeCnt = 3;
            int dDof = 3; // u v w
            int rDof = 3; //θx θy θz
            int lDof = 2;  // deltaλ0、deltaλ1
            int additionalDof = lDof;
            int localROffset = dDof;
            int localDof = dDof + rDof;
            int localLOffset = (int)(elemNodeCnt * localDof);
            int dNodeCnt = (int)World.GetNodeCount(dQuantityId);
            int rNodeCnt = (int)World.GetNodeCount(rQuantityId);
            int lNodeCnt = (int)World.GetNodeCount(lQuantityId);
            int rOffset = dNodeCnt * dDof;
            int lOffset = rOffset + rNodeCnt * rDof;

            System.Diagnostics.Debug.Assert(AdditionalValueIds.Count == 5);
            uint uValueId = AdditionalValueIds[0];
            uint vnValueId = AdditionalValueIds[1];
            uint v1ValueId = AdditionalValueIds[2];
            uint v2ValueId = AdditionalValueIds[3];
            uint lValueId = AdditionalValueIds[4];
            FieldValue uFV = World.GetFieldValue(uValueId);
            FieldValue vnFV = World.GetFieldValue(vnValueId);
            FieldValue v1FV = World.GetFieldValue(v1ValueId);
            FieldValue v2FV = World.GetFieldValue(v2ValueId);
            FieldValue lambdaFV = World.GetFieldValue(lValueId);

            TriangleFE dTriFE = World.GetTriangleFE(dQuantityId, feId);
            TriangleFE rTriFE = World.GetTriangleFE(rQuantityId, feId);
            TriangleFE lTriFE = World.GetTriangleFE(lQuantityId, feId);
            uint maId = dTriFE.MaterialId;
            if (!World.IsMaterialId(maId))
            {
                return;
            }
            Material ma0 = World.GetMaterial(maId);
            System.Diagnostics.Debug.Assert(ma0 is MITCStVenantThicknessStretchPlateMaterial);

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
            int[] lCoIds = lTriFE.NodeCoordIds;
            uint lElemNodeCnt = lTriFE.NodeCount;
            // 暫定で要素マトリクスをこのメソッド内に直接記述
            System.Diagnostics.Debug.Assert(lElemNodeCnt == 1); // constant次要素
            int[] lNodes = new int[lElemNodeCnt];
            for (int iNode = 0; iNode < lElemNodeCnt; iNode++)
            {
                int coId = lCoIds[iNode];
                int nodeId = World.Coord2Node(lQuantityId, coId);
                lNodes[iNode] = nodeId;
            }
            OpenTK.Vector3d[] dPts = new OpenTK.Vector3d[dElemNodeCnt];
            for (int iNode = 0; iNode < dElemNodeCnt; iNode++)
            {
                int coId = dCoIds[iNode];
                double[] coord = World.GetCoord(dQuantityId, coId);
                dPts[iNode] = new OpenTK.Vector3d(coord[0], coord[1], coord[2]);
            }
            double[] curUg = new double[dElemNodeCnt * localDof + additionalDof];
            for (int iNode = 0; iNode < dElemNodeCnt; iNode++)
            {
                int nodeId = dNodes[iNode];
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
                int nodeId = rNodes[iNode];
                if (nodeId == -1)
                {
                    continue;
                }
                for (int iDof = 0; iDof < rDof; iDof++)
                {
                    curUg[iNode * localDof + iDof + localROffset] = U[nodeId * rDof + iDof + rOffset];
                }
            }
            for (int iNode = 0; iNode < lElemNodeCnt; iNode++)
            {
                int nodeId = lNodes[iNode];
                if (nodeId == -1)
                {
                    continue;
                }
                for (int iDof = 0; iDof < lDof; iDof++)
                {
                    curUg[iDof + localLOffset] = U[nodeId * lDof + iDof + lOffset];
                }
            }

            var ma = ma0 as MITCStVenantThicknessStretchPlateMaterial;
            double h = ma.Thickness;
            double rho = ma.MassDensity;
            double lambda = ma.LameLambda;
            double mu = ma.LameMu;

            //------------------------------------------
            // local
            double[] Qe0;
            double[] Qe2;
            IvyFEM.Lapack.DoubleMatrix Ke1;
            IvyFEM.Lapack.DoubleMatrix Ke2;
            CalcMITCStVenantThicknessStretchPateKe(
                lambda, mu,
                feId, dTriFE, dCoIds, h, dPts,
                curUg,
                uFV, vnFV, v1FV, v2FV, lambdaFV,
                out Qe0, out Qe2, out Ke1, out Ke2);
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
                            double k1Value = Ke1[row * localDof + rowDof, col * localDof + colDof];
                            double k2Value = Ke2[row * localDof + rowDof, col * localDof + colDof];
                            A[rowNodeId * dDof + rowDof, colNodeId * dDof + colDof] +=
                                k1Value + k2Value;
                            B[rowNodeId * dDof + rowDof] += k2Value * U[colNodeId * dDof + colDof];
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
                            double k1Value = Ke1[row * localDof + rowDof, col * localDof + colDof + localROffset];
                            double k2Value = Ke2[row * localDof + rowDof, col * localDof + colDof + localROffset];
                            A[rowNodeId * dDof + rowDof, colNodeId * rDof + colDof + rOffset] +=
                                k1Value + k2Value;
                            B[rowNodeId * dDof + rowDof] += k2Value * U[colNodeId * rDof + colDof + rOffset];
                        }
                    }
                }
                // stretch
                for (int col = 0; col < lElemNodeCnt; col++)
                {
                    int colNodeId = lNodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    for (int rowDof = 0; rowDof < dDof; rowDof++)
                    {
                        for (int colDof = 0; colDof < lDof; colDof++)
                        {
                            double k1Value = Ke1[row * localDof + rowDof, colDof + localLOffset];
                            double k2Value = Ke2[row * localDof + rowDof, colDof + localLOffset];
                            A[rowNodeId * dDof + rowDof, colNodeId * lDof + colDof + lOffset] +=
                                k1Value + k2Value;
                            B[rowNodeId * dDof + rowDof] += k2Value * U[colNodeId * lDof + colDof + lOffset];
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
                            double k1Value = Ke1[row * localDof + rowDof + localROffset, col * localDof + colDof];
                            double k2Value = Ke2[row * localDof + rowDof + localROffset, col * localDof + colDof];
                            A[rowNodeId * rDof + rowDof + rOffset, colNodeId * dDof + colDof] +=
                                k1Value + k2Value;
                            B[rowNodeId * rDof + rowDof + rOffset] += k2Value * U[colNodeId * dDof + colDof];
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
                            double k1Value =
                                Ke1[row * localDof + rowDof + localROffset, col * localDof + colDof + localROffset];
                            double k2Value =
                                Ke2[row * localDof + rowDof + localROffset, col * localDof + colDof + localROffset];
                            A[rowNodeId * rDof + rowDof + rOffset, colNodeId * rDof + colDof + rOffset] +=
                                k1Value + k2Value;
                            B[rowNodeId * rDof + rowDof + rOffset] += k2Value * U[colNodeId * rDof + colDof + rOffset];
                        }
                    }
                }
                // stretch
                for (int col = 0; col < lElemNodeCnt; col++)
                {
                    int colNodeId = lNodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    for (int rowDof = 0; rowDof < rDof; rowDof++)
                    {
                        for (int colDof = 0; colDof < lDof; colDof++)
                        {
                            double k1Value =
                                Ke1[row * localDof + rowDof + localROffset, colDof + localLOffset];
                            double k2Value =
                                Ke2[row * localDof + rowDof + localROffset, colDof + localLOffset];
                            A[rowNodeId * rDof + rowDof + rOffset, colNodeId * lDof + colDof + lOffset] +=
                                k1Value + k2Value;
                            B[rowNodeId * rDof + rowDof + rOffset] += k2Value * U[colNodeId * lDof + colDof + lOffset];
                        }
                    }
                }
            }
            // stretch
            for (int row = 0; row < lElemNodeCnt; row++)
            {
                int rowNodeId = lNodes[row];
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
                    for (int rowDof = 0; rowDof < lDof; rowDof++)
                    {
                        for (int colDof = 0; colDof < dDof; colDof++)
                        {
                            double k1Value = Ke1[rowDof + localLOffset, col * localDof + colDof];
                            double k2Value = Ke2[rowDof + localLOffset, col * localDof + colDof];
                            A[rowNodeId * lDof + rowDof + lOffset, colNodeId * dDof + colDof] +=
                                k1Value + k2Value;
                            B[rowNodeId * lDof + rowDof + lOffset] += k2Value * U[colNodeId * dDof + colDof];
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
                    for (int rowDof = 0; rowDof < lDof; rowDof++)
                    {
                        for (int colDof = 0; colDof < rDof; colDof++)
                        {
                            double k1Value =
                                Ke1[rowDof + localLOffset, col * localDof + colDof + localROffset];
                            double k2Value =
                                Ke2[rowDof + localLOffset, col * localDof + colDof + localROffset];
                            A[rowNodeId * lDof + rowDof + lOffset, colNodeId * rDof + colDof + rOffset] +=
                                k1Value + k2Value;
                            B[rowNodeId * lDof + rowDof + lOffset] += k2Value * U[colNodeId * rDof + colDof + rOffset];
                        }
                    }
                }
                // stretch
                for (int col = 0; col < lElemNodeCnt; col++)
                {
                    int colNodeId = lNodes[col];
                    if (colNodeId == -1)
                    {
                        continue;
                    }
                    for (int rowDof = 0; rowDof < lDof; rowDof++)
                    {
                        for (int colDof = 0; colDof < lDof; colDof++)
                        {
                            double k1Value =
                                Ke1[rowDof + localLOffset, colDof + localLOffset];
                            double k2Value =
                                Ke2[rowDof + localLOffset, colDof + localLOffset];
                            A[rowNodeId * lDof + rowDof + lOffset, colNodeId * lDof + colDof + lOffset] += 
                                k1Value + k2Value;
                            B[rowNodeId * lDof + rowDof + lOffset] += k2Value * U[colNodeId * lDof + colDof + lOffset];
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
                    B[rowNodeId * dDof + rowDof] += -Qe0[row * localDof + rowDof] - Qe2[row * localDof + rowDof];
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
                    B[rowNodeId * rDof + rowDof + rOffset] +=
                        -Qe0[row * localDof + rowDof + localROffset] - Qe2[row * localDof + rowDof + localROffset];
                }
            }
            for (int row = 0; row < lElemNodeCnt; row++)
            {
                int rowNodeId = lNodes[row];
                if (rowNodeId == -1)
                {
                    continue;
                }
                for (int rowDof = 0; rowDof < lDof; rowDof++)
                {
                    B[rowNodeId * lDof + rowDof + lOffset] += -Qe0[rowDof + localLOffset] - Qe2[rowDof + localLOffset];
                }
            }
        }
    }
}
