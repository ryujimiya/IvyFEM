using System;
using System.Collections.Generic;
using System.Drawing.Imaging;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class FEWorld
    {
        public IMesher Mesh { get; set; } = null;
        public uint Dimension => Mesh != null ? Mesh.Dimension : 0;
        private IList<double> VertexCoords = new List<double>();
        private ObjectArray<Material> MaterialArray = new ObjectArray<Material>();
        private Dictionary<uint, uint> CadEdge2Material = new Dictionary<uint, uint>();
        private Dictionary<uint, uint> CadLoop2Material = new Dictionary<uint, uint>();
        private Dictionary<uint, uint> CadSolid2Material = new Dictionary<uint, uint>();

        private ObjectArray<FieldValue> FieldValueArray = new ObjectArray<FieldValue>();

        public TriangleIntegrationPointCount TriIntegrationPointCount { get; set; } =
            TriangleIntegrationPointCount.Point7;
        public TetrahedronIntegrationPointCount TetIntegrationPointCount { get; set; } =
            TetrahedronIntegrationPointCount.Point5;
        // 2D
        public double RotAngle { get; set; } = 0.0;
        public double[] RotOrigin { get; set; } = null;

        private IList<FEWorldQuantity> Quantitys = new List<FEWorldQuantity>();

        public FEWorld()
        {

        }

        public uint AddQuantity(uint dof, uint feOrder, FiniteElementType feType, int idBaseOffset = 0)
        {
            uint id = (uint)Quantitys.Count;
            var quantity = new FEWorldQuantity(id, Dimension, dof, feOrder, feType, idBaseOffset);
            Quantitys.Add(quantity);
            return id;
        }

        public void ClearQuantity()
        {
            Quantitys.Clear();
        }

        public int GetQuantityCount()
        {
            return Quantitys.Count;
        }

        public void IsPortOnly(uint quantityId, bool value)
        {
            Quantitys[(int)quantityId].IsPortOnly = value;
        }

        public bool IsPointOnly(uint quantityId)
        {
            return Quantitys[(int)quantityId].IsPortOnly;
        }

        public void Clear()
        {
            Mesh = null;
            MaterialArray.Clear();
            CadEdge2Material.Clear();
            CadLoop2Material.Clear();
            CadSolid2Material.Clear();
            FieldValueArray.Clear();

            ClearElements();

            Quantitys.Clear();
        }

        private void ClearElements()
        {
            VertexCoords.Clear();
            foreach (var quantity in Quantitys)
            {
                quantity.ClearElements();
            }
        }

        public uint GetVertexCoordCount()
        {
            return (uint)VertexCoords.Count / Dimension;
        }

        public double[] GetVertexCoord(int coId)
        {
            double[] coord = _GetVertexCoord(coId);

            // 回転移動
            if (Dimension == 2)
            {
                coord = CadUtils2D.GetRotCoord2D(coord, RotAngle, RotOrigin);
            }
            else
            {
                // TODO:
            }
            return coord;
        }

        private double[] _GetVertexCoord(int coId)
        {
            System.Diagnostics.Debug.Assert(coId * Dimension + (Dimension - 1) < VertexCoords.Count);
            double[] coord = new double[Dimension];
            for (int iDim = 0; iDim < Dimension; iDim++)
            {
                coord[iDim] = VertexCoords[(int)(coId * Dimension + iDim)];
            }
            return coord;
        }

        public uint GetFEOrder(uint quantityId)
        {
            return Quantitys[(int)quantityId].FEOrder;
        }

        public IList<FieldFixedCad> GetZeroFieldFixedCads(uint quantityId)
        {
            return Quantitys[(int)quantityId].ZeroFieldFixedCads;
        }

        public void SetZeroFieldFixedCads(uint quantityId, IList<FieldFixedCad> fixedCads)
        {
            Quantitys[(int)quantityId].ZeroFieldFixedCads = fixedCads;
        }

        public IList<FieldFixedCad> GetFieldFixedCads(uint quantityId)
        {
            return Quantitys[(int)quantityId].FieldFixedCads;
        }

        public void SetFieldFixedCads(uint quantityId, IList<FieldFixedCad> fixedCads)
        {
            Quantitys[(int)quantityId].FieldFixedCads = fixedCads;
        }

        public IList<FieldFixedCad> GetForceFieldFixedCads(uint quantityId)
        {
            return Quantitys[(int)quantityId].ForceFieldFixedCads;
        }

        public void SetForceFieldFixedCads(uint quantityId, IList<FieldFixedCad> fixedCads)
        {
            Quantitys[(int)quantityId].ForceFieldFixedCads = fixedCads;
        }

        public int GetMultipointConstraintCount(uint quantityId)
        {
            return Quantitys[(int)quantityId].GetMultipointConstraintCount();
        }

        public int AddMultipointConstraint(uint quantityId, MultipointConstraint mpc)
        {
            int index = Quantitys[(int)quantityId].AddMultipointConstraint(mpc);
            return index;
        }

        public MultipointConstraint GetMultipointConstraint(uint quantityId, int index)
        {
            return Quantitys[(int)quantityId].GetMultipointConstraint(index);
        }

        public void ClearMultipointConstraint(uint quantityId)
        {
            Quantitys[(int)quantityId].ClearMultipointConstraint();
        }

        public IList<uint> GetContactSlaveCadIds(uint quantityId)
        {
            return Quantitys[(int)quantityId].ContactSlaveCadIds;
        }

        public void SetContactSlaveCadIds(uint quantityId, IList<uint> slaveIds)
        {
            Quantitys[(int)quantityId].ContactSlaveCadIds = slaveIds;
        }

        public IList<uint> GetContactMasterCadIds(uint quantityId)
        {
            return Quantitys[(int)quantityId].ContactMasterCadIds;
        }

        public void SetContactMasterCadIds(uint quantityId, IList<uint> masterIds)
        {
            Quantitys[(int)quantityId].ContactMasterCadIds = masterIds;
        }

        public uint GetCoordCount(uint quantityId)
        {
            return Quantitys[(int)quantityId].GetCoordCount();
        }

        public double[] GetCoord(uint quantityId, int coId)
        {
            return Quantitys[(int)quantityId].GetCoord(coId, RotAngle, RotOrigin);
        }

        public uint GetNodeCount(uint quantityId)
        {
            return Quantitys[(int)quantityId].GetNodeCount();
        }

        public int GetOffset(uint quantityId)
        {
            int cnt = 0;
            for (uint tmpId = 0; tmpId < quantityId; tmpId++)
            {
                int quantityDof = (int)GetDof(tmpId);
                int quantityNodeCnt = (int)GetNodeCount(tmpId);
                cnt += quantityDof * quantityNodeCnt;
            }
            return cnt;
        }

        public int Coord2Node(uint quantityId, int coId)
        {
            return Quantitys[(int)quantityId].Coord2Node(coId);
        }

        public int Node2Coord(uint quantityId, int nodeId)
        {
            return Quantitys[(int)quantityId].Node2Coord(nodeId);
        }

        public uint GetEdgeNodeCount(uint quantityId)
        {
            return Quantitys[(int)quantityId].GetEdgeNodeCount();
        }

        public int[] GetEdgeCoordIds(uint quantityId, int edgeId)
        {
            return Quantitys[(int)quantityId].GetEdgeCoordIds(edgeId);
        }

        public int GetEdgeIdFromCoords(uint quantityId, int coId1, int coId2, out bool isReverse)
        {
            return Quantitys[(int)quantityId].GetEdgeIdFromCoords(coId1, coId2, out isReverse);
        }

        public int Edge2EdgeNode(uint quantityId, int edgeId)
        {
            return Quantitys[(int)quantityId].Edge2EdgeNode(edgeId);
        }

        public int EdgeNode2Edge(uint quantityId, int nodeId)
        {
            return Quantitys[(int)quantityId].EdgeNode2Edge(nodeId);
        }

        public void SetIncidentPortId(uint quantityId, int incidentPortId)
        {
            Quantitys[(int)quantityId].IncidentPortId = incidentPortId;
        }

        public int GetIncidentPortId(uint quantityId)
        {
            return Quantitys[(int)quantityId].IncidentPortId;
        }

        public void SetIncidentModeId(uint quantityId, int incidentModeId)
        {
            Quantitys[(int)quantityId].IncidentModeId = incidentModeId;
        }

        public int GetIncidentModeId(uint quantityId)
        {
            return Quantitys[(int)quantityId].IncidentModeId;
        }

        public IList<PortCondition> GetPortConditions(uint quantityId)
        {
            return Quantitys[(int)quantityId].PortConditions;
        }

        public uint GetPortCount(uint quantityId)
        {
            return Quantitys[(int)quantityId].GetPortCount();
        }

        public uint GetPortNodeCount(uint quantityId, uint portId)
        {
            return Quantitys[(int)quantityId].GetPortNodeCount(portId);
        }

        public IList<int> GetPortCoIds(uint quantityId, uint portId)
        {
            return Quantitys[(int)quantityId].GetPortCoIds(this, portId);
        }

        public double GetPortLineLength(uint quantityId, uint portId)
        {
            return Quantitys[(int)quantityId].GetPortLineLength(this, portId, RotAngle, RotOrigin);
        }

        public int PortCoord2Node(uint quantityId, uint portId, int coId)
        {
            return Quantitys[(int)quantityId].PortCoord2Node(portId, coId);
        }

        public int PortNode2Coord(uint quantityId, uint portId, int nodeId)
        {
            return Quantitys[(int)quantityId].PortNode2Coord(portId, nodeId);
        }

        public uint GetPortEdgeNodeCount(uint quantityId, uint portId)
        {
            return Quantitys[(int)quantityId].GetPortEdgeNodeCount(portId);
        }

        public int PortEdge2EdgeNode(uint quantityId, uint portId, int edgeId)
        {
            return Quantitys[(int)quantityId].PortEdge2EdgeNode(portId, edgeId);
        }

        public int PortEdgeNode2Edge(uint quantityId, uint portId, int nodeId)
        {
            return Quantitys[(int)quantityId].PortEdgeNode2Edge(portId, nodeId);
        }

        public IList<int> GetPeriodicPortBcCoIds(uint quantityId, uint portId, uint bcIndex)
        {
            return Quantitys[(int)quantityId].GetPeriodicPortBcCoIds(portId, bcIndex);
        }

        public uint GetDof(uint quantityId)
        {
            return Quantitys[(int)quantityId].Dof;
        }

        public IList<uint> GetMaterialIds()
        {
            return MaterialArray.GetObjectIds();
        }

        public bool IsMaterialId(uint maId)
        {
            return MaterialArray.IsObjectId(maId);
        }

        public Material GetMaterial(uint maId)
        {
            System.Diagnostics.Debug.Assert(MaterialArray.IsObjectId(maId));
            return MaterialArray.GetObject(maId);
        }

        public uint AddMaterial(Material material)
        {
            uint freeId = MaterialArray.GetFreeObjectId();
            uint maId = MaterialArray.AddObject(freeId, material);
            System.Diagnostics.Debug.Assert(maId == freeId);
            return maId;
        }

        public void ClearMaterial()
        {
            MaterialArray.Clear();
        }

        public void SetCadEdgeMaterial(uint eCadId, uint maId)
        {
            if (!CadEdge2Material.ContainsKey(eCadId))
            {
                CadEdge2Material.Add(eCadId, maId);
            }
            else
            {
                CadEdge2Material[eCadId] = maId;
            }
        }

        public uint GetCadEdgeMaterial(uint eCadId)
        {
            System.Diagnostics.Debug.Assert(CadEdge2Material.ContainsKey(eCadId));
            uint maId = CadEdge2Material[eCadId];
            return maId;
        }

        public void ClearCadEdgeMaterial()
        {
            CadEdge2Material.Clear();
        }

        public void SetCadLoopMaterial(uint lCadId, uint maId)
        {
            if (!CadLoop2Material.ContainsKey(lCadId))
            {
                CadLoop2Material.Add(lCadId, maId);
            }
            else
            {
                CadLoop2Material[lCadId] = maId;
            }
        }

        public uint GetCadLoopMaterial(uint lCadId)
        {
            System.Diagnostics.Debug.Assert(CadLoop2Material.ContainsKey(lCadId));
            uint maId = CadLoop2Material[lCadId];
            return maId;
        }

        public void ClearCadLoopMaterial()
        {
            CadLoop2Material.Clear();
        }

        public void SetCadSolidMaterial(uint sCadId, uint maId)
        {
            if (!CadSolid2Material.ContainsKey(sCadId))
            {
                CadSolid2Material.Add(sCadId, maId);
            }
            else
            {
                CadSolid2Material[sCadId] = maId;
            }
        }

        public uint GetCadSolidMaterial(uint lCadId)
        {
            System.Diagnostics.Debug.Assert(CadSolid2Material.ContainsKey(lCadId));
            uint maId = CadSolid2Material[lCadId];
            return maId;
        }

        public void ClearCadSolidMaterial()
        {
            CadSolid2Material.Clear();
        }

        public IList<int> GetCoordIdsFromCadId(uint quantityId, uint cadId, CadElementType cadElemType)
        {
            return Quantitys[(int)quantityId].GetCoordIdsFromCadId(this, cadId, cadElemType);
        }

        public IList<FieldFixedCad> GetFixedCadsFromCoord(uint quantityId, int coId)
        {
            return Quantitys[(int)quantityId].GetFixedCadsFromCoord(coId);
        }

        public IList<FieldFixedCad> GetForceFixedCadsFromCoord(uint quantityId, int coId)
        {
            return Quantitys[(int)quantityId].GetForceFixedCadsFromCoord(coId);
        }

        public IList<MultipointConstraint> GetMultipointConstraintFromCoord(uint quantityId, int coId)
        {
            return Quantitys[(int)quantityId].GetMultipointConstraintFromCoord(coId);
        }

        public uint GetLineFEIdFromMesh(uint quantityId, uint meshId, uint iElem)
        {
            return Quantitys[(int)quantityId].GetLineFEIdFromMesh(meshId, iElem);
        }

        public uint GetTriangleFEIdFromMesh(uint quantityId, uint meshId, uint iElem)
        {
            return Quantitys[(int)quantityId].GetTriangleFEIdFromMesh(meshId, iElem);
        }

        public uint GetTetrahedronFEIdFromMesh(uint quantityId, uint meshId, uint iElem)
        {
            return Quantitys[(int)quantityId].GetTetrahedronFEIdFromMesh(meshId, iElem);
        }

        public IList<uint> GetLineFEIdsFromCoord(uint quantityId, int coId)
        {
            return Quantitys[(int)quantityId].GetLineFEIdsFromCoord(coId);
        }

        public IList<uint> GetTriangleFEIdsFromCoord(uint quantityId, int coId)
        {
            return Quantitys[(int)quantityId].GetTriangleFEIdsFromCoord(coId);
        }

        public IList<uint> GetTriangleFEIdsFromEdgeCoord(uint quantityId, int coId1, int coId2)
        {
            return Quantitys[(int)quantityId].GetTriangleFEIdsFromEdgeCoord(coId1, coId2);
        }

        public IList<uint> GetTetrahedronFEIdsFromCoord(uint quantityId, int coId)
        {
            return Quantitys[(int)quantityId].GetTetrahedronFEIdsFromCoord(coId);
        }

        public IList<uint> GetLineFEIdsFromEdgeCadId(uint quantityId, uint eId)
        {
            return Quantitys[(int)quantityId].GetLineFEIdsFromEdgeCadId(this, eId);
        }

        public IList<uint> GetTriangleFEIdsFromLoopCadId(uint quantityId, uint lId)
        {
            return Quantitys[(int)quantityId].GetTriangleFEIdsFromLoopCadId(this, lId);
        }

        public IList<uint> GetTetrahedronFEIdsFromSolidCadId(uint quantityId, uint lId)
        {
            return Quantitys[(int)quantityId].GetTetrahedronFEIdsFromSolidCadId(this, lId);
        }

        public IList<uint> GetLineFEIds(uint quantityId)
        {
            return Quantitys[(int)quantityId].GetLineFEIds();
        }

        public LineFE GetLineFE(uint quantityId, uint feId)
        {
            return Quantitys[(int)quantityId].GetLineFE(feId);
        }

        public IList<uint> GetPortLineFEIds(uint quantityId, uint portId)
        {
            return Quantitys[(int)quantityId].GetPortLineFEIds(portId);
        }

        public IList<uint> GetPeriodicPortLineFEIds(uint quantityId, uint portId, uint bcIndex)
        {
            return Quantitys[(int)quantityId].GetPeriodicPortLineFEIds(portId, bcIndex);
        }

        public IList<uint> GetPeriodicPortTriangleFEIds(uint quantityId, uint portId)
        {
            return Quantitys[(int)quantityId].GetPeriodicPortTriangleFEIds(portId);
        }

        public IList<uint> GetPortTriangleFEIds(uint quantityId, uint portId)
        {
            return Quantitys[(int)quantityId].GetPortTriangleFEIds(portId);
        }

        public IList<uint> GetContactSlaveFEIds(uint quantityId)
        {
            return Quantitys[(int)quantityId].GetContactSlaveFEIds();
        }

        public IList<uint> GetContactMasterFEIds(uint quantityId)
        {
            return Quantitys[(int)quantityId].GetContactMasterFEIds();
        }

        public IList<uint> GetTriangleFEIds(uint quantityId)
        {
            return Quantitys[(int)quantityId].GetTriangleFEIds();
        }

        public TriangleFE GetTriangleFE(uint quantityId, uint feId)
        {
            return Quantitys[(int)quantityId].GetTriangleFE(feId);
        }

        public IList<uint> GetTetrahedronFEIds(uint quantityId)
        {
            return Quantitys[(int)quantityId].GetTetrahedronFEIds();
        }

        public TetrahedronFE GetTetrahedronFE(uint quantityId, uint feId)
        {
            return Quantitys[(int)quantityId].GetTetrahedronFE(feId);
        }

        public void MakeElements()
        {
            ClearElements();

            System.Diagnostics.Debug.Assert(Mesh != null);

            Mesh.GetCoords(out VertexCoords);

            if (Is2DOrPlate())
            {
                foreach (var quantity in Quantitys)
                {
                    quantity.MakeElements2D(this, VertexCoords, CadLoop2Material, CadEdge2Material);
                }
            }
            else
            {
                foreach (var quantity in Quantitys)
                {
                    quantity.MakeElements3D(this, VertexCoords, CadSolid2Material, CadLoop2Material, CadEdge2Material);
                }
            }
        }

        private bool Is2DOrPlate()
        {
            bool is2D = true;
            if (Mesh is Mesher2D)
            {
                is2D = true;
            }
            else if (Mesh is Mesher3D)
            {
                // 3Dでも梁やシェルの場合は2Dのルーチンで処理
                var mesher3D = Mesh as Mesher3D;
                Cad3D cad = mesher3D.Cad;
                IList<uint> sIds = cad.GetElementIds(CadElementType.Solid);
                is2D = sIds.Count > 0 ? false : true;
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return is2D;
        }

        public bool IsFieldValueId(uint valueId)
        {
            return FieldValueArray.IsObjectId(valueId);
        }

        public IList<uint> GetFieldValueIds()
        {
            return FieldValueArray.GetObjectIds();
        }

        public FieldValue GetFieldValue(uint valueId)
        {
            System.Diagnostics.Debug.Assert(FieldValueArray.IsObjectId(valueId));
            return FieldValueArray.GetObject(valueId);
        }

        public void ClearFieldValue()
        {
            FieldValueArray.Clear();
        }

        public uint AddFieldValue(FieldValueType fieldType, FieldDerivativeType derivativeType,
            uint quantityId, bool isBubble, FieldShowType showType)
        {
            FieldValueNodeType nodeType = isBubble ? FieldValueNodeType.Bubble : FieldValueNodeType.Node;
            uint valueId = AddFieldValue(fieldType, derivativeType, quantityId, nodeType, showType);
            return valueId;
        }

        public uint AddFieldValue(FieldValueType fieldType, FieldDerivativeType derivativeType,
            uint quantityId, FieldValueNodeType nodeType, FieldShowType showType)
        {
            uint pointCnt = 0;
            if (nodeType == FieldValueNodeType.Node)
            {
                pointCnt = GetCoordCount(quantityId);
            }
            else if (nodeType == FieldValueNodeType.Bubble)
            {
                if (Is2DOrPlate())
                {
                    // 2D
                    pointCnt = (uint)GetTriangleFEIds(quantityId).Count;
                }
                else
                {
                    // 3D
                    pointCnt = (uint)GetTetrahedronFEIds(quantityId).Count;
                }
            }
            else if (nodeType == FieldValueNodeType.ElementNode)
            {
                if (Is2DOrPlate())
                {
                    // 2D
                    IList<uint> feIds = GetTriangleFEIds(quantityId);
                    uint elemCnt = (uint)feIds.Count;
                    System.Diagnostics.Debug.Assert(elemCnt > 0);
                    TriangleFE triFE0 = GetTriangleFE(quantityId, feIds[0]);
                    uint elemNodeCnt = triFE0.NodeCount;
                    pointCnt = elemCnt * elemNodeCnt;
                }
                else
                {
                    // 3D
                    IList<uint> feIds = GetTetrahedronFEIds(quantityId);
                    uint elemCnt = (uint)feIds.Count;
                    System.Diagnostics.Debug.Assert(elemCnt > 0);
                    TetrahedronFE tetFE0 = GetTetrahedronFE(quantityId, feIds[0]);
                    uint elemNodeCnt = tetFE0.NodeCount;
                    pointCnt = elemCnt * elemNodeCnt;
                }
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            FieldValue fv = new FieldValue(quantityId, fieldType, derivativeType,
                nodeType, showType, pointCnt);

            uint freeId = FieldValueArray.GetFreeObjectId();
            uint valueId = FieldValueArray.AddObject(freeId, fv);
            System.Diagnostics.Debug.Assert(valueId == freeId);
            return valueId;
        }

        public void UpdateFieldValueValuesFromNodeValues(
            uint valueId, FieldDerivativeType dt, double[] nodeValues)
        {
            System.Diagnostics.Debug.Assert(FieldValueArray.IsObjectId(valueId));
            FieldValue fv = FieldValueArray.GetObject(valueId);
            uint quantityId = fv.QuantityId;
            uint dof = fv.Dof;
            System.Diagnostics.Debug.Assert(fv.Dof == GetDof(quantityId));
            double[] values = fv.GetDoubleValues(dt);
            uint coCnt = GetCoordCount(quantityId);
            uint offsetNode = 0;
            for (uint qId = 0; qId < quantityId; qId++)
            {
                offsetNode += GetNodeCount(qId) * GetDof(qId);
            }

            for (int coId = 0; coId < coCnt; coId++)
            {
                int nodeId = Coord2Node(quantityId, coId);
                if (nodeId == -1)
                {
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        values[coId * dof + iDof] = 0;
                    }
                }
                else
                {
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        values[coId * dof + iDof] = nodeValues[offsetNode + nodeId * dof + iDof];
                    }
                }
            }
        }

        public void UpdateFieldValueValuesFromNodeValues(
            uint valueId, FieldDerivativeType dt, System.Numerics.Complex[] nodeValues)
        {
            System.Diagnostics.Debug.Assert(FieldValueArray.IsObjectId(valueId));
            FieldValue fv = FieldValueArray.GetObject(valueId);
            uint quantityId = fv.QuantityId;
            uint dof = fv.Dof;
            System.Diagnostics.Debug.Assert(fv.Dof == GetDof(quantityId));
            System.Numerics.Complex[] values = fv.GetComplexValues(dt);
            uint coCnt = GetCoordCount(quantityId);
            uint offsetNode = 0;
            for (uint qId = 0; qId < quantityId; qId++)
            {
                offsetNode += GetNodeCount(qId) * GetDof(qId);
            }
            for (int coId = 0; coId < coCnt; coId++)
            {
                int nodeId = Coord2Node(quantityId, coId);
                if (nodeId == -1)
                {
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        values[coId * dof + iDof] = 0;
                    }
                }
                else
                {
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        values[coId * dof + iDof] = nodeValues[offsetNode + nodeId * dof + iDof];
                    }
                }
            }
        }

        public void UpdateBubbleFieldValueValuesFromNodeValues(
            uint valueId, FieldDerivativeType dt, double[] nodeValues)
        {
            if (Is2DOrPlate())
            {
                UpdateBubbleFieldValueValuesFromNodeValues2D(valueId, dt, nodeValues);
            }
            else
            {
                UpdateBubbleFieldValueValuesFromNodeValues3D(valueId, dt, nodeValues);
            }
        }

        private void UpdateBubbleFieldValueValuesFromNodeValues2D(
            uint valueId, FieldDerivativeType dt, double[] nodeValues)
        {
            System.Diagnostics.Debug.Assert(FieldValueArray.IsObjectId(valueId));
            FieldValue fv = FieldValueArray.GetObject(valueId);
            uint quantityId = fv.QuantityId;
            uint dof = fv.Dof;
            System.Diagnostics.Debug.Assert(fv.Dof == GetDof(quantityId));
            double[] values = fv.GetDoubleValues(dt);
            uint coCnt = GetCoordCount(quantityId);
            uint offsetNode = 0;
            for (uint qId = 0; qId < quantityId; qId++)
            {
                offsetNode += GetNodeCount(qId) * GetDof(qId);
            }

            IList<uint> feIds = GetTriangleFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = GetTriangleFE(quantityId, feId);
                int[] coIds = triFE.NodeCoordIds;
                uint elemNodeCnt = triFE.NodeCount;
                double[] bubbleValue = new double[dof];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = coIds[iNode];
                    int nodeId = Coord2Node(quantityId, coId);
                    if (nodeId == -1)
                    {
                        //for (int iDof = 0; iDof < dof; iDof++)
                        //{
                        //    bubbleValue[iDof] += 0;
                        //}
                    }
                    else
                    {
                        for (int iDof = 0; iDof < dof; iDof++)
                        {
                            bubbleValue[iDof] += nodeValues[offsetNode + nodeId * dof + iDof];
                        }
                    }
                }
                for (int iDof = 0; iDof < dof; iDof++)
                {
                    bubbleValue[iDof] /= elemNodeCnt; 
                }

                for (int iDof = 0; iDof < dof; iDof++)
                {
                    values[(feId - 1) * dof + iDof] = bubbleValue[iDof];
                }                
            }
        }

        private void UpdateBubbleFieldValueValuesFromNodeValues3D(
            uint valueId, FieldDerivativeType dt, double[] nodeValues)
        {
            System.Diagnostics.Debug.Assert(FieldValueArray.IsObjectId(valueId));
            FieldValue fv = FieldValueArray.GetObject(valueId);
            uint quantityId = fv.QuantityId;
            uint dof = fv.Dof;
            System.Diagnostics.Debug.Assert(fv.Dof == GetDof(quantityId));
            double[] values = fv.GetDoubleValues(dt);
            uint coCnt = GetCoordCount(quantityId);
            uint offsetNode = 0;
            for (uint qId = 0; qId < quantityId; qId++)
            {
                offsetNode += GetNodeCount(qId) * GetDof(qId);
            }

            IList<uint> feIds = GetTetrahedronFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = GetTetrahedronFE(quantityId, feId);
                int[] coIds = tetFE.NodeCoordIds;
                uint elemNodeCnt = tetFE.NodeCount;
                double[] bubbleValue = new double[dof];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = coIds[iNode];
                    int nodeId = Coord2Node(quantityId, coId);
                    if (nodeId == -1)
                    {
                        //for (int iDof = 0; iDof < dof; iDof++)
                        //{
                        //    bubbleValue[iDof] += 0;
                        //}
                    }
                    else
                    {
                        for (int iDof = 0; iDof < dof; iDof++)
                        {
                            bubbleValue[iDof] += nodeValues[offsetNode + nodeId * dof + iDof];
                        }
                    }
                }
                for (int iDof = 0; iDof < dof; iDof++)
                {
                    bubbleValue[iDof] /= elemNodeCnt;
                }

                for (int iDof = 0; iDof < dof; iDof++)
                {
                    values[(feId - 1) * dof + iDof] = bubbleValue[iDof];
                }
            }
        }

        public void UpdateBubbleFieldValueValuesFromNodeValues(
            uint valueId, FieldDerivativeType dt, System.Numerics.Complex[] nodeValues)
        {
            if (Is2DOrPlate())
            {
                UpdateBubbleFieldValueValuesFromNodeValues2D(valueId, dt, nodeValues);
            }
            else
            {
                UpdateBubbleFieldValueValuesFromNodeValues3D(valueId, dt, nodeValues);
            }
        }

        private void UpdateBubbleFieldValueValuesFromNodeValues2D(
            uint valueId, FieldDerivativeType dt, System.Numerics.Complex[] nodeValues)
        {
            System.Diagnostics.Debug.Assert(FieldValueArray.IsObjectId(valueId));
            FieldValue fv = FieldValueArray.GetObject(valueId);
            uint quantityId = fv.QuantityId;
            uint dof = fv.Dof;
            System.Diagnostics.Debug.Assert(fv.Dof == GetDof(quantityId));
            System.Numerics.Complex[] values = fv.GetComplexValues(dt);
            uint coCnt = GetCoordCount(quantityId);
            uint offsetNode = 0;
            for (uint qId = 0; qId < quantityId; qId++)
            {
                offsetNode += GetNodeCount(qId) * GetDof(qId);
            }

            IList<uint> feIds = GetTriangleFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = GetTriangleFE(quantityId, feId);
                int[] coIds = triFE.NodeCoordIds;
                uint elemNodeCnt = triFE.NodeCount;
                System.Numerics.Complex[] bubbleValue = new System.Numerics.Complex[dof];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = coIds[iNode];
                    int nodeId = Coord2Node(quantityId, coId);
                    if (nodeId == -1)
                    {
                        //for (int iDof = 0; iDof < dof; iDof++)
                        //{
                        //    bubbleValue[iDof] += 0;
                        //}
                    }
                    else
                    {
                        for (int iDof = 0; iDof < dof; iDof++)
                        {
                            bubbleValue[iDof] += nodeValues[offsetNode + nodeId * dof + iDof];
                        }
                    }
                }
                for (int iDof = 0; iDof < dof; iDof++)
                {
                    bubbleValue[iDof] /= elemNodeCnt;
                }

                for (int iDof = 0; iDof < dof; iDof++)
                {
                    values[(feId - 1) * dof + iDof] = bubbleValue[iDof];
                }
            }
        }

        private void UpdateBubbleFieldValueValuesFromNodeValues3D(
            uint valueId, FieldDerivativeType dt, System.Numerics.Complex[] nodeValues)
        {
            System.Diagnostics.Debug.Assert(FieldValueArray.IsObjectId(valueId));
            FieldValue fv = FieldValueArray.GetObject(valueId);
            uint quantityId = fv.QuantityId;
            uint dof = fv.Dof;
            System.Diagnostics.Debug.Assert(fv.Dof == GetDof(quantityId));
            System.Numerics.Complex[] values = fv.GetComplexValues(dt);
            uint coCnt = GetCoordCount(quantityId);
            uint offsetNode = 0;
            for (uint qId = 0; qId < quantityId; qId++)
            {
                offsetNode += GetNodeCount(qId) * GetDof(qId);
            }

            IList<uint> feIds = GetTetrahedronFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = GetTetrahedronFE(quantityId, feId);
                int[] coIds = tetFE.NodeCoordIds;
                uint elemNodeCnt = tetFE.NodeCount;
                System.Numerics.Complex[] bubbleValue = new System.Numerics.Complex[dof];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = coIds[iNode];
                    int nodeId = Coord2Node(quantityId, coId);
                    if (nodeId == -1)
                    {
                        //for (int iDof = 0; iDof < dof; iDof++)
                        //{
                        //    bubbleValue[iDof] += 0;
                        //}
                    }
                    else
                    {
                        for (int iDof = 0; iDof < dof; iDof++)
                        {
                            bubbleValue[iDof] += nodeValues[offsetNode + nodeId * dof + iDof];
                        }
                    }
                }
                for (int iDof = 0; iDof < dof; iDof++)
                {
                    bubbleValue[iDof] /= elemNodeCnt;
                }

                for (int iDof = 0; iDof < dof; iDof++)
                {
                    values[(feId - 1) * dof + iDof] = bubbleValue[iDof];
                }
            }
        }

        public void UpdateFieldValueValuesFromCoordValues(
            uint valueId, FieldDerivativeType dt, double[] coordValues)
        {
            System.Diagnostics.Debug.Assert(FieldValueArray.IsObjectId(valueId));
            FieldValue fv = FieldValueArray.GetObject(valueId);
            //uint quantityId = fv.QuantityId;
            //uint dof = fv.Dof;
            double[] values = fv.GetDoubleValues(dt);
            System.Diagnostics.Debug.Assert(values.Length == coordValues.Length);
            coordValues.CopyTo(values, 0);
        }

        public void UpdateFieldValueValuesFromCoordValues(
            uint valueId, FieldDerivativeType dt, System.Numerics.Complex[] coordValues)
        {
            System.Diagnostics.Debug.Assert(FieldValueArray.IsObjectId(valueId));
            FieldValue fv = FieldValueArray.GetObject(valueId);
            //uint quantityId = fv.QuantityId;
            //uint dof = fv.Dof;
            System.Numerics.Complex[] values = fv.GetComplexValues(dt);
            System.Diagnostics.Debug.Assert(values.Length == coordValues.Length);
            coordValues.CopyTo(values, 0);
        }

        public void UpdateBubbleFieldValueValuesFromCoordValues(
            uint valueId, FieldDerivativeType dt, double[] coordValues)
        {
            if (Is2DOrPlate())
            {
                UpdateBubbleFieldValueValuesFromCoordValues2D(valueId, dt, coordValues);
            }
            else
            {
                UpdateBubbleFieldValueValuesFromCoordValues3D(valueId, dt, coordValues);
            }
        }

        private void UpdateBubbleFieldValueValuesFromCoordValues2D(
            uint valueId, FieldDerivativeType dt, double[] coordValues)
        {
            System.Diagnostics.Debug.Assert(FieldValueArray.IsObjectId(valueId));
            FieldValue fv = FieldValueArray.GetObject(valueId);
            uint quantityId = fv.QuantityId;
            uint dof = fv.Dof;
            double[] values = fv.GetDoubleValues(dt);
            uint coCnt = GetCoordCount(quantityId);
            System.Diagnostics.Debug.Assert(coCnt * dof == coordValues.Length);

            IList<uint> feIds = GetTriangleFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = GetTriangleFE(quantityId, feId);
                int[] coIds = triFE.NodeCoordIds;
                uint elemNodeCnt = triFE.NodeCount;
                double[] bubbleValue = new double[dof];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = coIds[iNode];
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        bubbleValue[iDof] += coordValues[coId * dof + iDof];
                    }
                }
                for (int iDof = 0; iDof < dof; iDof++)
                {
                    bubbleValue[iDof] /= elemNodeCnt;
                }

                for (int iDof = 0; iDof < dof; iDof++)
                {
                    values[(feId - 1) * dof + iDof] = bubbleValue[iDof];
                }
            }
        }

        private void UpdateBubbleFieldValueValuesFromCoordValues3D(
            uint valueId, FieldDerivativeType dt, double[] coordValues)
        {
            System.Diagnostics.Debug.Assert(FieldValueArray.IsObjectId(valueId));
            FieldValue fv = FieldValueArray.GetObject(valueId);
            uint quantityId = fv.QuantityId;
            uint dof = fv.Dof;
            double[] values = fv.GetDoubleValues(dt);
            uint coCnt = GetCoordCount(quantityId);
            System.Diagnostics.Debug.Assert(coCnt * dof == coordValues.Length);

            IList<uint> feIds = GetTetrahedronFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = GetTetrahedronFE(quantityId, feId);
                int[] coIds = tetFE.NodeCoordIds;
                uint elemNodeCnt = tetFE.NodeCount;
                double[] bubbleValue = new double[dof];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = coIds[iNode];
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        bubbleValue[iDof] += coordValues[coId * dof + iDof];
                    }
                }
                for (int iDof = 0; iDof < dof; iDof++)
                {
                    bubbleValue[iDof] /= elemNodeCnt;
                }

                for (int iDof = 0; iDof < dof; iDof++)
                {
                    values[(feId - 1) * dof + iDof] = bubbleValue[iDof];
                }
            }
        }

        public void UpdateBubbleFieldValueValuesFromCoordValues(
            uint valueId, FieldDerivativeType dt, System.Numerics.Complex[] coordValues)
        {
            if (Is2DOrPlate())
            {
                UpdateBubbleFieldValueValuesFromCoordValues2D(valueId, dt, coordValues);
            }
            else
            {
                UpdateBubbleFieldValueValuesFromCoordValues3D(valueId, dt, coordValues);
            }
        }

        private void UpdateBubbleFieldValueValuesFromCoordValues2D(
            uint valueId, FieldDerivativeType dt, System.Numerics.Complex[] coordValues)
        {
            System.Diagnostics.Debug.Assert(FieldValueArray.IsObjectId(valueId));
            FieldValue fv = FieldValueArray.GetObject(valueId);
            uint quantityId = fv.QuantityId;
            uint dof = fv.Dof;
            System.Numerics.Complex[] values = fv.GetComplexValues(dt);
            uint coCnt = GetCoordCount(quantityId);
            System.Diagnostics.Debug.Assert(coCnt * dof == coordValues.Length);

            IList<uint> feIds = GetTriangleFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = GetTriangleFE(quantityId, feId);
                int[] coIds = triFE.NodeCoordIds;
                uint elemNodeCnt = triFE.NodeCount;
                System.Numerics.Complex[] bubbleValue = new System.Numerics.Complex[dof];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = coIds[iNode];
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        bubbleValue[iDof] += coordValues[coId * dof + iDof];
                    }
                }
                for (int iDof = 0; iDof < dof; iDof++)
                {
                    bubbleValue[iDof] /= elemNodeCnt;
                }

                for (int iDof = 0; iDof < dof; iDof++)
                {
                    values[(feId - 1) * dof + iDof] = bubbleValue[iDof];
                }
            }
        }

        private void UpdateBubbleFieldValueValuesFromCoordValues3D(
            uint valueId, FieldDerivativeType dt, System.Numerics.Complex[] coordValues)
        {
            System.Diagnostics.Debug.Assert(FieldValueArray.IsObjectId(valueId));
            FieldValue fv = FieldValueArray.GetObject(valueId);
            uint quantityId = fv.QuantityId;
            uint dof = fv.Dof;
            System.Numerics.Complex[] values = fv.GetComplexValues(dt);
            uint coCnt = GetCoordCount(quantityId);
            System.Diagnostics.Debug.Assert(coCnt * dof == coordValues.Length);

            IList<uint> feIds = GetTetrahedronFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TetrahedronFE triFE = GetTetrahedronFE(quantityId, feId);
                int[] coIds = triFE.NodeCoordIds;
                uint elemNodeCnt = triFE.NodeCount;
                System.Numerics.Complex[] bubbleValue = new System.Numerics.Complex[dof];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = coIds[iNode];
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        bubbleValue[iDof] += coordValues[coId * dof + iDof];
                    }
                }
                for (int iDof = 0; iDof < dof; iDof++)
                {
                    bubbleValue[iDof] /= elemNodeCnt;
                }

                for (int iDof = 0; iDof < dof; iDof++)
                {
                    values[(feId - 1) * dof + iDof] = bubbleValue[iDof];
                }
            }
        }

        private double[] AddDisplacement(FE fe, int iNode, double[] co)
        {
            int dim = co.Length;
            double[] curCo = new double[dim];
            co.CopyTo(curCo, 0);
            if (fe.Displacements != null)
            {
                double[] u = fe.Displacements[iNode];
                for (int iDim = 0; iDim < dim; iDim++)
                {
                    curCo[iDim] += u[iDim];
                }
            }
            return curCo;
        }

        public IList<uint> GetTriangleFEsWithPointInside(uint quantityId, double[] coord)
        {
            IList<uint> feIds = new List<uint>();
            if (coord.Length == 2)
            {
                feIds = GetTriangleFEWithPointInside2D(quantityId, coord);
            }
            else if (coord.Length == 3)
            {
                // 3D plate
                feIds = GetTriangleFEWithPointInside3D(quantityId, coord);
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return feIds;
        }

        private IList<uint> GetTriangleFEWithPointInside2D(uint quantityId, double[] coord)
        {
            System.Diagnostics.Debug.Assert(coord.Length == 2);
            IList<uint> hitFEIds = new List<uint>();

            OpenTK.Vector2d pt = new OpenTK.Vector2d(coord[0], coord[1]);
            IList<uint> feIds = GetTriangleFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = GetTriangleFE(quantityId, feId);
                int[] vertexCoIds = triFE.VertexCoordIds;
                System.Diagnostics.Debug.Assert(vertexCoIds.Length == 3);
                OpenTK.Vector2d[] v = new OpenTK.Vector2d[vertexCoIds.Length];
                for (int i = 0; i < vertexCoIds.Length; i++)
                {
                    int coId = vertexCoIds[i];
                    double[] vCo = GetCoord(quantityId, coId);
                    vCo = AddDisplacement(triFE, i, vCo);
                    v[i] = new OpenTK.Vector2d(vCo[0], vCo[1]);
                }
                bool isInside = CadUtils2D.IsPointInsideTriangle(pt, v[0], v[1], v[2]);
                if (isInside)
                {
                    hitFEIds.Add(feId);
                }
            }
            return hitFEIds;
        }

        private IList<uint> GetTriangleFEWithPointInside3D(uint quantityId, double[] coord)
        {
            System.Diagnostics.Debug.Assert(coord.Length == 3);
            IList<uint> hitFEIds = new List<uint>();

            OpenTK.Vector3d pt = new OpenTK.Vector3d(coord[0], coord[1], coord[2]);
            IList<uint> feIds = GetTriangleFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = GetTriangleFE(quantityId, feId);
                int[] vertexCoIds = triFE.VertexCoordIds;
                System.Diagnostics.Debug.Assert(vertexCoIds.Length == 3);
                OpenTK.Vector3d[] v = new OpenTK.Vector3d[vertexCoIds.Length];
                for (int i = 0; i < vertexCoIds.Length; i++)
                {
                    int coId = vertexCoIds[i];
                    double[] vCo = GetCoord(quantityId, coId);
                    vCo = AddDisplacement(triFE, i, vCo);
                    v[i] = new OpenTK.Vector3d(vCo[0], vCo[1], vCo[2]);
                }
                bool isInside = CadUtils3D.IsPointInsideTriangle(pt, v[0], v[1], v[2]);
                if (isInside)
                {
                    hitFEIds.Add(feId);
                }
            }
            return hitFEIds;
        }

        public IList<uint> GetTetrahedronFEsWithPointInside(uint quantityId, double[] coord)
        {
            System.Diagnostics.Debug.Assert(coord.Length == 3);
            IList<uint> hitFEIds = new List<uint>();

            OpenTK.Vector3d pt = new OpenTK.Vector3d(coord[0], coord[1], coord[2]);
            IList<uint> feIds = GetTetrahedronFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = GetTetrahedronFE(quantityId, feId);
                int[] vertexCoIds = tetFE.VertexCoordIds;
                System.Diagnostics.Debug.Assert(vertexCoIds.Length == 4);
                OpenTK.Vector3d[] v = new OpenTK.Vector3d[vertexCoIds.Length];
                for (int i = 0; i < vertexCoIds.Length; i++)
                {
                    int coId = vertexCoIds[i];
                    double[] vCo = GetCoord(quantityId, coId);
                    vCo = AddDisplacement(tetFE, i, vCo);
                    v[i] = new OpenTK.Vector3d(vCo[0], vCo[1], vCo[2]);
                }
                bool isInside = CadUtils3D.IsPointInsideTetrahedron(pt, v[0], v[1], v[2], v[3]);
                if (isInside)
                {
                    hitFEIds.Add(feId);
                }
            }
            return hitFEIds;
        }

        public double[] GetDoublePointValueFromNodeValues(
            uint quantityId, double[] coord, double[] nodeValues)
        {
            double[] value = null;
            if (Is2DOrPlate())
            {
                value = GetDoublePointValueFromNodeValues2D(quantityId, coord, nodeValues);
            }
            else
            {
                value = GetDoublePointValueFromNodeValues3D(quantityId, coord, nodeValues);
            }
            return value;
        }

        private double[] GetDoublePointValueFromNodeValues2D(
            uint quantityId, double[] coord, double[] nodeValues)
        {
            double[] retValue = null;
            uint dof = GetDof(quantityId);
            int offset = GetOffset(quantityId);

            IList<uint> feIds = GetTriangleFEsWithPointInside(quantityId, coord);
            if (feIds.Count == 0)
            {
                return retValue;
            }
            int cnt = 0;
            retValue = new double[dof];
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = GetTriangleFE(quantityId, feId);
                uint elemNodeCnt = triFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = triFE.NodeCoordIds[iNode];
                    int nodeId = Coord2Node(quantityId, coId);
                    nodes[iNode] = nodeId;
                }
                double[] values = new double[elemNodeCnt * dof];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int nodeId = nodes[iNode];
                    if (nodeId == -1)
                    {
                        continue;
                    }
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        values[iNode * dof + iDof] = nodeValues[nodeId * dof + iDof + offset];
                    }
                }

                double[] L = triFE.Coord2L(coord);
                double[] N = triFE.CalcN(L);

                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        retValue[iDof] += N[iNode] * values[iNode * dof + iDof];
                    }
                }
                cnt++;
            }
            for (int iDof = 0; iDof < dof; iDof++)
            {
                retValue[iDof] /= cnt;
            }

            return retValue;
        }

        private double[] GetDoublePointValueFromNodeValues3D(
            uint quantityId, double[] coord, double[] nodeValues)
        {
            double[] retValue = null;
            uint dof = GetDof(quantityId);
            int offset = GetOffset(quantityId);
 
            IList<uint> feIds = GetTetrahedronFEsWithPointInside(quantityId, coord);
            if (feIds.Count == 0)
            {
                return retValue;
            }
            int cnt = 0;
            retValue = new double[dof];
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = GetTetrahedronFE(quantityId, feId);
                uint elemNodeCnt = tetFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = tetFE.NodeCoordIds[iNode];
                    int nodeId = Coord2Node(quantityId, coId);
                    nodes[iNode] = nodeId;
                }
                double[] values = new double[elemNodeCnt * dof];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int nodeId = nodes[iNode];
                    if (nodeId == -1)
                    {
                        continue;
                    }
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        values[iNode * dof + iDof] = nodeValues[nodeId * dof + iDof + offset];
                    }
                }

                double[] L = tetFE.Coord2L(coord);
                double[] N = tetFE.CalcN(L);

                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        retValue[iDof] += N[iNode] * values[iNode * dof + iDof];
                    }
                }
                cnt++; 
            }
            for (int iDof = 0; iDof < dof; iDof++)
            {
                retValue[iDof] /= cnt;
            }

            return retValue;
        }

        public System.Numerics.Complex[] GetComplexPointValueFromNodeValues(
            uint quantityId, double[] coord, System.Numerics.Complex[] nodeValues)
        {
            System.Numerics.Complex[] value = null;
            if (Is2DOrPlate())
            {
                value = GetComplexPointValueFromNodeValues2D(quantityId, coord, nodeValues);
            }
            else
            {
                value = GetComplexPointValueFromNodeValues3D(quantityId, coord, nodeValues);
            }
            return value;
        }

        private System.Numerics.Complex[] GetComplexPointValueFromNodeValues2D(
            uint quantityId, double[] coord, System.Numerics.Complex[] nodeValues)
        {
            System.Numerics.Complex[] retValue = null;
            uint dof = GetDof(quantityId);
            int offset = GetOffset(quantityId);

            IList<uint> feIds = GetTriangleFEsWithPointInside(quantityId, coord);
            if (feIds.Count == 0)
            {
                return retValue;
            }
            int cnt = 0;
            retValue = new System.Numerics.Complex[dof];
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = GetTriangleFE(quantityId, feId);
                uint elemNodeCnt = triFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = triFE.NodeCoordIds[iNode];
                    int nodeId = Coord2Node(quantityId, coId);
                    nodes[iNode] = nodeId;
                }
                System.Numerics.Complex[] values = new System.Numerics.Complex[elemNodeCnt * dof];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int nodeId = nodes[iNode];
                    if (nodeId == -1)
                    {
                        continue;
                    }
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        values[iNode * dof + iDof] = nodeValues[nodeId * dof + iDof + offset];
                    }
                }

                double[] L = triFE.Coord2L(coord);
                double[] N = triFE.CalcN(L);

                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        retValue[iDof] += N[iNode] * values[iNode * dof + iDof];
                    }
                }
                cnt++;
            }
            for (int iDof = 0; iDof < dof; iDof++)
            {
                retValue[iDof] /= cnt;
            }

            return retValue;
        }

        private System.Numerics.Complex[] GetComplexPointValueFromNodeValues3D(
            uint quantityId, double[] coord, System.Numerics.Complex[] nodeValues)
        {
            System.Numerics.Complex[] retValue = null;
            uint dof = GetDof(quantityId);
            int offset = GetOffset(quantityId);

            IList<uint> feIds = GetTetrahedronFEsWithPointInside(quantityId, coord);
            if (feIds.Count == 0)
            {
                return retValue;
            }
            int cnt = 0;
            retValue = new System.Numerics.Complex[dof];
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = GetTetrahedronFE(quantityId, feId);
                uint elemNodeCnt = tetFE.NodeCount;
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = tetFE.NodeCoordIds[iNode];
                    int nodeId = Coord2Node(quantityId, coId);
                    nodes[iNode] = nodeId;
                }
                System.Numerics.Complex[] values = new System.Numerics.Complex[elemNodeCnt * dof];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int nodeId = nodes[iNode];
                    if (nodeId == -1)
                    {
                        continue;
                    }
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        values[iNode * dof + iDof] = nodeValues[nodeId * dof + iDof + offset];
                    }
                }

                double[] L = tetFE.Coord2L(coord);
                double[] N = tetFE.CalcN(L);

                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        retValue[iDof] += N[iNode] * values[iNode * dof + iDof];
                    }
                }
                cnt++;
            }
            for (int iDof = 0; iDof < dof; iDof++)
            {
                retValue[iDof] /= cnt;
            }

            return retValue;
        }

        public void UpdateFEDisplacements(uint quantityId, double[] U)
        {
            UpdateLineFEDisplacements(quantityId, U);
            UpdateTriangleFEDisplacements(quantityId, U);
            UpdateTetrahedronFEDisplacements(quantityId, U);
        }

        private void UpdateLineFEDisplacements(uint quantityId, double[] U)
        {
            uint quantityId0 = 0;
            int dof = (int)GetDof(quantityId0);
            IList<uint> feIds = GetLineFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                LineFE lineFE = GetLineFE(quantityId, feId);
                LineFE lineFE0 = GetLineFE(quantityId0, feId);
                uint elemNodeCnt = lineFE.NodeCount;
                uint elemNodeCnt0 = lineFE.NodeCount;
                System.Diagnostics.Debug.Assert(elemNodeCnt == elemNodeCnt0);
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = lineFE.NodeCoordIds[iNode];
                    int nodeId = Coord2Node(quantityId, coId);
                    nodes[iNode] = nodeId;
                }
                int[] nodes0 = new int[elemNodeCnt0];
                for (int iNode = 0; iNode < elemNodeCnt0; iNode++)
                {
                    int coId = lineFE0.NodeCoordIds[iNode];
                    int nodeId = Coord2Node(quantityId0, coId);
                    nodes0[iNode] = nodeId;
                }
                double[][] displacements = new double[elemNodeCnt][];
                for (int iNode = 0; iNode < elemNodeCnt0; iNode++)
                {
                    double[] u = new double[dof];
                    int nodeId0 = nodes0[iNode];
                    if (nodeId0 == -1)
                    {
                        for (int iDof = 0; iDof < dof; iDof++)
                        {
                            u[iDof] = 0;
                        }
                    }
                    else
                    {
                        for (int iDof = 0; iDof < dof; iDof++)
                        {
                            u[iDof] = U[nodeId0 * dof + iDof];
                        }
                    }
                    displacements[iNode] = u;
                }
                lineFE.SetDisplacements(displacements);
            }
        }

        private void UpdateTriangleFEDisplacements(uint quantityId, double[] U)
        {
            uint quantityId0 = 0;
            int dof = (int)GetDof(quantityId0);
            IList<uint> feIds = GetTriangleFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = GetTriangleFE(quantityId, feId);
                TriangleFE triFE0 = GetTriangleFE(quantityId0, feId);
                uint elemNodeCnt = triFE.NodeCount;
                uint elemNodeCnt0 = triFE0.NodeCount;
                System.Diagnostics.Debug.Assert(elemNodeCnt == elemNodeCnt0);
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = triFE.NodeCoordIds[iNode];
                    int nodeId = Coord2Node(quantityId, coId);
                    nodes[iNode] = nodeId;
                }
                int[] nodes0 = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt0; iNode++)
                {
                    int coId = triFE0.NodeCoordIds[iNode];
                    int nodeId = Coord2Node(quantityId0, coId);
                    nodes0[iNode] = nodeId;
                }
                double[][] displacements = new double[elemNodeCnt][];
                for (int iNode = 0; iNode < elemNodeCnt0; iNode++)
                {
                    double[] u = new double[dof];
                    int nodeId0 = nodes0[iNode];
                    if (nodeId0 == -1)
                    {
                        for (int iDof = 0; iDof < dof; iDof++)
                        {
                            u[iDof] = 0;
                        }
                    }
                    else
                    {
                        for (int iDof = 0; iDof < dof; iDof++)
                        {
                            u[iDof] = U[nodeId0 * dof + iDof];
                        }
                    }
                    displacements[iNode] = u;
                }
                triFE.SetDisplacements(displacements);
            }
        }

        private void UpdateTetrahedronFEDisplacements(uint quantityId, double[] U)
        {
            uint quantityId0 = 0;
            int dof = (int)GetDof(quantityId0);
            IList<uint> feIds = GetTetrahedronFEIds(quantityId);
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = GetTetrahedronFE(quantityId, feId);
                TetrahedronFE tetFE0 = GetTetrahedronFE(quantityId0, feId);
                uint elemNodeCnt = tetFE.NodeCount;
                uint elemNodeCnt0 = tetFE0.NodeCount;
                System.Diagnostics.Debug.Assert(elemNodeCnt == elemNodeCnt0);
                int[] nodes = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt; iNode++)
                {
                    int coId = tetFE.NodeCoordIds[iNode];
                    int nodeId = Coord2Node(quantityId, coId);
                    nodes[iNode] = nodeId;
                }
                int[] nodes0 = new int[elemNodeCnt];
                for (int iNode = 0; iNode < elemNodeCnt0; iNode++)
                {
                    int coId = tetFE0.NodeCoordIds[iNode];
                    int nodeId = Coord2Node(quantityId0, coId);
                    nodes0[iNode] = nodeId;
                }
                double[][] displacements = new double[elemNodeCnt][];
                for (int iNode = 0; iNode < elemNodeCnt0; iNode++)
                {
                    double[] u = new double[dof];
                    int nodeId0 = nodes0[iNode];
                    if (nodeId0 == -1)
                    {
                        for (int iDof = 0; iDof < dof; iDof++)
                        {
                            u[iDof] = 0;
                        }
                    }
                    else
                    {
                        for (int iDof = 0; iDof < dof; iDof++)
                        {
                            u[iDof] = U[nodeId0 * dof + iDof];
                        }
                    }
                    displacements[iNode] = u;
                }
                tetFE.SetDisplacements(displacements);
            }
        }

        public void ClearFEDisplacements(uint quantityId)
        {
            {
                IList<uint> feIds = GetLineFEIds(quantityId);
                foreach (uint feId in feIds)
                {
                    LineFE lineFE = GetLineFE(quantityId, feId);
                    lineFE.SetDisplacements(null);
                }
            }
            {
                IList<uint> feIds = GetTriangleFEIds(quantityId);
                foreach (uint feId in feIds)
                {
                    TriangleFE triFE = GetTriangleFE(quantityId, feId);
                    triFE.SetDisplacements(null);
                }
            }
            {
                IList<uint> feIds = GetTetrahedronFEIds(quantityId);
                foreach (uint feId in feIds)
                {
                    TetrahedronFE tetFE = GetTetrahedronFE(quantityId, feId);
                    tetFE.SetDisplacements(null);
                }
            }
        }

        public IList<LineFE> MakeLineElementsForDraw(uint quantityId)
        {
            return Quantitys[(int)quantityId].MakeLineElementsForDraw(this);
        }
    }
}
