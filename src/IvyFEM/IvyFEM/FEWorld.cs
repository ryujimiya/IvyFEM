using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class FEWorld
    {
        public Mesher2D Mesh { get; set; } = null;
        public uint Dimension { get; private set; } = 2;
        private IList<double> VertexCoords = new List<double>();
        private ObjectArray<Material> MaterialArray = new ObjectArray<Material>();
        private Dictionary<uint, uint> CadEdge2Material = new Dictionary<uint, uint>();
        private Dictionary<uint, uint> CadLoop2Material = new Dictionary<uint, uint>();

        public int IncidentPortId { get; set; } = -1;
        public int IncidentModeId { get; set; } = -1;
        public IList<IList<uint>> PortEIdss { get; } = new List<IList<uint>>();
        private ObjectArray<FieldValue> FieldValueArray = new ObjectArray<FieldValue>();

        public TriangleIntegrationPointCount TriIntegrationPointCount { get; set; } = TriangleIntegrationPointCount.Point3;
        private IList<FEWorldQuantity> Quantitys = new List<FEWorldQuantity>();

        public FEWorld()
        {

        }

        public uint AddQuantity(uint dof, uint feOrder)
        {
            uint id = (uint)Quantitys.Count;
            var quantity = new FEWorldQuantity(id, Dimension, dof, feOrder);
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

        public void Clear()
        {
            Mesh = null;
            MaterialArray.Clear();
            CadEdge2Material.Clear();
            CadLoop2Material.Clear();
            IncidentPortId = -1;
            IncidentModeId = -1;
            foreach (var portEIds in PortEIdss)
            {
                portEIds.Clear();
            }
            PortEIdss.Clear();
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

        public IList<FieldFixedCad> GetFieldFixedCads(uint quantityId)
        {
            return Quantitys[(int)quantityId].FieldFixedCads;
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

        public IList<uint> GetContactSlaveEIds(uint quantityId)
        {
            return Quantitys[(int)quantityId].ContactSlaveEIds;
        }

        public IList<uint> GetContactMasterEIds(uint quantityId)
        {
            return Quantitys[(int)quantityId].ContactMasterEIds;
        }

        public uint GetCoordCount(uint quantityId)
        {
            return Quantitys[(int)quantityId].GetCoordCount();
        }

        public double[] GetCoord(uint quantityId, int coId)
        {
            return Quantitys[(int)quantityId].GetCoord(coId);
        }

        public uint GetNodeCount(uint quantityId)
        {
            return Quantitys[(int)quantityId].GetNodeCount();
        }

        public int Coord2Node(uint quantityId, int coId)
        {
            return Quantitys[(int)quantityId].Coord2Node(coId);
        }

        public int Node2Coord(uint quantityId, int nodeId)
        {
            return Quantitys[(int)quantityId].Node2Coord(nodeId);
        }

        public uint GetPortCount()
        {
            return (uint)PortEIdss.Count;
        }

        public uint GetPortNodeCount(uint quantityId, uint portId)
        {
            return Quantitys[(int)quantityId].GetPortNodeCount(portId);
        }

        public int PortCoord2Node(uint quantityId, uint portId, int coId)
        {
            return Quantitys[(int)quantityId].PortCoord2Node(portId, coId);
        }

        public int PortNode2Coord(uint quantityId, uint portId, int nodeId)
        {
            return Quantitys[(int)quantityId].PortNode2Coord(portId, nodeId);
        }

        public uint GetDof(uint quantityId)
        {
            return Quantitys[(int)quantityId].Dof;
        }

        public IList<uint> GetMaterialIds()
        {
            return MaterialArray.GetObjectIds();
        }

        public Material GetMaterial(uint maId)
        {
            System.Diagnostics.Debug.Assert(MaterialArray.IsObjectId(maId));
            return MaterialArray.GetObject(maId);
        }

        public uint AddMaterial(Material material)
        {
            uint freeId = MaterialArray.GetFreeObjectId();
            uint maId = MaterialArray.AddObject(new KeyValuePair<uint, Material>(freeId, material));
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

        public IList<int> GetCoordIdsFromCadId(uint quantityId, uint cadId, CadElementType cadElemType)
        {
            return Quantitys[(int)quantityId].GetCoordIdsFromCadId(this, cadId, cadElemType);
        }

        public IList<FieldFixedCad> GetFixedCadsFromCoord(uint quantityId, int coId)
        {
            return Quantitys[(int)quantityId].GetFixedCadsFromCoord(coId);
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

        public IList<uint> GetContactSlaveLineFEIds(uint quantityId)
        {
            return Quantitys[(int)quantityId].GetContactSlaveLineFEIds();
        }

        public IList<uint> GetContactMasterLineFEIds(uint quantityId)
        {
            return Quantitys[(int)quantityId].GetContactMasterLineFEIds();
        }

        public IList<uint> GetTriangleFEIds(uint quantityId)
        {
            return Quantitys[(int)quantityId].GetTriangleFEIds();
        }

        public TriangleFE GetTriangleFE(uint quantityId, uint feId)
        {
            return Quantitys[(int)quantityId].GetTriangleFE(feId);
        }

        public void MakeElements()
        {
            ClearElements();

            System.Diagnostics.Debug.Assert(Mesh != null);

            Mesh.GetCoords(out VertexCoords);

            foreach (var quantity in Quantitys)
            {
                quantity.MakeElements(this, VertexCoords, CadLoop2Material, CadEdge2Material);
            }
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
            uint quantityId, uint dof, bool isBubble, FieldShowType showType)
        {
            uint pointCnt = 0;
            if (isBubble)
            {
                pointCnt = (uint)GetTriangleFEIds(quantityId).Count;
            }
            else
            {
                pointCnt = GetCoordCount(quantityId);
            }
            FieldValue fv = new FieldValue();
            fv.Type = fieldType;
            fv.DerivativeType = derivativeType;
            fv.IsBubble = isBubble;
            fv.ShowType = showType;
            fv.QuantityId = quantityId;
            fv.AllocValues(dof, pointCnt);

            uint freeId = FieldValueArray.GetFreeObjectId();
            uint valueId = FieldValueArray.AddObject(new KeyValuePair<uint, FieldValue>(freeId, fv));
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

        public IList<LineFE> MakeBoundOfElements(uint quantityId)
        {
            return Quantitys[(int)quantityId].MakeBoundOfElements(this);
        }

    }
}
