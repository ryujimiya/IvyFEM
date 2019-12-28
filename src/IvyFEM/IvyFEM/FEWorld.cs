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

        private ObjectArray<FieldValue> FieldValueArray = new ObjectArray<FieldValue>();

        public TriangleIntegrationPointCount TriIntegrationPointCount { get; set; } =
            TriangleIntegrationPointCount.Point7;
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

        public void Clear()
        {
            Mesh = null;
            MaterialArray.Clear();
            CadEdge2Material.Clear();
            CadLoop2Material.Clear();
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
            coord = GetRotCoord(coord, RotAngle, RotOrigin);
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

        /// <summary>
        /// 回転移動する
        /// </summary>
        /// <param name="srcPt"></param>
        /// <param name="rotAngle"></param>
        /// <param name="rotOrigin"></param>
        /// <returns></returns>
        public static double[] GetRotCoord(double[] srcPt, double rotAngle, double[] rotOrigin = null)
        {
            double[] destPt = new double[2];
            double x0 = 0;
            double y0 = 0;
            if (rotOrigin != null)
            {
                x0 = rotOrigin[0];
                y0 = rotOrigin[1];
            }
            double cosA = Math.Cos(rotAngle);
            double sinA = Math.Sin(rotAngle);
            destPt[0] = cosA * (srcPt[0] - x0) + sinA * (srcPt[1] - y0);
            destPt[1] = -sinA * (srcPt[0] - x0) + cosA * (srcPt[1] - y0);
            return destPt;
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
            return Quantitys[(int)quantityId].GetCoord(coId, RotAngle, RotOrigin);
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

        public IList<uint> GetLineFEIdsFromEdgeCadId(uint quantityId, uint eId)
        {
            return Quantitys[(int)quantityId].GetLineFEIdsFromEdgeCadId(this, eId);
        }

        public IList<uint> GetTriangleFEIdsFromLoopCadId(uint quantityId, uint lId)
        {
            return Quantitys[(int)quantityId].GetTriangleFEIdsFromLoopCadId(this, lId);
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
            uint quantityId, bool isBubble, FieldShowType showType)
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
            FieldValue fv = new FieldValue(quantityId, fieldType, derivativeType,
                isBubble, showType, pointCnt);

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

        public void UpdateBubbleFieldValueValuesFromNodeValues(
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

        public void UpdateBubbleFieldValueValuesFromCoordValues(
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

        public IList<LineFE> MakeBoundOfElements(uint quantityId)
        {
            return Quantitys[(int)quantityId].MakeBoundOfElements(this);
        }

    }
}
