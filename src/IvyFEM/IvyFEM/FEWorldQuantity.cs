using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class FEWorldQuantity
    {
        public uint Id { get; set; } = 0;
        public uint Dimension { get; set; } = 0;
        public uint Dof { get; set; } = 1;
        public uint FEOrder { get; set; } = 1;
        public int IdBaseOffset { get; set; } = 0;
        public FiniteElementType FEType { get; set; } = FiniteElementType.ScalarLagrange;
        public bool IsPortOnly { get; set; } = false; // ポート境界のみ節点番号を振る

        public IList<FieldFixedCad> ZeroFieldFixedCads { get; set; } = new List<FieldFixedCad>();
        public IList<FieldFixedCad> FieldFixedCads { get; set; } = new List<FieldFixedCad>();
        public IList<FieldFixedCad> ForceFieldFixedCads { get; set; } = new List<FieldFixedCad>();
        private Dictionary<int, IList<FieldFixedCad>> Co2FixedCads = new Dictionary<int, IList<FieldFixedCad>>();
        private Dictionary<int, IList<FieldFixedCad>> Co2ForceFixedCads = new Dictionary<int, IList<FieldFixedCad>>();
        private IList<double> Coords = new List<double>();
        private IList<string> Edges = new List<string>();
        private IList<MultipointConstraint> MultipointConstraints = new List<MultipointConstraint>();
        public IList<uint> ContactSlaveEIds { get; set; } = new List<uint>();
        public IList<uint> ContactMasterEIds { get; set; } = new List<uint>();
        private Dictionary<int, IList<MultipointConstraint>> Co2MultipointConstraints =
            new Dictionary<int, IList<MultipointConstraint>>();
        public int IncidentPortId { get; set; } = -1;
        public int IncidentModeId { get; set; } = -1;
        public IList<PortCondition> PortConditions { get; private set; } = new List<PortCondition>();
        private IList<Dictionary<int, int>> PortCo2Nodes = new List<Dictionary<int, int>>();
        private IList<Dictionary<int, int>> PortNode2Cos = new List<Dictionary<int, int>>();
        private IList<IList<IList<int>>> PeriodicPortBcCosss = new List<IList<IList<int>>>();
        private IList<IList<uint>> PortLineFEIdss = new List<IList<uint>>();
        private IList<IList<IList<uint>>> PeriodicPortLineFEIdsss = new List<IList<IList<uint>>>();
        private IList<IList<uint>> PeriodicPortTriangleFEIdss = new List<IList<uint>>();
        private IList<IList<uint>> PortTriangleFEIdss = new List<IList<uint>>();
        private Dictionary<int, int> Co2Node = new Dictionary<int, int>();
        private Dictionary<int, int> Node2Co = new Dictionary<int, int>();
        private Dictionary<int, int> Edge2ENode = new Dictionary<int, int>();
        private Dictionary<int, int> ENode2Edge = new Dictionary<int, int>();
        private Dictionary<string, uint> Mesh2LineFE = new Dictionary<string, uint>();
        private Dictionary<string, uint> Mesh2TriangleFE = new Dictionary<string, uint>();
        private Dictionary<string, uint> Mesh2TetrahedronFE = new Dictionary<string, uint>();
        private Dictionary<int, IList<uint>> Co2LineFE = new Dictionary<int, IList<uint>>();
        private Dictionary<int, IList<uint>> Co2TriangleFE = new Dictionary<int, IList<uint>>();
        private Dictionary<string, IList<uint>> EdgeCos2TriangleFE = new Dictionary<string, IList<uint>>();
        private Dictionary<int, IList<uint>> Co2TetrahedronFE = new Dictionary<int, IList<uint>>();
        private IList<uint> ContactSlaveLineFEIds = new List<uint>();
        private IList<uint> ContactMasterLineFEIds = new List<uint>();
        private ObjectArray<LineFE> LineFEArray = new ObjectArray<LineFE>();
        private ObjectArray<TriangleFE> TriangleFEArray = new ObjectArray<TriangleFE>();
        private ObjectArray<TetrahedronFE> TetrahedronFEArray = new ObjectArray<TetrahedronFE>();

        public FEWorldQuantity(
            uint id, uint dimension, uint dof, uint feOrder, FiniteElementType feType, int idBaseOffset)
        {
            Id = id;
            Dimension = dimension;
            Dof = dof;
            FEOrder = feOrder;
            FEType = feType;
            IdBaseOffset = idBaseOffset;
        }

        public void Clear()
        {
            IncidentPortId = -1;
            IncidentModeId = -1;
            PortConditions.Clear();

            ClearElements();
        }

        public void ClearElements()
        {
            Coords.Clear();
            Edges.Clear();
            Co2FixedCads.Clear();
            Co2ForceFixedCads.Clear();
            Co2MultipointConstraints.Clear();
            foreach (var portCo2Node in PortCo2Nodes)
            {
                portCo2Node.Clear();
            }
            PortCo2Nodes.Clear();
            foreach (var portNode2Co in PortNode2Cos)
            {
                portNode2Co.Clear();
            }
            PortNode2Cos.Clear();
            PeriodicPortBcCosss.Clear();
            Co2Node.Clear();
            Node2Co.Clear();
            Edge2ENode.Clear();
            ENode2Edge.Clear();
            Mesh2LineFE.Clear();
            Mesh2TriangleFE.Clear();
            Mesh2TetrahedronFE.Clear();
            Co2LineFE.Clear();
            Co2TriangleFE.Clear();
            Co2TetrahedronFE.Clear();
            EdgeCos2TriangleFE.Clear();
            ContactSlaveLineFEIds.Clear();
            ContactMasterLineFEIds.Clear();
            LineFEArray.Clear();
            PortLineFEIdss.Clear();
            PeriodicPortLineFEIdsss.Clear();
            PeriodicPortTriangleFEIdss.Clear();
            PortTriangleFEIdss.Clear();
            TriangleFEArray.Clear();
        }

        public int GetMultipointConstraintCount()
        {
            return MultipointConstraints.Count;
        }

        public int AddMultipointConstraint(MultipointConstraint mpc)
        {
            MultipointConstraints.Add(mpc);
            int index = MultipointConstraints.Count - 1;
            return index;
        }

        public MultipointConstraint GetMultipointConstraint(int index)
        {
            return MultipointConstraints[index];
        }

        public void ClearMultipointConstraint()
        {
            MultipointConstraints.Clear();
        }

        internal uint GetCoordCount()
        {
            return (uint)Coords.Count / Dimension;
        }

        public double[] GetCoord(int coId, double rotAngle, double[] rotOrigin)
        {
            double[] coord = _GetCoord(coId);
            if (Dimension == 2)
            {
                coord = CadUtils2D.GetRotCoord2D(coord, rotAngle, rotOrigin);
            }
            else if (Dimension == 3)
            {
                // TODO:
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return coord;
        }

        private double[] _GetCoord(int coId)
        {
            System.Diagnostics.Debug.Assert(coId * Dimension + (Dimension - 1) < Coords.Count);
            double[] coord = new double[Dimension];
            for (int iDim = 0; iDim < Dimension; iDim++)
            {
                coord[iDim] = Coords[(int)(coId * Dimension + iDim)];
            }
            return coord;
        }

        public uint GetNodeCount()
        {
            return (uint)Co2Node.Count;
        }

        public int Coord2Node(int coId)
        {
            if (!Co2Node.ContainsKey(coId))
            {
                return -1;
            }
            return Co2Node[coId];
        }

        public int Node2Coord(int nodeId)
        {
            if (!Node2Co.ContainsKey(nodeId))
            {
                return -1;
            }
            return Node2Co[nodeId];
        }

        public uint GetPortCount()
        {
            return (uint)PortConditions.Count;
        }

        public uint GetPortNodeCount(uint portId)
        {
            System.Diagnostics.Debug.Assert(portId < PortCo2Nodes.Count);
            return (uint)PortCo2Nodes[(int)portId].Count;
        }

        public IList<int> GetPortCoIds(FEWorld world, uint portId)
        {
            System.Diagnostics.Debug.Assert(portId < PortConditions.Count);
            IList<int> coIds = new List<int>();

            IList<uint> eIds = PortConditions[(int)portId].EIds;
            foreach (uint eId in eIds)
            {
                IList<int> tmpCoIds = GetCoordIdsFromCadId(world, eId, CadElementType.Edge);
                foreach (int tmpCoId in tmpCoIds)
                {
                    if (!coIds.Contains(tmpCoId))
                    {
                        coIds.Add(tmpCoId);
                    }
                }
            }
            return coIds;
        }

        public double GetPortLineLength(FEWorld world, uint portId, double rotAngle, double[] rotOrigin)
        {
            System.Diagnostics.Debug.Assert(portId < PortConditions.Count);
            System.Diagnostics.Debug.Assert(PortConditions[(int)portId].CadElemType == CadElementType.Edge);
            if (PortConditions[(int)portId].CadElemType != CadElementType.Edge)
            {
                return 0.0;
            }

            IList<int> coIds = GetPortCoIds(world, portId);
            // Note: 始点、終点とは限らない
            int coId1 = coIds[0];
            double[] pt1 = GetCoord(coId1, rotAngle, rotOrigin);
            int coId2 = coIds[1];
            double[] pt2 = GetCoord(coId2, rotAngle, rotOrigin);
            double[] dir = CadUtils.GetDirection(pt1, pt2);
            uint dim = Dimension;

            double minX = double.MaxValue;
            double maxX = double.MinValue;
            foreach (int coId in coIds)
            {
                double[] pt = GetCoord(coId, rotAngle, rotOrigin);
                //pt - pt1
                double[] vec = new double[(int)dim];
                for (int idim = 0; idim < dim; idim++)
                {
                    vec[idim] = pt[idim] - pt1[idim];
                }
                double X = CadUtils.Dot(vec, dir);
                if (minX > X)
                {
                    minX = X;
                }
                if (maxX < X)
                {
                    maxX = X;
                }
            }
            double length = maxX - minX;
            System.Diagnostics.Debug.Assert(length > 0);
            return length;
        }

        public int PortCoord2Node(uint portId, int coId)
        {
            var portCo2Node = PortCo2Nodes[(int)portId];
            if (!portCo2Node.ContainsKey(coId))
            {
                return -1;
            }
            return portCo2Node[coId];
        }

        public int PortNode2Coord(uint portId, int nodeId)
        {
            var portNode2Co = PortNode2Cos[(int)portId];
            if (!portNode2Co.ContainsKey(nodeId))
            {
                return -1;
            }
            return portNode2Co[nodeId];
        }

        public IList<int> GetPeriodicPortBcCoIds(uint portId, uint bcIndex)
        {
            IList<int> coIds = PeriodicPortBcCosss[(int)portId][(int)bcIndex];
            return coIds;
        }

        public uint GetEdgeCount()
        {
            return (uint)Edges.Count;
        }

        public int[] GetEdgeCoordIds(int edgeId)
        {
            string edgeKey = Edges[edgeId];
            string[] tokens = edgeKey.Split('_');
            System.Diagnostics.Debug.Assert(tokens.Length == 2);
            int coId1 = int.Parse(tokens[0]);
            int coId2 = int.Parse(tokens[1]);
            int[] edgeCoordIds = { coId1, coId2 };
            return edgeCoordIds;
        }

        public int GetEdgeIdFromCoords(int coId1, int coId2, out bool isReverse)
        {
            isReverse = (coId1 > coId2);
            string edgeKey = isReverse ?
                string.Format("{0}_{1}", coId2, coId1) :
                string.Format("{0}_{1}", coId1, coId2);
            int edgeId = Edges.IndexOf(edgeKey);
            return edgeId;
        }

        public uint GetEdgeNodeCount()
        {
            return (uint)Edge2ENode.Count;
        }

        public int Edge2EdgeNode(int edgeId)
        {
            if (!Edge2ENode.ContainsKey(edgeId))
            {
                return -1;
            }
            return Edge2ENode[edgeId];
        }

        public int EdgeNode2Edge(int edgeNodeId)
        {
            if (!ENode2Edge.ContainsKey(edgeNodeId))
            {
                return -1;
            }
            return ENode2Edge[edgeNodeId];
        }

        public IList<int> GetCoordIdsFromCadId(FEWorld world, uint cadId, CadElementType cadElemType)
        {
            var mesh = world.Mesh;
            IList<int> coIds = null;
            if (cadElemType == CadElementType.Vertex)
            {
                uint meshId = mesh.GetIdFromCadId(cadId, cadElemType);
                uint elemCnt;
                MeshType meshType;
                int loc;
                uint cadIdTmp;
                mesh.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadIdTmp);
                MeshType dummyMeshType;
                int[] vertexs;
                mesh.GetConnectivity(meshId, out dummyMeshType, out vertexs);
                System.Diagnostics.Debug.Assert(meshType == dummyMeshType);

                coIds = vertexs.ToList();
            }
            else if (cadElemType == CadElementType.Edge)
            {
                coIds = new List<int>();
                IList<uint> feIds = LineFEArray.GetObjectIds();
                foreach (uint feId in feIds)
                {
                    LineFE lineFE = LineFEArray.GetObject(feId);
                    uint cadIdTmp;
                    {
                        uint meshId = lineFE.MeshId;
                        uint elemCnt;
                        MeshType meshType;
                        int loc;
                        mesh.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadIdTmp);
                        System.Diagnostics.Debug.Assert(meshType == MeshType.Bar);
                    }
                    if (cadIdTmp == cadId)
                    {
                        foreach (int coId in lineFE.NodeCoordIds)
                        {
                            if (coIds.IndexOf(coId) == -1)
                            {
                                coIds.Add(coId);
                            }
                        }
                    }
                }
            }
            else if (cadElemType == CadElementType.Loop)
            {
                coIds = new List<int>();
                IList<uint> feIds = TriangleFEArray.GetObjectIds();
                foreach (uint feId in feIds)
                {
                    TriangleFE triFE = TriangleFEArray.GetObject(feId);
                    uint cadIdTmp;
                    {
                        uint meshId = triFE.MeshId;
                        uint elemCnt;
                        MeshType meshType;
                        int loc;
                        mesh.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadIdTmp);
                        System.Diagnostics.Debug.Assert(meshType == MeshType.Tri);
                    }
                    if (cadIdTmp == cadId)
                    {
                        foreach (int coId in triFE.NodeCoordIds)
                        {
                            if (coIds.IndexOf(coId) == -1)
                            {
                                coIds.Add(coId);
                            }
                        }
                    }
                }
            }
            else if (cadElemType == CadElementType.Solid)
            {
                coIds = new List<int>();
                IList<uint> feIds = TetrahedronFEArray.GetObjectIds();
                foreach (uint feId in feIds)
                {
                    TetrahedronFE tetFE = TetrahedronFEArray.GetObject(feId);
                    uint cadIdTmp;
                    {
                        uint meshId = tetFE.MeshId;
                        uint elemCnt;
                        MeshType meshType;
                        int loc;
                        mesh.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadIdTmp);
                        System.Diagnostics.Debug.Assert(meshType == MeshType.Tet);
                    }
                    if (cadIdTmp == cadId)
                    {
                        foreach (int coId in tetFE.NodeCoordIds)
                        {
                            if (coIds.IndexOf(coId) == -1)
                            {
                                coIds.Add(coId);
                            }
                        }
                    }
                }
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
                throw new InvalidOperationException();
            }

            return coIds;
        }

        public IList<int[]> GetEdgeCoordIdssFromCadId(FEWorld world, uint cadId)
        {
            var mesh = world.Mesh;
            IList<int[]> edgeCoIdss = null;
            {
                edgeCoIdss = new List<int[]>();
                IList<uint> feIds = LineFEArray.GetObjectIds();
                foreach (uint feId in feIds)
                {
                    LineFE lineFE = LineFEArray.GetObject(feId);
                    uint cadIdTmp;
                    {
                        uint meshId = lineFE.MeshId;
                        uint elemCnt;
                        MeshType meshType;
                        int loc;
                        mesh.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadIdTmp);
                        System.Diagnostics.Debug.Assert(meshType == MeshType.Bar);
                    }
                    if (cadIdTmp == cadId)
                    {
                        var nodeCoIds = lineFE.NodeCoordIds;
                        if (lineFE.Order == 1)
                        {
                            int coId1 = nodeCoIds[0];
                            int coId2 = nodeCoIds[1];
                            int[] edgeCoIds = coId1 < coId2 ?
                                new int[] { coId1, coId2 } :
                                new int[] { coId2, coId1 };
                            edgeCoIdss.Add(edgeCoIds);
                        }
                        else if (lineFE.Order == 2)
                        {
                            int coId1 = nodeCoIds[0];
                            int coId2 = nodeCoIds[1];
                            int coId3 = nodeCoIds[2];
                            {
                                int[] edgeCoIds = coId1 < coId3 ?
                                    new int[] { coId1, coId3 } :
                                    new int[] { coId3, coId1 };
                                edgeCoIdss.Add(edgeCoIds);
                            }
                            {
                                int[] edgeCoIds = coId3 < coId2 ?
                                    new int[] { coId3, coId2 } :
                                    new int[] { coId2, coId3 };
                                edgeCoIdss.Add(edgeCoIds);
                            }
                        }
                        else
                        {
                            // TODO:
                            System.Diagnostics.Debug.Assert(false);
                        }
                    }
                }
            }
            return edgeCoIdss;
        }

        public IList<FieldFixedCad> GetFixedCadsFromCoord(int coId)
        {
            if (!Co2FixedCads.ContainsKey(coId))
            {
                return new List<FieldFixedCad>();
            }
            return Co2FixedCads[coId];
        }

        public IList<FieldFixedCad> GetForceFixedCadsFromCoord(int coId)
        {
            if (!Co2ForceFixedCads.ContainsKey(coId))
            {
                return new List<FieldFixedCad>();
            }
            return Co2ForceFixedCads[coId];
        }

        public IList<MultipointConstraint> GetMultipointConstraintFromCoord(int coId)
        {
            if (!Co2MultipointConstraints.ContainsKey(coId))
            {
                return null;
            }
            return Co2MultipointConstraints[coId];
        }

        public uint GetLineFEIdFromMesh(uint meshId, uint iElem)
        {
            string key = meshId + "_" + iElem;
            if (Mesh2LineFE.ContainsKey(key))
            {
                uint feId = Mesh2LineFE[key];
                return feId;
            }
            return 0;
        }

        public uint GetTriangleFEIdFromMesh(uint meshId, uint iElem)
        {
            string key = meshId + "_" + iElem;
            if (Mesh2TriangleFE.ContainsKey(key))
            {
                uint feId = Mesh2TriangleFE[key];
                return feId;
            }
            return 0;
        }

        public uint GetTetrahedronFEIdFromMesh(uint meshId, uint iElem)
        {
            string key = meshId + "_" + iElem;
            if (Mesh2TetrahedronFE.ContainsKey(key))
            {
                uint feId = Mesh2TetrahedronFE[key];
                return feId;
            }
            return 0;
        }

        public IList<uint> GetLineFEIdsFromCoord(int coId)
        {
            IList<uint> feIds = new List<uint>();
            if (Co2LineFE.ContainsKey(coId))
            {
                feIds = Co2LineFE[coId];
            }
            return feIds;
        }

        public IList<uint> GetTriangleFEIdsFromCoord(int coId)
        {
            IList<uint> feIds = new List<uint>();
            if (Co2TriangleFE.ContainsKey(coId))
            {
                feIds = Co2TriangleFE[coId];
            }
            return feIds;
        }

        public IList<uint> GetTetrahedronFEIdsFromCoord(int coId)
        {
            IList<uint> feIds = new List<uint>();
            if (Co2TetrahedronFE.ContainsKey(coId))
            {
                feIds = Co2TetrahedronFE[coId];
            }
            return feIds;
        }

        public IList<uint> GetTriangleFEIdsFromEdgeCoord(int coId1, int coId2)
        {
            IList<uint> feIds = new List<uint>();
            int v1 = coId1;
            int v2 = coId2;
            if (v1 > v2)
            {
                int tmp = v1;
                v1 = v2;
                v2 = tmp;
            }
            string edgeKey = v1 + "_" + v2;

            if (EdgeCos2TriangleFE.ContainsKey(edgeKey))
            {
                feIds = EdgeCos2TriangleFE[edgeKey];
            }
            return feIds;
        }

        public IList<uint> GetLineFEIdsFromEdgeCadId(FEWorld world, uint eId)
        {
            IList<uint> retFEIds = new List<uint>();

            IList<uint> feIds = GetLineFEIds();
            foreach (uint feId in feIds)
            {
                LineFE lineFE = GetLineFE(feId);
                uint meshId = lineFE.MeshId;
                int meshElemId = lineFE.MeshElemId;
                uint cadId;
                uint elemCount;
                MeshType meshType;
                int loc;
                world.Mesh.GetMeshInfo(meshId, out elemCount, out meshType, out loc, out cadId);
                System.Diagnostics.Debug.Assert(meshType == MeshType.Bar);

                if (cadId == eId)
                {
                    retFEIds.Add(feId);
                }
            }
            return retFEIds;
        }

        public IList<uint> GetTriangleFEIdsFromLoopCadId(FEWorld world, uint lId)
        {
            IList<uint> retFEIds = new List<uint>();

            IList<uint> feIds = GetTriangleFEIds();
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = GetTriangleFE(feId);
                uint meshId = triFE.MeshId;
                int meshElemId = triFE.MeshElemId;
                uint cadId;
                uint elemCount;
                MeshType meshType;
                int loc;
                world.Mesh.GetMeshInfo(meshId, out elemCount, out meshType, out loc, out cadId);
                System.Diagnostics.Debug.Assert(meshType == MeshType.Tri);

                if (cadId == lId)
                {
                    retFEIds.Add(feId);
                }
            }
            return retFEIds;
        }

        public IList<uint> GetTetrahedronFEIdsFromSolidCadId(FEWorld world, uint sId)
        {
            IList<uint> retFEIds = new List<uint>();

            IList<uint> feIds = GetTriangleFEIds();
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = GetTetrahedronFE(feId);
                uint meshId = tetFE.MeshId;
                int meshElemId = tetFE.MeshElemId;
                uint cadId;
                uint elemCount;
                MeshType meshType;
                int loc;
                world.Mesh.GetMeshInfo(meshId, out elemCount, out meshType, out loc, out cadId);
                System.Diagnostics.Debug.Assert(meshType == MeshType.Tet);

                if (cadId == sId)
                {
                    retFEIds.Add(feId);
                }
            }
            return retFEIds;
        }

        public IList<uint> GetLineFEIds()
        {
            return LineFEArray.GetObjectIds();
        }

        public LineFE GetLineFE(uint feId)
        {
            System.Diagnostics.Debug.Assert(LineFEArray.IsObjectId(feId));
            return LineFEArray.GetObject(feId);
        }

        public IList<uint> GetPortLineFEIds(uint portId)
        {
            System.Diagnostics.Debug.Assert(portId < PortLineFEIdss.Count);
            return PortLineFEIdss[(int)portId];
        }

        public IList<uint> GetPeriodicPortLineFEIds(uint portId, uint bcIndex)
        {
            System.Diagnostics.Debug.Assert(portId < PeriodicPortLineFEIdsss.Count);
            return PeriodicPortLineFEIdsss[(int)portId][(int)bcIndex];
        }

        public IList<uint> GetPeriodicPortTriangleFEIds(uint portId)
        {
            System.Diagnostics.Debug.Assert(portId < PeriodicPortTriangleFEIdss.Count);
            return PeriodicPortTriangleFEIdss[(int)portId];
        }

        public IList<uint> GetPortTriangleFEIds(uint portId)
        {
            System.Diagnostics.Debug.Assert(portId < PortTriangleFEIdss.Count);
            return PortTriangleFEIdss[(int)portId];
        }

        public IList<uint> GetContactSlaveLineFEIds()
        {
            return ContactSlaveLineFEIds;
        }

        public IList<uint> GetContactMasterLineFEIds()
        {
            return ContactMasterLineFEIds;
        }

        public IList<uint> GetTriangleFEIds()
        {
            return TriangleFEArray.GetObjectIds();
        }

        public TriangleFE GetTriangleFE(uint feId)
        {
            System.Diagnostics.Debug.Assert(TriangleFEArray.IsObjectId(feId));
            return TriangleFEArray.GetObject(feId);
        }

        public IList<uint> GetTetrahedronFEIds()
        {
            return TetrahedronFEArray.GetObjectIds();
        }

        public TetrahedronFE GetTetrahedronFE(uint feId)
        {
            System.Diagnostics.Debug.Assert(TetrahedronFEArray.IsObjectId(feId));
            return TetrahedronFEArray.GetObject(feId);
        }

        private IList<int> GetZeroCoordIds(FEWorld world)
        {
            IList<int> zeroCoIds = new List<int>();

            foreach (var fixedCad in ZeroFieldFixedCads)
            {
                IList<int> coIds = GetCoordIdsFromCadId(world, fixedCad.CadId, fixedCad.CadElemType);
                foreach (int coId in coIds)
                {
                    if (!zeroCoIds.Contains(coId))
                    {
                        zeroCoIds.Add(coId);
                    }
                }
            }
            return zeroCoIds;
        }

        private IList<int> GetZeroEdgeIds(FEWorld world)
        {
            IList<int> zeroEdgeIds = new List<int>();

            foreach (var fixedCad in ZeroFieldFixedCads)
            {
                if (fixedCad.CadElemType != CadElementType.Edge)
                {
                    continue;
                }
                IList<int[]> edgeCoIdss = GetEdgeCoordIdssFromCadId(world, fixedCad.CadId);
                foreach (int[] edgeCoIds in edgeCoIdss)
                {
                    // 並び替えはすでに済んでいる
                    int coId1 = edgeCoIds[0];
                    int coId2 = edgeCoIds[1];
                    bool isReverse;
                    int edgeId = GetEdgeIdFromCoords(coId1, coId2, out isReverse);
                    if (!zeroEdgeIds.Contains(edgeId))
                    {
                        zeroEdgeIds.Add(edgeId);
                    }
                }
            }

            return zeroEdgeIds;
        }

        private void MakeCo2FixedCads(
            FEWorld world, IList<FieldFixedCad> fieldFixedCads, Dictionary<int, IList<FieldFixedCad>> co2FixedCads)
        {
            co2FixedCads.Clear();
            foreach (var fixedCad in fieldFixedCads)
            {
                IList<int> coIds = GetCoordIdsFromCadId(world, fixedCad.CadId, fixedCad.CadElemType);
                foreach (int coId in coIds)
                {
                    IList<FieldFixedCad> fixedCads = null;
                    if (!co2FixedCads.ContainsKey(coId))
                    {
                        fixedCads = new List<FieldFixedCad>();
                        co2FixedCads[coId] = fixedCads;
                    }
                    else
                    {
                        fixedCads = co2FixedCads[coId];
                    }
                    if (fixedCads.IndexOf(fixedCad) == -1)
                    {
                        // 同じ変数の拘束条件がすでにあるかチェック
                        {
                            IList<FieldFixedCad> sameFixedCads = new List<FieldFixedCad>();
                            foreach (FieldFixedCad tmp in fixedCads)
                            {
                                bool isSameTarget = false;
                                foreach (uint iDof in tmp.FixedDofIndexs)
                                {
                                    if (fixedCad.FixedDofIndexs.Contains(iDof))
                                    {
                                        isSameTarget = true;
                                        break;
                                    }
                                }
                                if (isSameTarget)
                                {
                                    sameFixedCads.Add(tmp);
                                }
                            }
                            if (sameFixedCads.Count > 0)
                            {
                                foreach (FieldFixedCad tmp in sameFixedCads)
                                {
                                    fixedCads.Remove(tmp);
                                }
                            }
                        }

                        fixedCads.Add(fixedCad);
                    }
                }
            }
        }

        private void SetDistributedFixedCadCoords(FEWorld world, IList<FieldFixedCad> fieldFixedCads)
        {
            foreach (var fixedCad in fieldFixedCads)
            {
                IList<int> coIds = GetCoordIdsFromCadId(world, fixedCad.CadId, fixedCad.CadElemType);
                if (fixedCad is DistributedFieldFixedCad)
                {
                    DistributedFieldFixedCad dist = fixedCad as DistributedFieldFixedCad;
                    dist.InitCoordIds(coIds);
                }
            }
        }

        private void MakeCo2MultipointConstraints(FEWorld world)
        {
            foreach (MultipointConstraint mpConstraint in MultipointConstraints)
            {
                var fixedCads = mpConstraint.FixedCads;
                foreach (var fixedCad in fixedCads)
                {
                    IList<int> coIds = GetCoordIdsFromCadId(world, fixedCad.CadId, fixedCad.CadElemType);
                    foreach (int coId in coIds)
                    {
                        IList<MultipointConstraint> mpConstraints = null;
                        if (!Co2MultipointConstraints.ContainsKey(coId))
                        {
                            Co2MultipointConstraints[coId] = new List<MultipointConstraint>();
                        }
                        mpConstraints = Co2MultipointConstraints[coId];
                        if (mpConstraints.IndexOf(mpConstraint) == -1)
                        {
                            mpConstraints.Add(mpConstraint);
                        }
                    }
                }
            }
        }

        private void MakeNode2CoFromCo2Node()
        {
            // 逆参照
            foreach (var portCo2Node in PortCo2Nodes)
            {
                var portNode2Co = new Dictionary<int, int>();
                PortNode2Cos.Add(portNode2Co);
                foreach (var pair in portCo2Node)
                {
                    int tmpPortCoId = pair.Key;
                    int tmpPortNodeId = pair.Value;
                    portNode2Co[tmpPortNodeId] = tmpPortCoId;
                }
            }
            foreach (var pair in Co2Node)
            {
                int tmpCoId = pair.Key;
                int tmpNodeId = pair.Value;
                Node2Co[tmpNodeId] = tmpCoId;
            }
        }

        private void MakeENode2EdgeFromEdge2ENode()
        {
            // 逆参照
            foreach (var pair in Edge2ENode)
            {
                int tmpEdgeId = pair.Key;
                int tmpEdgeNodeId = pair.Value;
                ENode2Edge[tmpEdgeNodeId] = tmpEdgeId;
            }
        }

        public IList<LineFE> MakeLineElementsForDraw(FEWorld world)
        {
            IList<LineFE> lineFEs = new List<LineFE>();
            HashSet<string> edges = new HashSet<string>();

            {
                var tetFEIds = GetTetrahedronFEIds();
                foreach (uint feId in tetFEIds)
                {
                    TetrahedronFE tetFE = GetTetrahedronFE(feId);
                    System.Diagnostics.Debug.Assert(tetFE.Order == FEOrder);
                    int[][] vertexCoIds =
                    {
                        new int[] { tetFE.VertexCoordIds[0], tetFE.VertexCoordIds[1] },
                        new int[] { tetFE.VertexCoordIds[1], tetFE.VertexCoordIds[2] },
                        new int[] { tetFE.VertexCoordIds[2], tetFE.VertexCoordIds[0] },
                        new int[] { tetFE.VertexCoordIds[0], tetFE.VertexCoordIds[3] },
                        new int[] { tetFE.VertexCoordIds[1], tetFE.VertexCoordIds[3] },
                        new int[] { tetFE.VertexCoordIds[2], tetFE.VertexCoordIds[3] }
                    };
                    int[][] nodeCoIds = null;
                    if (tetFE.Order == 1)
                    {
                        int[][] nodeCoIds1 =
                        {
                            new int[] { tetFE.NodeCoordIds[0], tetFE.NodeCoordIds[1] },
                            new int[] { tetFE.NodeCoordIds[1], tetFE.NodeCoordIds[2] },
                            new int[] { tetFE.NodeCoordIds[2], tetFE.NodeCoordIds[0] },
                            new int[] { tetFE.NodeCoordIds[0], tetFE.NodeCoordIds[3] },
                            new int[] { tetFE.NodeCoordIds[1], tetFE.NodeCoordIds[3] },
                            new int[] { tetFE.NodeCoordIds[2], tetFE.NodeCoordIds[3] }
                        };
                        nodeCoIds = nodeCoIds1;
                    }
                    else if (tetFE.Order == 2)
                    {
                        int[][] nodeCoIds2 =
                        {
                            new int[] { tetFE.NodeCoordIds[0], tetFE.NodeCoordIds[1], tetFE.NodeCoordIds[4] },
                            new int[] { tetFE.NodeCoordIds[1], tetFE.NodeCoordIds[2], tetFE.NodeCoordIds[5] },
                            new int[] { tetFE.NodeCoordIds[2], tetFE.NodeCoordIds[0], tetFE.NodeCoordIds[6] },
                            new int[] { tetFE.NodeCoordIds[0], tetFE.NodeCoordIds[3], tetFE.NodeCoordIds[7] },
                            new int[] { tetFE.NodeCoordIds[1], tetFE.NodeCoordIds[3], tetFE.NodeCoordIds[8] },
                            new int[] { tetFE.NodeCoordIds[2], tetFE.NodeCoordIds[3], tetFE.NodeCoordIds[9] }
                        };
                        nodeCoIds = nodeCoIds2;
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    for (int iEdge = 0; iEdge < 6; iEdge++)
                    {
                        int v1 = vertexCoIds[iEdge][0];
                        int v2 = vertexCoIds[iEdge][1];
                        if (v1 > v2)
                        {
                            int tmp = v1;
                            v1 = v2;
                            v2 = tmp;
                        }
                        string edgeKey = v1 + "_" + v2;
                        if (edges.Contains(edgeKey))
                        {
                            continue;
                        }
                        else
                        {
                            edges.Add(edgeKey);
                        }
                        var newLineFE = new LineFE((int)FEOrder, FEType);
                        newLineFE.World = world;
                        newLineFE.QuantityId = (int)this.Id;
                        newLineFE.SetVertexCoordIds(vertexCoIds[iEdge]);
                        newLineFE.SetNodeCoordIds(nodeCoIds[iEdge]);
                        // MeshId等は対応するものがないのでセットしない
                        lineFEs.Add(newLineFE);
                    }
                }
            }
            {
                var triFEIds = GetTriangleFEIds();
                foreach (uint feId in triFEIds)
                {
                    TriangleFE triFE = GetTriangleFE(feId);
                    System.Diagnostics.Debug.Assert(triFE.Order == FEOrder);
                    int[][] vertexCoIds =
                    {
                        new int[] { triFE.VertexCoordIds[0], triFE.VertexCoordIds[1] },
                        new int[] { triFE.VertexCoordIds[1], triFE.VertexCoordIds[2] },
                        new int[] { triFE.VertexCoordIds[2], triFE.VertexCoordIds[0] }
                    };
                    int[][] nodeCoIds = null;
                    if (triFE.Order == 1)
                    {
                        int[][] nodeCoIds1 =
                        {
                            new int[] { triFE.NodeCoordIds[0], triFE.NodeCoordIds[1] },
                            new int[] { triFE.NodeCoordIds[1], triFE.NodeCoordIds[2] },
                            new int[] { triFE.NodeCoordIds[2], triFE.NodeCoordIds[0] }
                        };
                        nodeCoIds = nodeCoIds1;
                    }
                    else if (triFE.Order == 2)
                    {
                        int[][] nodeCoIds2 =
                        {
                            new int[] { triFE.NodeCoordIds[0], triFE.NodeCoordIds[1], triFE.NodeCoordIds[3] },
                            new int[] { triFE.NodeCoordIds[1], triFE.NodeCoordIds[2], triFE.NodeCoordIds[4] },
                            new int[] { triFE.NodeCoordIds[2], triFE.NodeCoordIds[0], triFE.NodeCoordIds[5] }
                        };
                        nodeCoIds = nodeCoIds2;
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    for (int iEdge = 0; iEdge < 3; iEdge++)
                    {
                        int v1 = vertexCoIds[iEdge][0];
                        int v2 = vertexCoIds[iEdge][1];
                        if (v1 > v2)
                        {
                            int tmp = v1;
                            v1 = v2;
                            v2 = tmp;
                        }
                        string edgeKey = v1 + "_" + v2;
                        if (edges.Contains(edgeKey))
                        {
                            continue;
                        }
                        else
                        {
                            edges.Add(edgeKey);
                        }
                        var newLineFE = new LineFE((int)FEOrder, FEType);
                        newLineFE.World = world;
                        newLineFE.QuantityId = (int)this.Id;
                        newLineFE.SetVertexCoordIds(vertexCoIds[iEdge]);
                        newLineFE.SetNodeCoordIds(nodeCoIds[iEdge]);
                        // MeshId等は対応するものがないのでセットしない
                        lineFEs.Add(newLineFE);
                    }
                }
            }

            {
                // Loopに属さない線要素を追加
                var lineFEIds = GetLineFEIds();
                foreach (uint feId in lineFEIds)
                {
                    LineFE lineFE = GetLineFE(feId);
                    int v1 = lineFE.VertexCoordIds[0];
                    int v2 = lineFE.VertexCoordIds[1];
                    if (v1 > v2)
                    {
                        int tmp = v1;
                        v1 = v2;
                        v2 = tmp;
                    }
                    string edgeKey = v1 + "_" + v2;
                    if (edges.Contains(edgeKey))
                    {
                        continue;
                    }
                    else
                    {
                        edges.Add(edgeKey);
                    }
                    var newLineFE = new LineFE((int)FEOrder, FEType);
                    newLineFE.World = world;
                    newLineFE.QuantityId = (int)this.Id;
                    newLineFE.SetVertexCoordIds(lineFE.VertexCoordIds.ToArray());
                    newLineFE.SetNodeCoordIds(lineFE.NodeCoordIds.ToArray());
                    // MeshId等は対応するものがないのでセットしない
                    lineFEs.Add(lineFE);
                }
            }

            return lineFEs;
        }
    }
}
