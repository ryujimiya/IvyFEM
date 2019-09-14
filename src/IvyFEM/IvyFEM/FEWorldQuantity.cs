using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class FEWorldQuantity
    {
        public uint Id { get; set; } = 0;
        public uint Dimension { get; set; } = 0;
        public uint Dof { get; set; } = 1;
        public uint FEOrder { get; set; } = 1;
        public int IdBaseOffset { get; set; } = 0;
        public FiniteElementType FEType { get; set; } = FiniteElementType.ScalarLagrange;
        public IList<FieldFixedCad> ZeroFieldFixedCads { get; private set; } = new List<FieldFixedCad>();
        public IList<FieldFixedCad> FieldFixedCads { get; private set; } = new List<FieldFixedCad>();
        private Dictionary<int, IList<FieldFixedCad>> Co2FixedCads = new Dictionary<int, IList<FieldFixedCad>>();
        private IList<double> Coords = new List<double>();
        private IList<MultipointConstraint> MultipointConstraints = new List<MultipointConstraint>();
        public IList<uint> ContactSlaveEIds { get; set; } = new List<uint>();
        public IList<uint> ContactMasterEIds { get; set; } = new List<uint>();
        private Dictionary<int, IList<MultipointConstraint>> Co2MultipointConstraints =
            new Dictionary<int, IList<MultipointConstraint>>();
        public int IncidentPortId { get; set; } = -1;
        public int IncidentModeId { get; set; } = -1;
        public IList<PortCondition> PortConditions { get; private set; } =new List<PortCondition>();
        private IList<Dictionary<int, int>> PortCo2Nodes = new List<Dictionary<int, int>>();
        private IList<Dictionary<int, int>> PortNode2Cos = new List<Dictionary<int, int>>();
        private IList<IList<Dictionary<int, int>>> PeriodicPortBcCo2Nodess = new List<IList<Dictionary<int, int>>>();
        private IList<IList<Dictionary<int, int>>> PeriodicPortBcNode2Coss = new List<IList<Dictionary<int, int>>>();
        private IList<IList<uint>> PortLineFEIdss = new List<IList<uint>>();
        private IList<IList<IList<uint>>> PeriodicPortLineFEIdsss = new List<IList<IList<uint>>>();
        private IList<IList<uint>> PeriodicPortTriangleFEIdss = new List<IList<uint>>();
        private Dictionary<int, int> Co2Node = new Dictionary<int, int>();
        private Dictionary<int, int> Node2Co = new Dictionary<int, int>();
        private Dictionary<string, uint> Mesh2LineFE = new Dictionary<string, uint>();
        private Dictionary<string, uint> Mesh2TriangleFE = new Dictionary<string, uint>();
        private Dictionary<int, IList<uint>> Co2LineFE = new Dictionary<int, IList<uint>>();
        private Dictionary<int, IList<uint>> Co2TriangleFE = new Dictionary<int, IList<uint>>();
        private Dictionary<string, IList<uint>> EdgeCos2TriangleFE = new Dictionary<string, IList<uint>>();
        private IList<uint> ContactSlaveLineFEIds = new List<uint>();
        private IList<uint> ContactMasterLineFEIds = new List<uint>();
        private ObjectArray<LineFE> LineFEArray = new ObjectArray<LineFE>();
        private ObjectArray<TriangleFE> TriangleFEArray = new ObjectArray<TriangleFE>();

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
            Co2FixedCads.Clear();
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
            PeriodicPortBcCo2Nodess.Clear();
            PeriodicPortBcNode2Coss.Clear();
            Co2Node.Clear();
            Node2Co.Clear();
            Mesh2LineFE.Clear();
            Mesh2TriangleFE.Clear();
            Co2TriangleFE.Clear();
            Co2LineFE.Clear();
            EdgeCos2TriangleFE.Clear();
            ContactSlaveLineFEIds.Clear();
            ContactMasterLineFEIds.Clear();
            LineFEArray.Clear();
            foreach (var portLineFEIds in PortLineFEIdss)
            {
                portLineFEIds.Clear();
            }
            PortLineFEIdss.Clear();
            PeriodicPortLineFEIdsss.Clear();
            PeriodicPortTriangleFEIdss.Clear();
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
            coord = FEWorld.GetRotCoord(coord, rotAngle, rotOrigin);
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
            var periodicPortBcCo2Node = PeriodicPortBcCo2Nodess[(int)portId][(int)bcIndex];
            IList<int> coIds = periodicPortBcCo2Node.Keys.ToList();
            return coIds;
        }

        public int PeriodicPortBcCoord2Node(uint portId, uint bcIndex, int coId)
        {
            var periodicPortBcCo2Node = PeriodicPortBcCo2Nodess[(int)portId][(int)bcIndex];
            if (!periodicPortBcCo2Node.ContainsKey(coId))
            {
                return -1;
            }
            return periodicPortBcCo2Node[coId];
        }

        public int PeriodicPortBcNode2Coord(uint portId, uint bcIndex, int nodeId)
        {
            var periodicPortBcNode2Co = PeriodicPortBcNode2Coss[(int)portId][(int)bcIndex];
            if (!periodicPortBcNode2Co.ContainsKey(nodeId))
            {
                return -1;
            }
            return periodicPortBcNode2Co[nodeId];
        }

        public IList<int> GetCoordIdsFromCadId(FEWorld world, uint cadId, CadElementType cadElemType)
        {
            Mesher2D mesh = world.Mesh;
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
            else
            {
                throw new InvalidOperationException();
            }

            return coIds;
        }

        public IList<FieldFixedCad> GetFixedCadsFromCoord(int coId)
        {
            if (!Co2FixedCads.ContainsKey(coId))
            {
                return new List<FieldFixedCad>();
            }
            return Co2FixedCads[coId];
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
                uint meshId= lineFE.MeshId;
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

        public void MakeElements(
            FEWorld world,
            IList<double> vertexCoords,
            Dictionary<uint, uint> cadLoop2Material,
            Dictionary<uint, uint> cadEdge2Material)
        {
            ClearElements();

            // 座標、三角形要素と線要素を生成する
            MakeCoordsAndElements(
                world, vertexCoords, cadLoop2Material, cadEdge2Material);

            // 多点拘束の座標生成
            MakeCo2MultipointConstraints(world);

            // Note: 線要素生成後でないと点を特定できない
            IList<int> zeroCoordIds = GetZeroCoordIds(world);

            // 多点拘束の対象外の節点を除外する
            if (Co2MultipointConstraints.Count > 0)
            {
                for (int coId = 0; coId < Coords.Count; coId++)
                {
                    if (!Co2MultipointConstraints.ContainsKey(coId))
                    {
                        if (zeroCoordIds.IndexOf(coId) == -1)
                        {
                            zeroCoordIds.Add(coId);
                        }
                    }
                }
            }

            // Note: 三角形要素生成後でないと特定できない
            MakeCo2FixedCads(world);
            SetDistributedFixedCadCoords(world);

            // 頂点→三角形要素のマップと辺→三角形要素のマップ作成
            MakeCo2AndEdgeCos2TriangleFE();

            // ポート上の線要素の節点ナンバリング
            NumberPortNodes(world, zeroCoordIds);
            SetDistributedPortCoords(world);

            // 三角形要素の節点ナンバリング
            NumberTriangleNodes(world, zeroCoordIds);

            // 接触解析のMaster/Slave線要素を準備する
            SetupContactMasterSlaveLineElements(world);

            // 節点→座標のマップ作成
            MakeNode2CoFromCo2Node();

            // 頂点→線要素のマップ
            MakeCo2LineFE();
        }

        // 座標、三角形要素と線要素を生成する
        private void MakeCoordsAndElements(
            FEWorld world,
            IList<double> vertexCoords,
            Dictionary<uint, uint> cadLoop2Material,
            Dictionary<uint, uint> cadEdge2Material)
        {
            Mesher2D mesh = world.Mesh;
            System.Diagnostics.Debug.Assert(mesh != null);

            if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 1)
            {
                Coords = new List<double>(vertexCoords);
            }
            else if (FEType == FiniteElementType.ScalarBell && FEOrder == 5)
            {
                Coords = new List<double>(vertexCoords);
            }
            else if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 2)
            {
                Coords = new List<double>(vertexCoords);
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }

            IList<uint> meshIds = mesh.GetIds();

            //////////////////////////////////////////////////
            // 領域の三角形要素
            // まず要素を作る
            // この順番で生成した要素は隣接していない
            Dictionary<string, IList<int>> edge2MidPt = new Dictionary<string, IList<int>>();
            foreach (uint meshId in meshIds)
            {
                uint elemCnt;
                MeshType meshType;
                int loc;
                uint cadId;
                mesh.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadId);
                if (meshType != MeshType.Tri)
                {
                    continue;
                }

                if (!cadLoop2Material.ContainsKey(cadId))
                {
                    throw new IndexOutOfRangeException();
                }
                uint maId = cadLoop2Material[cadId];

                int elemVertexCnt = 3;
                int elemNodeCnt = 0;
                if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 1)
                {
                    elemNodeCnt = 3;
                }
                else if (FEType == FiniteElementType.ScalarBell && FEOrder == 5)
                {
                    elemNodeCnt = 3;
                }
                else if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 2)
                {
                    elemNodeCnt = 6;
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
                MeshType dummyMeshType;
                int[] vertexs;
                mesh.GetConnectivity(meshId, out dummyMeshType, out vertexs);
                System.Diagnostics.Debug.Assert(meshType == dummyMeshType);
                System.Diagnostics.Debug.Assert(elemVertexCnt * elemCnt == vertexs.Length);

                for (int iElem = 0; iElem < elemCnt; iElem++)
                {
                    int[] vertexCoIds = new int[elemVertexCnt];
                    for (int iPt = 0; iPt < elemVertexCnt; iPt++)
                    {
                        int coId = vertexs[iElem * elemVertexCnt + iPt];
                        vertexCoIds[iPt] = coId;
                    }
                    int[] nodeCoIds = new int[elemNodeCnt];
                    if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 1)
                    {
                        System.Diagnostics.Debug.Assert(nodeCoIds.Length == vertexCoIds.Length);
                        vertexCoIds.CopyTo(nodeCoIds, 0);
                    }
                    else if (FEType == FiniteElementType.ScalarBell && FEOrder == 5)
                    {
                        System.Diagnostics.Debug.Assert(nodeCoIds.Length == vertexCoIds.Length);
                        vertexCoIds.CopyTo(nodeCoIds, 0);
                    }
                    else if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 2)
                    {
                        for (int i = 0; i < elemVertexCnt; i++)
                        {
                            nodeCoIds[i] = vertexCoIds[i];

                            {
                                int v1 = vertexCoIds[i];
                                int v2 = vertexCoIds[(i + 1) % elemVertexCnt];
                                if (v1 > v2)
                                {
                                    int tmp = v1;
                                    v1 = v2;
                                    v2 = tmp;
                                }
                                string edgeKey = v1 + "_" + v2;
                                int midPtCoId = -1;
                                if (edge2MidPt.ContainsKey(edgeKey))
                                {
                                    midPtCoId = edge2MidPt[edgeKey][0];
                                }
                                else
                                {
                                    double[] vPt1 = world.GetVertexCoord(v1);
                                    double[] vPt2 = world.GetVertexCoord(v2);
                                    double[] midPt = { (vPt1[0] + vPt2[0]) / 2.0, (vPt1[1] + vPt2[1]) / 2.0 };
                                    midPtCoId = (int)(Coords.Count / Dimension);
                                    Coords.Add(midPt[0]);
                                    Coords.Add(midPt[1]);
                                    var list = new List<int>();
                                    list.Add(midPtCoId);
                                    edge2MidPt[edgeKey] = list;
                                }

                                nodeCoIds[i + elemVertexCnt] = midPtCoId;
                            }
                        }
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }

                    TriangleFE fe = new TriangleFE((int)FEOrder, FEType);
                    fe.World = world;
                    fe.QuantityId = (int)this.Id;
                    fe.QuantityIdBaseOffset = IdBaseOffset;
                    fe.SetVertexCoordIds(vertexCoIds);
                    fe.SetNodeCoordIds(nodeCoIds);
                    fe.MaterialId = maId;
                    fe.MeshId = meshId;
                    fe.MeshElemId = iElem;
                    // 仮登録
                    uint freeId = TriangleFEArray.GetFreeObjectId();
                    uint feId = TriangleFEArray.AddObject(freeId, fe);
                    System.Diagnostics.Debug.Assert(feId == freeId);
                }
            }

            //////////////////////////////////////////////////
            // 境界の線要素
            foreach (uint meshId in meshIds)
            {
                uint elemCnt;
                MeshType meshType;
                int loc;
                uint cadId;
                mesh.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadId);
                if (meshType != MeshType.Bar)
                {
                    continue;
                }

                int elemVertexCnt = 2;
                int elemNodeCnt = 0;
                if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 1)
                {
                    elemNodeCnt = 2;
                }
                else if (FEType == FiniteElementType.ScalarBell && FEOrder == 5)
                {
                    // 暫定：Lagrange線要素で代用
                    elemNodeCnt = 2;
                }
                else if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 2)
                {
                    elemNodeCnt = 3;
                }
                else
                {
                    System.Diagnostics.Debug.Assert(false);
                }
                MeshType dummyMeshType;
                int[] vertexs;
                mesh.GetConnectivity(meshId, out dummyMeshType, out vertexs);
                System.Diagnostics.Debug.Assert(meshType == dummyMeshType);
                System.Diagnostics.Debug.Assert(elemVertexCnt * elemCnt == vertexs.Length);

                //System.Diagnostics.Debug.Assert(CadEdge2Material.ContainsKey(cadId));
                //if (!CadEdge2Material.ContainsKey(cadId))
                //{
                //    throw new IndexOutOfRangeException();
                //}
                // 未指定のマテリアルも許容する
                uint maId = cadEdge2Material.ContainsKey(cadId) ? cadEdge2Material[cadId] : 0;

                for (int iElem = 0; iElem < elemCnt; iElem++)
                {
                    int[] vertexCoIds = new int[elemVertexCnt];
                    for (int iPt = 0; iPt < elemVertexCnt; iPt++)
                    {
                        int coId = vertexs[iElem * elemVertexCnt + iPt];
                        vertexCoIds[iPt] = coId;
                    }
                    int[] nodeCoIds = new int[elemNodeCnt];
                    if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 1)
                    {
                        System.Diagnostics.Debug.Assert(nodeCoIds.Length == vertexCoIds.Length);
                        vertexCoIds.CopyTo(nodeCoIds, 0);
                    }
                    else if (FEType == FiniteElementType.ScalarBell && FEOrder == 5)
                    {
                        // 暫定：Lagrange線要素で代用
                        System.Diagnostics.Debug.Assert(nodeCoIds.Length == vertexCoIds.Length);
                        vertexCoIds.CopyTo(nodeCoIds, 0);
                    }
                    else if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 2)
                    {
                        for (int i = 0; i < 2; i++)
                        {
                            nodeCoIds[i] = vertexCoIds[i];
                        }
                        // 線要素上の中点
                        int v1 = vertexCoIds[0];
                        int v2 = vertexCoIds[1];
                        if (v1 > v2)
                        {
                            int tmp = v1;
                            v1 = v2;
                            v2 = tmp;
                        }
                        string edgeKey = v1 + "_" + v2;
                        if (!edge2MidPt.ContainsKey(edgeKey))
                        {
                            System.Diagnostics.Debug.Assert(false);
                        }
                        int midPtCoId = edge2MidPt[edgeKey][0];
                        nodeCoIds[2] = midPtCoId;
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }

                    LineFE fe = new LineFE((int)FEOrder, FEType);
                    fe.World = world;
                    fe.QuantityId = (int)this.Id;
                    fe.SetVertexCoordIds(vertexCoIds);
                    fe.SetNodeCoordIds(nodeCoIds);
                    fe.MaterialId = maId;
                    fe.MeshId = meshId;
                    fe.MeshElemId = iElem;
                    uint freeId = LineFEArray.GetFreeObjectId();
                    uint feId = LineFEArray.AddObject(freeId, fe);
                    System.Diagnostics.Debug.Assert(feId == freeId);

                    string key = string.Format(meshId + "_" + iElem);
                    Mesh2LineFE.Add(key, feId);
                }
            }
        }

        // ポート上の線要素の節点ナンバリング
        private void NumberPortNodes(FEWorld world, IList<int> zeroCoordIds)
        {
            Mesher2D mesh = world.Mesh;

            // ポート上の線要素の抽出と節点ナンバリング
            uint portCnt = GetPortCount();
            for (int portId = 0; portId < portCnt; portId++)
            {
                PortCondition portCondition = PortConditions[portId];
                if (portCondition.IsPeriodic)
                {
                    NumberPeriodicPortNodes(world, mesh, portCondition, zeroCoordIds);
                }
                else
                {
                    NumberNormalPortNodes(world, mesh, portCondition, zeroCoordIds);
                }
            }
        }

        // 通常のポート（境界のみ）
        private void NumberNormalPortNodes(
            FEWorld world, Mesher2D mesh, PortCondition portCondition, IList<int> zeroCoordIds)
        {
            System.Diagnostics.Debug.Assert(!portCondition.IsPeriodic);
            IList<uint> portEIds = portCondition.EIds;
            var lineFEIds = new List<uint>();
            PortLineFEIdss.Add(lineFEIds);

            var portCo2Node = new Dictionary<int, int>();
            PortCo2Nodes.Add(portCo2Node);
            int portNodeId = 0;

            IList<uint> feIds = LineFEArray.GetObjectIds();
            IList<int> portCoIds = new List<int>();
            foreach (var feId in feIds)
            {
                LineFE lineFE = LineFEArray.GetObject(feId);
                uint cadId;
                {
                    uint meshId = lineFE.MeshId;
                    uint elemCnt;
                    MeshType meshType;
                    int loc;
                    mesh.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadId);
                    System.Diagnostics.Debug.Assert(meshType == MeshType.Bar);
                }
                if (portEIds.Contains(cadId))
                {
                    // ポート上の線要素
                    lineFEIds.Add(feId);

                    int[] coIds = lineFE.NodeCoordIds;

                    foreach (int coId in coIds)
                    {
                        if (portCoIds.IndexOf(coId) == -1)
                        {
                            portCoIds.Add(coId);
                        }
                    }
                }
            }
            // 境界の方向順に節点番号を振る
            uint eId0 = portEIds[0];
            IList<int> sortedCoIds;
            SortPortCoIds(world, mesh, eId0, portCoIds, out sortedCoIds);
            foreach (int coId in sortedCoIds)
            {
                if (!portCo2Node.ContainsKey(coId) &&
                    zeroCoordIds.IndexOf(coId) == -1)
                {
                    portCo2Node[coId] = portNodeId;
                    portNodeId++;
                }
            }
        }

        // 周期構造のポート（2つの境界と内部領域）
        private void NumberPeriodicPortNodes(
            FEWorld world, Mesher2D mesh, PortCondition portCondition, IList<int> zeroCoordIds)
        {
            System.Diagnostics.Debug.Assert(portCondition.IsPeriodic);

            var portCo2Node = new Dictionary<int, int>();
            PortCo2Nodes.Add(portCo2Node);
            int portNodeId = 0;

            // 境界1、2、内部の順に番号を振る
            // 2つの境界
            IList<IList<uint>> portLineFEIdss = new List<IList<uint>>();
            PeriodicPortLineFEIdsss.Add(portLineFEIdss);
            IList<Dictionary<int, int>> bcCo2Nodes = new List<Dictionary<int, int>>();
            PeriodicPortBcCo2Nodess.Add(bcCo2Nodes);
            IList<uint>[] portEIdss = { portCondition.BcEIdsForPeriodic1, portCondition.BcEIdsForPeriodic2 };
            for (int bcIndex = 0; bcIndex < 2; bcIndex++)
            {
                IList<uint> portEIds = portEIdss[bcIndex];
                Dictionary<int, int> bcCo2Node = new Dictionary<int, int>();
                bcCo2Nodes.Add(bcCo2Node);
                var lineFEIds = new List<uint>();
                portLineFEIdss.Add(lineFEIds);

                IList<int> bcCoIds = new List<int>();
                IList<uint> feIds = LineFEArray.GetObjectIds();
                foreach (var feId in feIds)
                {
                    LineFE lineFE = LineFEArray.GetObject(feId);
                    uint cadId;
                    {
                        uint meshId = lineFE.MeshId;
                        uint elemCnt;
                        MeshType meshType;
                        int loc;
                        mesh.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadId);
                        System.Diagnostics.Debug.Assert(meshType == MeshType.Bar);
                    }
                    if (portEIds.Contains(cadId))
                    {
                        // ポート上の線要素
                        lineFEIds.Add(feId);

                        int[] coIds = lineFE.NodeCoordIds;
                        foreach (int coId in coIds)
                        {
                            if (bcCoIds.IndexOf(coId) == -1)
                            {
                                bcCoIds.Add(coId);
                            }
                        }
                    }
                }
                // 境界の方向順に節点番号を振る
                uint eId0 = portEIds[0];
                IList<int> sortedCoIds;
                SortPortCoIds(world, mesh, eId0, bcCoIds, out sortedCoIds);
                foreach (int coId in sortedCoIds)
                {
                    if (!portCo2Node.ContainsKey(coId) &&
                        zeroCoordIds.IndexOf(coId) == -1)
                    {
                        portCo2Node[coId] = portNodeId;
                        bcCo2Node[coId] = portNodeId;
                        portNodeId++;
                    }
                }
            }

            // 内部領域の三角形要素
            IList<uint> triFEIds = new List<uint>();
            PeriodicPortTriangleFEIdss.Add(triFEIds);
            IList<uint> lIds = portCondition.LoopIdsForPeriodic;
            foreach (uint lId in lIds)
            {
                IList<uint> tmpTriFEIds = world.GetTriangleFEIdsFromLoopCadId(Id, lId);
                foreach (uint feId in tmpTriFEIds)
                {
                    if (triFEIds.IndexOf(feId) == -1)
                    {
                        triFEIds.Add(feId);

                        TriangleFE triFE = TriangleFEArray.GetObject(feId);
                        int[] coIds = triFE.NodeCoordIds;

                        foreach (int coId in coIds)
                        {
                            if (!portCo2Node.ContainsKey(coId) &&
                                zeroCoordIds.IndexOf(coId) == -1)
                            {
                                portCo2Node[coId] = portNodeId;
                                portNodeId++;
                            }
                        }
                    }
                }
            }
        }

        private void SortPortCoIds(
            FEWorld world, Mesher2D mesh, uint eId0, IList<int> bcCoIds, out IList<int> sortedCoIds)
        {
            sortedCoIds = null;

            // 境界の方向順に節点番号を振る
            OpenTK.Vector2d sPt;
            {
                Edge2D e = mesh.Cad2D.GetEdge(eId0);
                sPt = e.GetVertexCoord(true); // 始点の座標
            }
            var coIdDistances = new List<KeyValuePair<int, double>>();
            foreach (int coId in bcCoIds)
            {
                double[] coord = world.GetCoord(Id, coId);
                double distance = OpenTK.Vector2d.Distance(sPt, new OpenTK.Vector2d(coord[0], coord[1]));
                coIdDistances.Add(new KeyValuePair<int, double>(coId, distance));
            }
            coIdDistances.Sort((a, b) =>
            {
                // 距離を比較
                double diff = a.Value - b.Value;
                // 昇順
                if (diff > 0)
                {
                    return 1;
                }
                else if (diff < 0)
                {
                    return -1;
                }
                return 0;
            });

            sortedCoIds = new List<int>();
            {
                foreach (var coIdDistance in coIdDistances)
                {
                    int coId = coIdDistance.Key;
                    if (sortedCoIds.IndexOf(coId) == -1)
                    {
                        sortedCoIds.Add(coId);
                    }
                }
            }
        }

        private void SetDistributedPortCoords(FEWorld world)
        {
            uint portCnt = GetPortCount();
            for (int portId = 0; portId < portCnt; portId++)
            {
                PortCondition portCondition = PortConditions[portId];
                IList<uint> portEIds = portCondition.EIds;
                if (portCondition is DistributedPortCondition)
                {
                    DistributedPortCondition dist = portCondition as DistributedPortCondition;
                    IList<int> coIds = new List<int>(); 
                    foreach (uint eId in portEIds)
                    {
                        IList<int> tmpCoIds = GetCoordIdsFromCadId(world, eId, CadElementType.Edge);
                        foreach (int coId in tmpCoIds)
                        {
                            if (coIds.Contains(coId))
                            {
                                continue;
                            }
                            coIds.Add(coId);
                        }
                    }
                    dist.InitCoordIds(coIds);
                }
            }
        }

        // 接触解析のMaster/Slave線要素を準備する
        private void SetupContactMasterSlaveLineElements(FEWorld world)
        {
            Mesher2D mesh = world.Mesh;

            IList<uint> feIds = LineFEArray.GetObjectIds();
            foreach (var feId in feIds)
            {
                LineFE lineFE = LineFEArray.GetObject(feId);
                uint cadId;
                {
                    uint meshId = lineFE.MeshId;
                    uint elemCnt;
                    MeshType meshType;
                    int loc;
                    mesh.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadId);
                    System.Diagnostics.Debug.Assert(meshType == MeshType.Bar);
                }
                if (ContactSlaveEIds.Contains(cadId))
                {
                    // Slave上の線要素
                    ContactSlaveLineFEIds.Add(feId);
                }
                if (ContactMasterEIds.Contains(cadId))
                {
                    // Master上の線要素
                    ContactMasterLineFEIds.Add(feId);
                }
            }
        }

        // 三角形要素の節点ナンバリング
        private void NumberTriangleNodes(FEWorld world, IList<int> zeroCoordIds)
        {
            Mesher2D mesh = world.Mesh;

            // ナンバリング
            int nodeId = 0;
            IList<uint> feIds = TriangleFEArray.GetObjectIds();
            foreach (uint feId in feIds)
            {
                TriangleFE fe = TriangleFEArray.GetObject(feId);
                int elemPtCnt = fe.NodeCoordIds.Length;
                int[] coIds = fe.NodeCoordIds;
                for (int iPt = 0; iPt < elemPtCnt; iPt++)
                {
                    int coId = coIds[iPt];
                    if (!Co2Node.ContainsKey(coId) &&
                        zeroCoordIds.IndexOf(coId) == -1)
                    {
                        Co2Node[coId] = nodeId;
                        nodeId++;
                    }
                }

                uint meshId = fe.MeshId;
                int iElem = fe.MeshElemId;
                uint elemCnt;
                MeshType meshType;
                int loc;
                uint cadId;
                mesh.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadId);
                System.Diagnostics.Debug.Assert(meshType == MeshType.Tri);
                var triArray = mesh.GetTriArrays();
                var tri = triArray[loc].Tris[iElem];
                tri.FEId = (int)feId;

                string key = string.Format(meshId + "_" + iElem);
                Mesh2TriangleFE.Add(key, feId);
            }
        }

        private void MakeCo2AndEdgeCos2TriangleFE()
        {
            Co2TriangleFE.Clear();
            IList<uint> feIds = GetTriangleFEIds();
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = GetTriangleFE(feId);
                // 節点→要素
                {
                    int[] coIds = triFE.NodeCoordIds;
                    for (int i = 0; i < coIds.Length; i++)
                    {
                        int coId = coIds[i];
                        IList<uint> targetIds = null;
                        if (Co2TriangleFE.ContainsKey(coId))
                        {
                            targetIds = Co2TriangleFE[coId];
                        }
                        else
                        {
                            targetIds = new List<uint>();
                            Co2TriangleFE[coId] = targetIds;
                        }
                        targetIds.Add(feId);
                    }
                }
                // 辺(頂点1-頂点2)→要素
                {
                    int[] coIds = triFE.VertexCoordIds;
                    for (int i = 0; i < coIds.Length; i++)
                    {
                        int v1 = coIds[i];
                        int v2 = coIds[(i + 1) % coIds.Length];
                        if (v1 > v2)
                        {
                            int tmp = v1;
                            v1 = v2;
                            v2 = tmp;
                        }
                        string edgeKey = v1 + "_" + v2;
                        IList<uint> targetIds = null;
                        if (EdgeCos2TriangleFE.ContainsKey(edgeKey))
                        {
                            targetIds = EdgeCos2TriangleFE[edgeKey];
                        }
                        else
                        {
                            targetIds = new List<uint>();
                            EdgeCos2TriangleFE[edgeKey] = targetIds;
                        }
                        targetIds.Add(feId);
                    }
                }
            }
        }

        private void MakeCo2LineFE()
        {
            Co2LineFE.Clear();
            IList<uint> feIds = GetLineFEIds();
            foreach (uint feId in feIds)
            {
                LineFE lineFE = GetLineFE(feId);
                // 節点→要素
                {
                    int[] coIds = lineFE.NodeCoordIds;
                    for (int i = 0; i < coIds.Length; i++)
                    {
                        int coId = coIds[i];
                        IList<uint> targetIds = null;
                        if (Co2LineFE.ContainsKey(coId))
                        {
                            targetIds = Co2LineFE[coId];
                        }
                        else
                        {
                            targetIds = new List<uint>();
                            Co2LineFE[coId] = targetIds;
                        }
                        targetIds.Add(feId);
                    }
                }
            }
        }

        private IList<int> GetZeroCoordIds(FEWorld world)
        {
            IList<int> zeroCoIds = new List<int>();

            foreach (var fixedCad in ZeroFieldFixedCads)
            {
                IList<int> coIds = GetCoordIdsFromCadId(world, fixedCad.CadId, fixedCad.CadElemType);
                foreach (int coId in coIds)
                {
                    zeroCoIds.Add(coId);
                }
            }
            return zeroCoIds;
        }

        private void MakeCo2FixedCads(FEWorld world)
        {
            Co2FixedCads.Clear();
            foreach (var fixedCad in FieldFixedCads)
            {
                IList<int> coIds = GetCoordIdsFromCadId(world, fixedCad.CadId, fixedCad.CadElemType);
                foreach (int coId in coIds)
                {
                    IList<FieldFixedCad> fixedCads = null;
                    if (!Co2FixedCads.ContainsKey(coId))
                    {
                        fixedCads = new List<FieldFixedCad>();
                        Co2FixedCads[coId] = fixedCads;
                    }
                    else
                    {
                        fixedCads = Co2FixedCads[coId];
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

        private void SetDistributedFixedCadCoords(FEWorld world)
        {
            foreach (var fixedCad in FieldFixedCads)
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
            foreach (var periodicPortBcCo2Nodes in PeriodicPortBcCo2Nodess)
            {
                IList<Dictionary<int, int>> periodicPortBcNode2Cos = new List<Dictionary<int, int>>();
                PeriodicPortBcNode2Coss.Add(periodicPortBcNode2Cos);
                System.Diagnostics.Debug.Assert(periodicPortBcCo2Nodes.Count == 2);
                foreach (var periodicPortBcCo2Node in periodicPortBcCo2Nodes)
                {
                    var periodicPortBcNode2Co = new Dictionary<int, int>();
                    periodicPortBcNode2Cos.Add(periodicPortBcNode2Co);
                    foreach (var pair in periodicPortBcCo2Node)
                    {
                        int tmpPortCoId = pair.Key;
                        int tmpPortNodeId = pair.Value;
                        periodicPortBcNode2Co[tmpPortNodeId] = tmpPortCoId;
                    }
                }
            }
        }

        public IList<LineFE> MakeBoundOfElements(FEWorld world)
        {
            IList<LineFE> boundOfTriangelFEs = new List<LineFE>();
            HashSet<string> edges = new HashSet<string>();

            var feIds = GetTriangleFEIds();
            foreach (uint feId in feIds)
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
                    var lineFE = new LineFE((int)FEOrder, FEType);
                    lineFE.World = world;
                    lineFE.QuantityId = (int)this.Id;
                    lineFE.SetVertexCoordIds(vertexCoIds[iEdge]);
                    lineFE.SetNodeCoordIds(nodeCoIds[iEdge]);
                    // MeshId等は対応するものがないのでセットしない
                    boundOfTriangelFEs.Add(lineFE);
                }
            }
            return boundOfTriangelFEs;
        }
    }
}
