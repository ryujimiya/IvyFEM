using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class FEWorldQuantity
    {
        public void MakeElements3D(
            FEWorld world,
            IList<double> vertexCoords,
            Dictionary<uint, uint> cadSolid2Material,
            Dictionary<uint, uint> cadLoop2Material,
            Dictionary<uint, uint> cadEdge2Material)
        {
            ClearElements();

            // 座標、四面体要素、三角形要素、線要素を生成する
            MakeCoordsAndElements3D(
                world, vertexCoords, cadSolid2Material, cadLoop2Material, cadEdge2Material);

            // 多点拘束の座標生成
            MakeCo2MultipointConstraints(world);

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

            MakeCo2FixedCads(world, FieldFixedCads, Co2FixedCads);
            MakeCo2FixedCads(world, ForceFieldFixedCads, Co2ForceFixedCads);
            SetDistributedFixedCadCoords(world, FieldFixedCads);
            SetDistributedFixedCadCoords(world, ForceFieldFixedCads);

            // 頂点→四面体要素のマップ作成
            MakeCo2TetrahedronFE3D();

            // ポート上の三角形要素の節点ナンバリング
            NumberPortNodes3D(world, zeroCoordIds);
            SetDistributedPortCoords(world);

            // 四面体要素の節点ナンバリング
            NumberTetrahedronNodes3D(world, zeroCoordIds);

            // 接触解析のMaster/Slave三角形要素を準備する
            SetupContactMasterSlaveTriangleElements3D(world);

            // 節点→座標のマップ作成
            MakeNode2CoFromCo2Node();

            // 頂点→三角形要素のマップ
            MakeCo2TriangleFE();
            // 頂点→線要素のマップ
            MakeCo2LineFE();
        }

        // 座標、四面体要素、三角形要素と線要素を生成する
        private void MakeCoordsAndElements3D(
            FEWorld world,
            IList<double> vertexCoords,
            Dictionary<uint, uint> cadSolid2Material,
            Dictionary<uint, uint> cadLoop2Material,
            Dictionary<uint, uint> cadEdge2Material)
        {
            var mesh = world.Mesh;
            System.Diagnostics.Debug.Assert(mesh != null);

            if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 1)
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
            // 領域の四面体要素
            // まず要素を作る
            // この順番で生成した要素は隣接していない
            Dictionary<string, IList<int>> edge2MidPt = new Dictionary<string, IList<int>>();
            uint[][] tetEdgeNodes = new uint[6][] {
                new uint[2] { 0, 1 },
                new uint[2] { 1, 2 },
                new uint[2] { 2, 0 },
                new uint[2] { 0, 3 },
                new uint[2] { 1, 3 },
                new uint[2] { 2, 3 }
            };
            foreach (uint meshId in meshIds)
            {
                uint elemCnt;
                MeshType meshType;
                int loc;
                uint cadId;
                mesh.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadId);
                if (meshType != MeshType.Tet)
                {
                    continue;
                }

                if (!cadSolid2Material.ContainsKey(cadId))
                {
                    throw new IndexOutOfRangeException();
                }
                uint maId = cadSolid2Material[cadId];

                int elemVertexCnt = 4;
                int elemNodeCnt = 0;
                if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 1)
                {
                    elemNodeCnt = 4;
                }
                else if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 2)
                {
                    elemNodeCnt = 10;
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
                    int[] elemVertexCoIds = new int[elemVertexCnt];
                    for (int iPt = 0; iPt < elemVertexCnt; iPt++)
                    {
                        int coId = vertexs[iElem * elemVertexCnt + iPt];
                        elemVertexCoIds[iPt] = coId;
                    }
                    int[] elemNodeCoIds = new int[elemNodeCnt];
                    if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 1)
                    {
                        System.Diagnostics.Debug.Assert(elemNodeCoIds.Length == elemVertexCoIds.Length);
                        elemVertexCoIds.CopyTo(elemNodeCoIds, 0);
                    }
                    else if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 2)
                    {
                        for (int i = 0; i < elemVertexCnt; i++)
                        {
                            elemNodeCoIds[i] = elemVertexCoIds[i];
                        }

                        for (int iEdge = 0; iEdge < 6; iEdge++)
                        {
                            uint[] edgeNode = tetEdgeNodes[iEdge];
                            int v1 = elemVertexCoIds[edgeNode[0]];
                            int v2 = elemVertexCoIds[edgeNode[1]];
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
                                uint dim = Dimension;
                                System.Diagnostics.Debug.Assert(vPt1.Length == dim);
                                System.Diagnostics.Debug.Assert(vPt2.Length == dim);
                                double[] midPt = new double[dim];
                                for (int idim = 0; idim < dim; idim++)
                                {
                                    midPt[idim] = (vPt1[idim] + vPt2[idim]) / 2.0;
                                }
                                midPtCoId = (int)(Coords.Count / dim);
                                for (int idim = 0; idim < dim; idim++)
                                {
                                    Coords.Add(midPt[idim]);
                                }
                                var list = new List<int>();
                                list.Add(midPtCoId);
                                edge2MidPt[edgeKey] = list;
                            }

                            elemNodeCoIds[iEdge + elemVertexCnt] = midPtCoId;
                        }
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }

                    TetrahedronFE fe = new TetrahedronFE((int)FEOrder, FEType);
                    fe.World = world;
                    fe.QuantityId = (int)this.Id;
                    fe.QuantityIdBaseOffset = IdBaseOffset;
                    fe.SetVertexCoordIds(elemVertexCoIds);
                    fe.SetNodeCoordIds(elemNodeCoIds);
                    fe.MaterialId = maId;
                    fe.MeshId = meshId;
                    fe.MeshElemId = iElem;
                    // 仮登録
                    uint freeId = TetrahedronFEArray.GetFreeObjectId();
                    uint feId = TetrahedronFEArray.AddObject(freeId, fe);
                    System.Diagnostics.Debug.Assert(feId == freeId);
                }
            }

            //////////////////////////////////////////////////
            // 境界の三角形要素
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

                int elemVertexCnt = 3;
                int elemNodeCnt = 0;
                if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 1)
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

                // 未指定のマテリアルも許容する
                uint maId = cadLoop2Material.ContainsKey(cadId) ? cadLoop2Material[cadId] : 0;

                for (int iElem = 0; iElem < elemCnt; iElem++)
                {
                    int[] elemVertexCoIds = new int[elemVertexCnt];
                    for (int iPt = 0; iPt < elemVertexCnt; iPt++)
                    {
                        int coId = vertexs[iElem * elemVertexCnt + iPt];
                        elemVertexCoIds[iPt] = coId;
                    }
                    uint triFEOrder = FEOrder;
                    int[] elemNodeCoIds = new int[elemNodeCnt];
                    if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 1)
                    {
                        System.Diagnostics.Debug.Assert(elemNodeCoIds.Length == elemVertexCoIds.Length);
                        elemVertexCoIds.CopyTo(elemNodeCoIds, 0);
                    }
                    else if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 2)
                    {
                        for (int i = 0; i < elemVertexCnt; i++)
                        {
                            elemNodeCoIds[i] = elemVertexCoIds[i];

                            {
                                int v1 = elemVertexCoIds[i];
                                int v2 = elemVertexCoIds[(i + 1) % elemVertexCnt];
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
                                    uint dim = Dimension;
                                    System.Diagnostics.Debug.Assert(vPt1.Length == dim);
                                    System.Diagnostics.Debug.Assert(vPt2.Length == dim);
                                    double[] midPt = new double[dim];
                                    for (int idim = 0; idim < dim; idim++)
                                    {
                                        midPt[idim] = (vPt1[idim] + vPt2[idim]) / 2.0;
                                    }
                                    midPtCoId = (int)(Coords.Count / dim);
                                    for (int idim = 0; idim < dim; idim++)
                                    {
                                        Coords.Add(midPt[idim]);
                                    }
                                    var list = new List<int>();
                                    list.Add(midPtCoId);
                                    edge2MidPt[edgeKey] = list;
                                }

                                elemNodeCoIds[i + elemVertexCnt] = midPtCoId;
                            }
                        }
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }

                    TriangleFE triFE = new TriangleFE((int)triFEOrder, FEType);
                    triFE.World = world;
                    triFE.QuantityId = (int)this.Id;
                    triFE.SetVertexCoordIds(elemVertexCoIds);
                    triFE.SetNodeCoordIds(elemNodeCoIds);
                    triFE.MaterialId = maId;
                    triFE.MeshId = meshId;
                    triFE.MeshElemId = iElem;
                    uint freeId = TriangleFEArray.GetFreeObjectId();
                    uint feId = TriangleFEArray.AddObject(freeId, triFE);
                    System.Diagnostics.Debug.Assert(feId == freeId);

                    string key = string.Format(meshId + "_" + iElem);
                    Mesh2TriangleFE.Add(key, feId);
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

                // 未指定のマテリアルも許容する
                uint maId = cadEdge2Material.ContainsKey(cadId) ? cadEdge2Material[cadId] : 0;

                for (int iElem = 0; iElem < elemCnt; iElem++)
                {
                    int[] elemVertexCoIds = new int[elemVertexCnt];
                    for (int iPt = 0; iPt < elemVertexCnt; iPt++)
                    {
                        int coId = vertexs[iElem * elemVertexCnt + iPt];
                        elemVertexCoIds[iPt] = coId;
                    }
                    uint lineFEOrder = FEOrder;
                    int[] elemNodeCoIds = new int[elemNodeCnt];
                    if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 1)
                    {
                        System.Diagnostics.Debug.Assert(elemNodeCoIds.Length == elemVertexCoIds.Length);
                        elemVertexCoIds.CopyTo(elemNodeCoIds, 0);
                    }
                    else if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 2)
                    {
                        for (int i = 0; i < 2; i++)
                        {
                            elemNodeCoIds[i] = elemVertexCoIds[i];
                        }
                        // 線要素上の中点
                        int v1 = elemVertexCoIds[0];
                        int v2 = elemVertexCoIds[1];
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
                            uint dim = Dimension;
                            System.Diagnostics.Debug.Assert(vPt1.Length == dim);
                            System.Diagnostics.Debug.Assert(vPt2.Length == dim);
                            double[] midPt = new double[dim];
                            for (int idim = 0; idim < dim; idim++)
                            {
                                midPt[idim] = (vPt1[idim] + vPt2[idim]) / 2.0;
                            }
                            midPtCoId = (int)(Coords.Count / dim);
                            for (int idim = 0; idim < dim; idim++)
                            {
                                Coords.Add(midPt[idim]);
                            }
                            var list = new List<int>();
                            list.Add(midPtCoId);
                            edge2MidPt[edgeKey] = list;
                        }
                        elemNodeCoIds[2] = midPtCoId;
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }

                    LineFE lineFE = new LineFE((int)lineFEOrder, FEType);
                    lineFE.World = world;
                    lineFE.QuantityId = (int)this.Id;
                    lineFE.SetVertexCoordIds(elemVertexCoIds);
                    lineFE.SetNodeCoordIds(elemNodeCoIds);
                    lineFE.MaterialId = maId;
                    lineFE.MeshId = meshId;
                    lineFE.MeshElemId = iElem;
                    uint freeId = LineFEArray.GetFreeObjectId();
                    uint feId = LineFEArray.AddObject(freeId, lineFE);
                    System.Diagnostics.Debug.Assert(feId == freeId);

                    string key = string.Format(meshId + "_" + iElem);
                    Mesh2LineFE.Add(key, feId);
                }
            }
        }

        // ポート上の三角形要素の節点ナンバリング
        private void NumberPortNodes3D(FEWorld world, IList<int> zeroCoordIds)
        {
            var mesh = world.Mesh;

            // ポート上の三角形要素の抽出と節点ナンバリング
            uint portCnt = GetPortCount();
            PortTriangleFEIdss.Clear();
            PortCo2Nodes.Clear();
            for (int portId = 0; portId < portCnt; portId++)
            {
                PortTriangleFEIdss.Add(new List<uint>());
                PortCo2Nodes.Add(new Dictionary<int, int>());
            }

            for (int portId = 0; portId < portCnt; portId++)
            {
                PortCondition portCondition = PortConditions[portId];
                NumberNormalPortNodes3D((uint)portId, world, mesh, portCondition, zeroCoordIds);
            }
        }

        // 通常のポート（境界のみ）
        private void NumberNormalPortNodes3D(
            uint portId, FEWorld world, IMesher mesh, PortCondition portCondition, IList<int> zeroCoordIds)
        {
            System.Diagnostics.Debug.Assert(!portCondition.IsPeriodic);
            IList<uint> portLIds = portCondition.LIds;
            var triFEIds = new List<uint>();
            PortTriangleFEIdss[(int)portId] = triFEIds;

            var portCo2Node = new Dictionary<int, int>();
            PortCo2Nodes[(int)portId] = portCo2Node;
            int portNodeId = 0;

            IList<uint> feIds = TriangleFEArray.GetObjectIds();
            IList<int> portCoIds = new List<int>();
            foreach (var feId in feIds)
            {
                TriangleFE triFE = TriangleFEArray.GetObject(feId);
                uint cadId;
                {
                    uint meshId = triFE.MeshId;
                    uint elemCnt;
                    MeshType meshType;
                    int loc;
                    mesh.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadId);
                    System.Diagnostics.Debug.Assert(meshType == MeshType.Tri);
                }
                if (portLIds.Contains(cadId))
                {
                    // ポート上の三角形要素
                    triFEIds.Add(feId);

                    int[] coIds = triFE.NodeCoordIds;

                    foreach (int coId in coIds)
                    {
                        if (portCoIds.IndexOf(coId) == -1)
                        {
                            portCoIds.Add(coId);
                        }
                    }
                }
            }
            // 起点からの距離で節点番号を振る
            uint lId0 = portLIds[0];
            IList<int> sortedCoIds;
            SortPortCoIds3D(world, mesh, lId0, portCoIds, out sortedCoIds);
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

        private void SortPortCoIds3D(
            FEWorld world, IMesher mesh, uint lId0, IList<int> bcCoIds, out IList<int> sortedCoIds)
        {
            sortedCoIds = null;

            uint dim = Dimension;
            // 始点からの距離で番号を振る
            double[] pt0;
            {
                uint meshId = mesh.GetIdFromCadId(lId0, CadElementType.Loop);
                MeshType meshType;
                uint elemCnt;
                int[] vertexs;
                mesh.GetConnectivity(meshId, out meshType, out vertexs);
                int coId = vertexs[0]; // 始点の座標
                pt0 = GetCoord(coId, world.RotAngle, world.RotOrigin);
            }

            var coIdDists = new List<KeyValuePair<int, double>>();
            foreach (int coId in bcCoIds)
            {
                double[] pt = world.GetCoord(Id, coId);
                double dist = CadUtils.GetDistance(pt0, pt); 
                coIdDists.Add(new KeyValuePair<int, double>(coId, dist));
            }
            coIdDists.Sort((a, b) =>
            {
                // 座標を比較
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
                foreach (var coIdDist in coIdDists)
                {
                    int coId = coIdDist.Key;
                    if (sortedCoIds.IndexOf(coId) == -1)
                    {
                        sortedCoIds.Add(coId);
                    }
                }
            }
        }

        // 四面体要素の節点ナンバリング
        private void NumberTetrahedronNodes3D(FEWorld world, IList<int> zeroCoordIds)
        {
            if (FEType == FiniteElementType.Edge)
            {
                return;
            }
            var mesh = world.Mesh;

            // ナンバリング
            int nodeId = 0;
            IList<uint> feIds = TetrahedronFEArray.GetObjectIds();
            foreach (uint feId in feIds)
            {
                TetrahedronFE fe = TetrahedronFEArray.GetObject(feId);
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
                System.Diagnostics.Debug.Assert(meshType == MeshType.Tet);

                string key = string.Format(meshId + "_" + iElem);
                Mesh2TetrahedronFE.Add(key, feId);
            }
        }

        private void MakeCo2TetrahedronFE3D()
        {
            Co2TetrahedronFE.Clear();
            IList<uint> feIds = GetTetrahedronFEIds();
            foreach (uint feId in feIds)
            {
                TetrahedronFE tetFE = GetTetrahedronFE(feId);
                // 節点→要素
                {
                    int[] coIds = tetFE.NodeCoordIds;
                    for (int i = 0; i < coIds.Length; i++)
                    {
                        int coId = coIds[i];
                        IList<uint> targetIds = null;
                        if (Co2TetrahedronFE.ContainsKey(coId))
                        {
                            targetIds = Co2TetrahedronFE[coId];
                        }
                        else
                        {
                            targetIds = new List<uint>();
                            Co2TetrahedronFE[coId] = targetIds;
                        }
                        targetIds.Add(feId);
                    }
                }
            }
        }

        // 接触解析のMaster/Slave三角形要素を準備する
        private void SetupContactMasterSlaveTriangleElements3D(FEWorld world)
        {
            var mesh = world.Mesh;

            IList<uint> feIds = TriangleFEArray.GetObjectIds();
            foreach (var feId in feIds)
            {
                TriangleFE triFE = TriangleFEArray.GetObject(feId);
                uint cadId;
                {
                    uint meshId = triFE.MeshId;
                    uint elemCnt;
                    MeshType meshType;
                    int loc;
                    mesh.GetMeshInfo(meshId, out elemCnt, out meshType, out loc, out cadId);
                    System.Diagnostics.Debug.Assert(meshType == MeshType.Tri);
                }
                if (ContactSlaveCadIds.Contains(cadId))
                {
                    // Slave上の三角形要素
                    ContactSlaveFEIds.Add(feId);
                }
                if (ContactMasterCadIds.Contains(cadId))
                {
                    // Master上の三角形要素
                    ContactMasterFEIds.Add(feId);
                }
            }
        }
    }
}
