using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class FEWorldQuantity
    {
        public void MakeElements2D(
            FEWorld world,
            IList<double> vertexCoords,
            Dictionary<uint, uint> cadLoop2Material,
            Dictionary<uint, uint> cadEdge2Material)
        {
            ClearElements();

            // 座標、三角形要素と線要素を生成する
            MakeCoordsAndElements2D(
                world, vertexCoords, cadLoop2Material, cadEdge2Material);

            // 多点拘束の座標生成
            MakeCo2MultipointConstraints(world);

            // Note: 線要素生成後でないと点を特定できない
            IList<int> zeroCoordIds = GetZeroCoordIds(world);
            IList<int> zeroEdgeIds = GetZeroEdgeIds(world);

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
            MakeCo2FixedCads(world, FieldFixedCads, Co2FixedCads);
            MakeCo2FixedCads(world, ForceFieldFixedCads, Co2ForceFixedCads);
            SetDistributedFixedCadCoords(world, FieldFixedCads);
            SetDistributedFixedCadCoords(world, ForceFieldFixedCads);

            // 頂点→三角形要素のマップと辺→三角形要素のマップ作成
            MakeCo2TriangleFE();
            MakeEdgeCos2TriangleFE2D();

            // ポート上の線要素の節点ナンバリング
            NumberPortNodes2D(world, zeroCoordIds);
            SetDistributedPortCoords(world);

            if (IsPortOnly)
            {
                NumberNodesPortOnly2D(world, zeroCoordIds);
            }
            else
            {
                // 三角形要素の節点ナンバリング
                NumberTriangleNodes2D(world, zeroCoordIds);
                NumberTriangleEdgeNodes2D(world, zeroEdgeIds);
                // Loopに属さない線要素の節点ナンバリング
                NumberLineNodesNotBelongToLoop2D(world, zeroCoordIds);
            }

            // 接触解析のMaster/Slave線要素を準備する
            SetupContactMasterSlaveLineElements2D(world);

            // 節点→座標のマップ作成
            MakeNode2CoFromCo2Node();
            MakeENode2EdgeFromEdge2ENode();

            // 頂点→線要素のマップ
            MakeCo2LineFE();
        }

        // 座標、三角形要素と線要素を生成する
        private void MakeCoordsAndElements2D(
            FEWorld world,
            IList<double> vertexCoords,
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
            else if (FEType == FiniteElementType.ScalarConstant && FEOrder == 0)
            {
                Coords = new List<double>(vertexCoords);
            }
            else if (FEType == FiniteElementType.ScalarBell && FEOrder == 5)
            {
                Coords = new List<double>(vertexCoords);
            }
            else if (FEType == FiniteElementType.ScalarHermite && FEOrder == 3)
            {
                Coords = new List<double>(vertexCoords);
            }
            else if (FEType == FiniteElementType.Edge && FEOrder == 1)
            {
                Coords = new List<double>(vertexCoords);
            }
            else if (FEType == FiniteElementType.Edge && FEOrder == 2)
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
                else if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 2)
                {
                    elemNodeCnt = 6;
                }
                else if (FEType == FiniteElementType.ScalarConstant && FEOrder == 0)
                {
                    elemNodeCnt = 1;
                }
                else if (FEType == FiniteElementType.ScalarBell && FEOrder == 5)
                {
                    elemNodeCnt = 3;
                }
                else if (FEType == FiniteElementType.ScalarHermite && FEOrder == 3)
                {
                    elemNodeCnt = 3;
                }
                else if (FEType == FiniteElementType.Edge && FEOrder == 1)
                {
                    elemNodeCnt = 3;
                }
                else if (FEType == FiniteElementType.Edge && FEOrder == 2)
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
                    else if (FEType == FiniteElementType.ScalarConstant && FEOrder == 0)
                    {
                        // 重心座標を求める
                        uint dim = Dimension;
                        double[] gCoord = new double[dim];
                        System.Diagnostics.Debug.Assert(elemVertexCoIds.Length == 3); // 三角形
                        for (int iVertex = 0; iVertex < elemVertexCoIds.Length; iVertex++)
                        {
                            int vCoId = elemVertexCoIds[iVertex];
                            double[] vCoord = GetCoord(vCoId, world.RotAngle, world.RotOrigin);
                            for (int iDim = 0; iDim < dim; iDim++)
                            {
                                gCoord[iDim] += vCoord[iDim];
                            }
                        }
                        for (int iDim = 0; iDim < dim; iDim++)
                        {
                            gCoord[iDim] /= 3.0;
                        }
                        int gCoId = (int)(Coords.Count / dim);
                        for (int iDim = 0; iDim < dim; iDim++)
                        {
                            Coords.Add(gCoord[iDim]);
                        }
                        elemNodeCoIds[0] = gCoId;
                    }
                    else if (FEType == FiniteElementType.ScalarBell && FEOrder == 5)
                    {
                        System.Diagnostics.Debug.Assert(elemNodeCoIds.Length == elemVertexCoIds.Length);
                        elemVertexCoIds.CopyTo(elemNodeCoIds, 0);
                    }
                    else if (FEType == FiniteElementType.ScalarHermite && FEOrder == 3)
                    {
                        // 暫定: Lagrange三角形要素で代用
                        System.Diagnostics.Debug.Assert(elemNodeCoIds.Length == elemVertexCoIds.Length);
                        elemVertexCoIds.CopyTo(elemNodeCoIds, 0);
                    }
                    else if (FEType == FiniteElementType.Edge && FEOrder == 1)
                    {
                        System.Diagnostics.Debug.Assert(elemNodeCoIds.Length == elemVertexCoIds.Length);
                        elemVertexCoIds.CopyTo(elemNodeCoIds, 0);
                    }
                    else if (FEType == FiniteElementType.Edge && FEOrder == 2)
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

                    TriangleFE fe = new TriangleFE((int)FEOrder, FEType);
                    fe.World = world;
                    fe.QuantityId = (int)this.Id;
                    fe.QuantityIdBaseOffset = IdBaseOffset;
                    fe.SetVertexCoordIds(elemVertexCoIds);
                    fe.SetNodeCoordIds(elemNodeCoIds);
                    if (FEType == FiniteElementType.Edge)
                    {
                        fe.SetEdgeCoordIdsFromNodeCoordIds();

                        // 辺の登録
                        int[][] edgeCoIdss = fe.EdgeCoordIdss;
                        foreach (int[] edgeCoIds in edgeCoIdss)
                        {
                            int coId1 = edgeCoIds[0];
                            int coId2 = edgeCoIds[1];
                            string edgeKey = coId1 < coId2 ?
                                string.Format("{0}_{1}", coId1, coId2) :
                                string.Format("{0}_{1}", coId2, coId1);
                            if (!Edges.Contains(edgeKey))
                            {
                                Edges.Add(edgeKey);
                            }
                        }
                    }
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
                else if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 2)
                {
                    elemNodeCnt = 3;
                }
                else if (FEType == FiniteElementType.ScalarConstant && FEOrder == 0)
                {
                    // 暫定：Lagrange線要素で代用
                    elemNodeCnt = 2;
                }
                else if (FEType == FiniteElementType.ScalarBell && FEOrder == 5)
                {
                    // 暫定：Lagrange線要素で代用
                    elemNodeCnt = 2;
                }
                else if (FEType == FiniteElementType.ScalarHermite && FEOrder == 3)
                {
                    // 暫定：Lagrange線要素で代用
                    elemNodeCnt = 2;
                }
                else if (FEType == FiniteElementType.Edge && FEOrder == 1)
                {
                    // 暫定: Lagrange線要素で代用
                    elemNodeCnt = 2;
                }
                else if (FEType == FiniteElementType.Edge && FEOrder == 2)
                {
                    // 暫定: Lagrange線要素で代用
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
                    else if (FEType == FiniteElementType.ScalarConstant && FEOrder == 0)
                    {
                        // 暫定：Lagrange線要素で代用
                        System.Diagnostics.Debug.Assert(elemNodeCoIds.Length == elemVertexCoIds.Length);
                        elemVertexCoIds.CopyTo(elemNodeCoIds, 0);
                    }
                    else if (FEType == FiniteElementType.ScalarBell && FEOrder == 5)
                    {
                        // 暫定：Lagrange線要素で代用
                        System.Diagnostics.Debug.Assert(elemNodeCoIds.Length == elemVertexCoIds.Length);
                        elemVertexCoIds.CopyTo(elemNodeCoIds, 0);
                    }
                    else if (FEType == FiniteElementType.ScalarHermite && FEOrder == 3)
                    {
                        // 暫定：Lagrange線要素で代用
                        System.Diagnostics.Debug.Assert(elemNodeCoIds.Length == elemVertexCoIds.Length);
                        elemVertexCoIds.CopyTo(elemNodeCoIds, 0);
                    }
                    else if (FEType == FiniteElementType.Edge && FEOrder == 1)
                    {
                        // 暫定：Lagrange線要素で代用
                        System.Diagnostics.Debug.Assert(elemNodeCoIds.Length == elemVertexCoIds.Length);
                        elemVertexCoIds.CopyTo(elemNodeCoIds, 0);
                    }
                    else if (FEType == FiniteElementType.Edge && FEOrder == 2)
                    {
                        // 暫定：Lagrange線要素で代用
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
                    if (FEType == FiniteElementType.Edge)
                    {
                        lineFE.SetEdgeCoordIdsFromNodeCoordIds();
                    }
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

        // ポート上の線要素の節点ナンバリング
        private void NumberPortNodes2D(FEWorld world, IList<int> zeroCoordIds)
        {
            var mesh = world.Mesh;

            // ポート上の線要素の抽出と節点ナンバリング
            uint portCnt = GetPortCount();
            for (int portId = 0; portId < portCnt; portId++)
            {
                PortLineFEIdss.Add(new List<uint>());
                PortCo2Nodes.Add(new Dictionary<int, int>());
                PeriodicPortLineFEIdsss.Add(new List<IList<uint>>());
                PeriodicPortBcCosss.Add(new List<IList<int>>());
                PeriodicPortTriangleFEIdss.Add(new List<uint>());
            }

            for (int portId = 0; portId < portCnt; portId++)
            {
                PortCondition portCondition = PortConditions[portId];
                if (portCondition.IsPeriodic)
                {
                    NumberPeriodicPortNodes2D((uint)portId, world, mesh, portCondition, zeroCoordIds);
                }
                else
                {
                    NumberNormalPortNodes2D((uint)portId, world, mesh, portCondition, zeroCoordIds);
                }
            }
        }

        // 通常のポート（境界のみ）
        private void NumberNormalPortNodes2D(
            uint portId, FEWorld world, IMesher mesh, PortCondition portCondition, IList<int> zeroCoordIds)
        {
            System.Diagnostics.Debug.Assert(!portCondition.IsPeriodic);
            IList<uint> portEIds = portCondition.EIds;
            var lineFEIds = new List<uint>();
            PortLineFEIdss[(int)portId] = lineFEIds;

            var portCo2Node = new Dictionary<int, int>();
            PortCo2Nodes[(int)portId] = portCo2Node;
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
            uint eId1 = portEIds[0];
            uint eId2 = portEIds[portEIds.Count - 1];
            IList<int> sortedCoIds;
            SortPortCoIds2D(world, mesh, eId1, eId2, portCoIds, out sortedCoIds);
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

        // 周期構造のポート（2つ(or4つ)の境界と内部領域）
        private void NumberPeriodicPortNodes2D(
            uint portId, FEWorld world, IMesher mesh, PortCondition portCondition, IList<int> zeroCoordIds)
        {
            System.Diagnostics.Debug.Assert(portCondition.IsPeriodic);

            var portCo2Node = new Dictionary<int, int>();
            PortCo2Nodes[(int)portId] = portCo2Node;
            int portNodeId = 0;

            // 境界1、2(、3、4)、内部の順に番号を振る
            // 2つの境界
            IList<IList<uint>> portLineFEIdss = new List<IList<uint>>();
            PeriodicPortLineFEIdsss[(int)portId] = portLineFEIdss;
            IList<IList<int>> bcCoss = new List<IList<int>>();
            PeriodicPortBcCosss[(int)portId] = bcCoss;
            IList<uint>[] portEIdss = {
                portCondition.BcEIdsForPeriodic1, portCondition.BcEIdsForPeriodic2,
                portCondition.BcEIdsForPeriodic3, portCondition.BcEIdsForPeriodic4
            };
            for (int bcIndex = 0; bcIndex < portEIdss.Length; bcIndex++)
            {
                IList<uint> portEIds = portEIdss[bcIndex];
                IList<int> bcCos = new List<int>();
                bcCoss.Add(bcCos);
                var lineFEIds = new List<uint>();
                portLineFEIdss.Add(lineFEIds);

                if (portEIds == null)
                {
                    // 上下の場合に限る
                    System.Diagnostics.Debug.Assert(bcIndex == 2 || bcIndex == 3);
                    continue;
                }

                IList<int> workBcCoIds = new List<int>();
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
                            if (workBcCoIds.IndexOf(coId) == -1)
                            {
                                workBcCoIds.Add(coId);
                            }
                        }
                    }
                }
                // 境界の方向順に節点番号を振る
                uint eId1 = portEIds[0];
                uint eId2 = portEIds[portEIds.Count - 1];
                IList<int> sortedCoIds;
                SortPortCoIds2D(world, mesh, eId1, eId2, workBcCoIds, out sortedCoIds);
                foreach (int coId in sortedCoIds)
                {
                    if (zeroCoordIds.IndexOf(coId) == -1)
                    {
                        if (!portCo2Node.ContainsKey(coId))
                        {
                            portCo2Node[coId] = portNodeId;
                            bcCos.Add(coId);
                            portNodeId++;
                        }
                        else
                        {
                            // すでにナンバリングされている境界と境界のコーナーの点
                            //int sharedNodeId = portCo2Node[coId];
                            bcCos.Add(coId);
                        }
                    }
                }
            }

            // 内部領域の三角形要素
            IList<uint> triFEIds = new List<uint>();
            PeriodicPortTriangleFEIdss[(int)portId] = triFEIds;
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

        private void SortPortCoIds2D(
            FEWorld world, IMesher mesh, uint eId1, uint eId2, IList<int> bcCoIds, out IList<int> sortedCoIds)
        {
            sortedCoIds = null;

            uint dim = Dimension;
            // 境界の方向順に節点番号を振る
            double[] pt1;
            {
                uint meshId = mesh.GetIdFromCadId(eId1, CadElementType.Edge);
                MeshType meshType;
                uint elemCnt;
                int[] vertexs;
                mesh.GetConnectivity(meshId, out meshType, out vertexs);
                int coId = vertexs[0]; // 始点の座標
                pt1 = GetCoord(coId, world.RotAngle, world.RotOrigin);
            }
            double[] pt2;
            {
                uint meshId = mesh.GetIdFromCadId(eId2, CadElementType.Edge);
                MeshType meshType;
                uint elemCnt;
                int[] vertexs;
                mesh.GetConnectivity(meshId, out meshType, out vertexs);
                int coId = vertexs[1]; // 終点の座標
                pt2 = GetCoord(coId, world.RotAngle, world.RotOrigin);
            }
            double[] dir = CadUtils.GetDirection(pt1, pt2);

            var coIdLineXs = new List<KeyValuePair<int, double>>();
            foreach (int coId in bcCoIds)
            {
                double[] pt = world.GetCoord(Id, coId);
                // pt - pt1
                double[] vec = new double[Dimension];
                for (int idim = 0; idim < dim; idim++)
                {
                    vec[idim] = pt[idim] - pt1[idim];
                }
                double lineX = CadUtils.Dot(dir, vec);
                coIdLineXs.Add(new KeyValuePair<int, double>(coId, lineX));
            }
            coIdLineXs.Sort((a, b) =>
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
                foreach (var coIdLineX in coIdLineXs)
                {
                    int coId = coIdLineX.Key;
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
                CadElementType cadElemType = portCondition.CadElemType;
                IList<uint> portCadIds = portCondition.CadIds;
                if (portCondition is DistributedPortCondition)
                {
                    DistributedPortCondition dist = portCondition as DistributedPortCondition;
                    IList<int> coIds = new List<int>();
                    foreach (uint cadId in portCadIds)
                    {
                        IList<int> tmpCoIds = GetCoordIdsFromCadId(world, cadId, cadElemType);
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

        // 三角形要素の節点ナンバリング
        private void NumberTriangleNodes2D(FEWorld world, IList<int> zeroCoordIds)
        {
            if (IsPortOnly)
            {
                return;
            }
            if (FEType == FiniteElementType.Edge)
            {
                return;
            }
            var mesh = world.Mesh;

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

                string key = string.Format(meshId + "_" + iElem);
                Mesh2TriangleFE.Add(key, feId);
            }
        }

        // 三角形要素の節点ナンバリング(辺節点）
        private void NumberTriangleEdgeNodes2D(FEWorld world, IList<int> zeroEdgeIds)
        {
            if (IsPortOnly)
            {
                return;
            }
            if (FEType != FiniteElementType.Edge)
            {
                return;
            }
            var mesh = world.Mesh;

            // ナンバリング
            int edgeNodeId = 0;
            IList<uint> feIds = TriangleFEArray.GetObjectIds();
            foreach (uint feId in feIds)
            {
                TriangleFE fe = TriangleFEArray.GetObject(feId);
                int elemEdgeCnt = fe.EdgeCoordIdss.Length;
                int[][] edgeCoIdss = fe.EdgeCoordIdss;
                for (int iEdge = 0; iEdge < elemEdgeCnt; iEdge++)
                {
                    int coId1 = edgeCoIdss[iEdge][0];
                    int coId2 = edgeCoIdss[iEdge][1];
                    bool isReverse;
                    int edgeId = GetEdgeIdFromCoords(coId1, coId2, out isReverse);
                    System.Diagnostics.Debug.Assert(edgeId != -1);

                    // FIXME: 
                    if (!Edge2ENode.ContainsKey(edgeId) &&
                        zeroEdgeIds.IndexOf(edgeId) == -1)
                    {
                        Edge2ENode[edgeId] = edgeNodeId;
                        edgeNodeId++;
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

                string key = string.Format(meshId + "_" + iElem);
                Mesh2TriangleFE.Add(key, feId);
            }
        }

        // ポートのみの場合の節点番号割り振り
        private void NumberNodesPortOnly2D(FEWorld world, IList<int> zeroCoordIds)
        {
            if (!IsPortOnly)
            {
                return;
            }
            var mesh = world.Mesh;

            // ポート上の線要素の抽出と節点ナンバリング
            uint portCnt = GetPortCount();

            // ナンバリング
            int nodeId = 0;
            for (int portId = 0; portId < portCnt; portId++)
            {
                PortCondition portCondition = PortConditions[portId];
                _NumberNodesPortOnly2D((uint)portId, world, mesh, portCondition, zeroCoordIds, ref nodeId);
            }
        }

        private void _NumberNodesPortOnly2D(
            uint portId, FEWorld world, IMesher mesh, PortCondition portCondition, IList<int> zeroCoordIds, ref int nodeId)
        {
            if (!IsPortOnly)
            {
                return;
            }
            System.Diagnostics.Debug.Assert(!portCondition.IsPeriodic);
            IList<uint> portEIds = portCondition.EIds;

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
            uint eId1 = portEIds[0];
            uint eId2 = portEIds[portEIds.Count - 1];
            IList<int> sortedCoIds;
            SortPortCoIds2D(world, mesh, eId1, eId2, portCoIds, out sortedCoIds);
            foreach (int coId in sortedCoIds)
            {
                if (!Co2Node.ContainsKey(coId) &&
                    zeroCoordIds.IndexOf(coId) == -1)
                {
                    Co2Node[coId] = nodeId;
                    nodeId++;
                }
            }
        }

        // Loopに属さない線要素の節点ナンバリング
        private void NumberLineNodesNotBelongToLoop2D(FEWorld world, IList<int> zeroCoordIds)
        {
            if (FEType == FiniteElementType.Edge)
            {
                return;
            }
            if (FEType == FiniteElementType.ScalarConstant)
            {
                return;
            }
            var mesh = world.Mesh;

            // ナンバリング
            // これまでに三角形要素の節点番号としてナンバリングした数を取得する
            int nodeId = Co2Node.Count;
            IList<uint> lineFEIds = LineFEArray.GetObjectIds();
            foreach (uint feId in lineFEIds)
            {
                LineFE fe = LineFEArray.GetObject(feId);
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
                System.Diagnostics.Debug.Assert(meshType == MeshType.Bar);

                string key = string.Format(meshId + "_" + iElem);
                if (!Mesh2LineFE.ContainsKey(key))
                {
                    Mesh2LineFE.Add(key, feId);
                }
            }
        }

        private void MakeCo2TriangleFE()
        {
            System.Diagnostics.Debug.Assert(Co2TriangleFE.Count == 0);
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
            }
        }

        private void MakeEdgeCos2TriangleFE2D()
        {
            System.Diagnostics.Debug.Assert(EdgeCos2TriangleFE.Count == 0);
            IList<uint> feIds = GetTriangleFEIds();
            foreach (uint feId in feIds)
            {
                TriangleFE triFE = GetTriangleFE(feId);
                if (FEType == FiniteElementType.Edge)
                {
                    // 辺(頂点2-3) [1次edge element]→要素
                    // 辺(頂点2-5, 5-3) [2次edge element]→要素
                    System.Diagnostics.Debug.Assert(FEOrder == 1 || FEOrder == 2);
                    int elemEdgeCnt = (int)triFE.EdgeCount;
                    for (int eIndex = 0; eIndex < elemEdgeCnt; eIndex++)
                    {
                        int[] coIds = triFE.EdgeCoordIdss[eIndex];
                        int v1 = coIds[0];
                        int v2 = coIds[1];
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
                else
                {
                    // Scalar elements
                    bool isUseVertex = false;
                    int elemNodeCnt = (int)triFE.NodeCount;
                    int[][] edgeNos = null;
                    if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 1)
                    {
                        edgeNos = new int[3][]
                        {
                            new int[] { 0, 1 },
                            new int[] { 1, 2 },
                            new int[] { 2, 0 }
                        };
                    }
                    else if (FEType == FiniteElementType.ScalarLagrange && FEOrder == 2)
                    {
                        edgeNos = new int[6][]
                        {
                            new int[] { 0, 3 },
                            new int[] { 3, 1 },
                            new int[] { 1, 4 },
                            new int[] { 4, 2 },
                            new int[] { 2, 5 },
                            new int[] { 5, 0 }
                        };
                    }
                    else if (FEType == FiniteElementType.ScalarConstant && FEOrder == 0)
                    {
                        // 暫定: ScalarLagrange 1次
                        isUseVertex = true;
                        elemNodeCnt = 3;
                        edgeNos = new int[3][]
                        {
                            new int[] { 0, 1 },
                            new int[] { 1, 2 },
                            new int[] { 2, 0 }
                        };
                    }
                    else if (FEType == FiniteElementType.ScalarBell && FEOrder == 5)
                    {
                        edgeNos = new int[3][]
                        {
                            new int[] { 0, 1 },
                            new int[] { 1, 2 },
                            new int[] { 2, 0 }
                        };
                    }
                    else if (FEType == FiniteElementType.ScalarHermite && FEOrder == 3)
                    {
                        edgeNos = new int[3][]
                        {
                            new int[] { 0, 1 },
                            new int[] { 1, 2 },
                            new int[] { 2, 0 }
                        };
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(false);
                    }
                    int elemEdgeCnt = edgeNos.Length;
                    for (int eIndex = 0; eIndex < elemEdgeCnt; eIndex++)
                    {
                        int iNode1 = edgeNos[eIndex][0];
                        int iNode2 = edgeNos[eIndex][1];
                        int coId1;
                        int coId2;
                        if (isUseVertex)
                        {
                            // 暫定:
                            coId1 = triFE.VertexCoordIds[iNode1];
                            coId2 = triFE.VertexCoordIds[iNode2];
                        }
                        else
                        {
                            // 通常
                            coId1 = triFE.NodeCoordIds[iNode1];
                            coId2 = triFE.NodeCoordIds[iNode2];
                        }

                        int v1 = coId1;
                        int v2 = coId2;
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
            System.Diagnostics.Debug.Assert(Co2LineFE.Count == 0);
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

        // 接触解析のMaster/Slave線要素を準備する
        private void SetupContactMasterSlaveLineElements2D(FEWorld world)
        {
            var mesh = world.Mesh;

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
                if (ContactSlaveCadIds.Contains(cadId))
                {
                    // Slave上の線要素
                    ContactSlaveFEIds.Add(feId);
                }
                if (ContactMasterCadIds.Contains(cadId))
                {
                    // Master上の線要素
                    ContactMasterFEIds.Add(feId);
                }
            }
        }
    }
}
