using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Mesher2D
    {
        public bool MakeMeshFrom3D(
            Cad3DToXY cad,
            IList<OpenTK.Vector2d> initialMeshVecs,
            IList<MeshVertex2D> initialMeshVertexs,
            IList<MeshBarArray2D> initialMeshBarArrays)
        {
            Cad = cad;

            IList<uint> meshingLoopIds = new List<uint>();
            IList<double> meshingLoopELens = new List<double>();
            foreach (var cadIdELen in MeshingLoopCadIdELens)
            {
                uint lId = cadIdELen.CadId;
                double eLen = cadIdELen.ELen;
                if (!cad.IsElementId(CadElementType.Loop, lId))
                {
                    continue;
                }
                meshingLoopIds.Add(lId);
                meshingLoopELens.Add(eLen);
            }
            IList<uint> meshingEdgeIds = new List<uint>();
            IList<double> meshingEdgeELens = new List<double>();
            foreach (var cadIdELen in MeshingEdgeCadIdELens)
            {
                uint eId = cadIdELen.CadId;
                double eLen = cadIdELen.ELen;
                if (!cad.IsElementId(CadElementType.Edge, eId))
                {
                    continue;
                }
                meshingEdgeIds.Add(eId);
                meshingEdgeELens.Add(eLen);
            }

            if (MeshingMode == 0)
            {
                return TessellationFrom3D(
                    cad, meshingLoopIds, initialMeshVecs, initialMeshVertexs, initialMeshBarArrays);
            }
            else if (MeshingMode == 1)
            {
                return MakeMeshElemLengthFrom3D(
                    cad,
                    meshingLoopIds, meshingLoopELens, meshingEdgeIds, meshingEdgeELens,
                    initialMeshVecs, initialMeshVertexs, initialMeshBarArrays);
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return false;
        }

        private bool TessellationFrom3D(
            Cad3DToXY cad,
            IList<uint> loopIds,
            IList<OpenTK.Vector2d> initialMeshVec2Ds,
            IList<MeshVertex2D> initialMeshVertex2Ds,
            IList<MeshBarArray2D> initialMeshBarArray2Ds)
        {
            //-------------------------------
            // 分割済みのメッシュ情報を復元する
            Dictionary<uint, uint> initialVertexVecId = new Dictionary<uint, uint>();
            {
                // VEC
                foreach (OpenTK.Vector2d vec in initialMeshVec2Ds)
                {
                    Vecs.Add(vec);
                }

                // VERTEX
                foreach (MeshVertex2D vertex in initialMeshVertex2Ds)
                {
                    uint addId = GetFreeObjectId();
                   {
                        MeshVertex2D tmpVer = new MeshVertex2D();
                        tmpVer.Id = addId;
                        tmpVer.VCadId = vertex.VCadId;
                        tmpVer.V = vertex.V;
                        Vertexs.Add(tmpVer);

                        initialVertexVecId[tmpVer.Id] = tmpVer.V;
                    }
                    {
                        int typeLocCnt = TypeLocs.Count;
                        for (int i = typeLocCnt; i < addId + 1; i++)
                        {
                            var typeLoc = new MeshTypeLoc { Type = 0, Loc = -1 };
                            TypeLocs.Add(typeLoc);
                        }
                        TypeLocs[(int)addId].Loc = Vertexs.Count - 1;
                        TypeLocs[(int)addId].Type = 0; // VERTEX
                    }
                }

                // EDGE
                foreach (MeshBarArray2D barArray in initialMeshBarArray2Ds)
                {
                    uint sVId;
                    uint eVId;
                    cad.GetEdgeVertexId(barArray.ECadId, out sVId, out eVId);
                    //---------------------
                    // EDGEの始点、終点のVERTEXの要素IDを取得
                    uint sMeshId;
                    uint eMeshId;
                    {
                        uint loc;
                        uint type;
                        if (!FindElemLocTypeFromCadIdType(CadElementType.Vertex, sVId, out loc, out type))
                        {
                            System.Diagnostics.Debug.Assert(false);
                        }
                        System.Diagnostics.Debug.Assert(type == 0 && loc < Vertexs.Count);
                        MeshVertex2D sVer = Vertexs[(int)loc];
                        sMeshId = sVer.Id;
                        System.Diagnostics.Debug.Assert(sVer.V == initialVertexVecId[sMeshId]);
                        if (!FindElemLocTypeFromCadIdType(CadElementType.Vertex, eVId, out loc, out type))
                        {
                            System.Diagnostics.Debug.Assert(false);
                        }
                        System.Diagnostics.Debug.Assert(type == 0 && loc < Vertexs.Count);
                        MeshVertex2D eVer = Vertexs[(int)loc];
                        eMeshId = eVer.Id;
                        System.Diagnostics.Debug.Assert(eVer.V == initialVertexVecId[eMeshId]);
                    }
                    //-----------------------

                    uint addId = GetFreeObjectId();

                    IList<MeshBar> tmpBars = new List<MeshBar>();
                    foreach (MeshBar bar in barArray.Bars)
                    {
                        MeshBar tmpBar = new MeshBar();
                        tmpBar.V[0] = bar.V[0];
                        tmpBar.V[1] = bar.V[1];
                        tmpBar.S2[0] = 0;
                        tmpBar.S2[1] = 0;
                        tmpBar.R2[0] = 0;
                        tmpBar.R2[1] = 0;
                        tmpBars.Add(bar);
                    }
                    {
                        MeshBarArray2D tmpBarArray = new MeshBarArray2D();
                        tmpBarArray.Id = addId;
                        tmpBarArray.ECadId = barArray.ECadId;
                        tmpBarArray.SEId[0] = sMeshId;
                        tmpBarArray.SEId[1] = eMeshId;
                        tmpBarArray.LRId[0] = 0;
                        tmpBarArray.LRId[1] = 0;
                        tmpBarArray.Bars = tmpBars;
                        BarArrays.Add(tmpBarArray);
                    }
                    {
                        int typeLocCnt = TypeLocs.Count;
                        for (int i = typeLocCnt; i < addId + 1; i++)
                        {
                            var typeLoc = new MeshTypeLoc { Type = 0, Loc = -1 };
                            TypeLocs.Add(typeLoc);
                        }
                        TypeLocs[(int)addId].Loc = BarArrays.Count - 1;
                        TypeLocs[(int)addId].Type = 1; // BAR
                    }
                }
            }
            //-------------------------------

            //-------------------------------
            // 分割はここから
            // VERTEX
            {
                IList<uint> vIds = cad.GetElementIds(CadElementType.Vertex);
                for (uint iV = 0; iV < vIds.Count; iV++)
                {
                    uint vId = vIds[(int)iV];
                    if (GetIdFromCadId(vId, CadElementType.Vertex) != 0)
                    {
                        // 登録済み
                        continue;
                    }
                    uint addId = GetFreeObjectId();
                    OpenTK.Vector2d vec = cad.GetVertexCoord(vId);
                    Vecs.Add(vec);
                    {
                        MeshVertex2D tmpVer = new MeshVertex2D();
                        tmpVer.Id = addId;
                        tmpVer.VCadId = vId;
                        tmpVer.V = (uint)(Vecs.Count - 1);
                        Vertexs.Add(tmpVer);
                    }
                    {
                        int typeLocCnt = TypeLocs.Count;
                        for (int i = typeLocCnt; i < addId + 1; i++)
                        {
                            var typeLoc = new MeshTypeLoc { Type = 0, Loc = -1 };
                            TypeLocs.Add(typeLoc);
                        }
                        TypeLocs[(int)addId].Loc = Vertexs.Count - 1;
                        TypeLocs[(int)addId].Type = 0; // VERTEX
                    }
                    System.Diagnostics.Debug.Assert(CheckMesh() == 0);
                }
            }

            // EDGE
            {
                IList<uint> eIds = cad.GetElementIds(CadElementType.Edge);
                for (uint iE = 0; iE < eIds.Count; iE++)
                {
                    uint eId = eIds[(int)iE];

                    if (GetIdFromCadId(eId, CadElementType.Edge) != 0)
                    {
                        // 登録済み
                        continue;
                    }

                    TessellateEdge(cad, eId);

                    System.Diagnostics.Debug.Assert(CheckMesh() == 0);
                }
            }

            // LOOP
            {
                for (uint iL = 0; iL < loopIds.Count; iL++)
                {
                    uint lId = loopIds[(int)iL];

                    TessellateLoop(cad, lId);

                    System.Diagnostics.Debug.Assert(CheckMesh() == 0);
                }
            }

            MakeIncludeRelation(cad);

            return true;
        }

        private bool MakeMeshElemLengthFrom3D(
            Cad3DToXY cad,
            IList<uint> loopIds, IList<double> loopELens, IList<uint> edgeIds, IList<double> edgeELens,
            IList<OpenTK.Vector2d> initialMeshVec2Ds,
            IList<MeshVertex2D> initialMeshVertex2Ds,
            IList<MeshBarArray2D> initialMeshBarArray2Ds)
        {
            System.Diagnostics.Debug.Assert(loopIds.Count == loopELens.Count);
            System.Diagnostics.Debug.Assert(edgeIds.Count == edgeELens.Count);

            ClearMeshData();

            // ループに使われている頂点
            IList<uint> vtxFlgs = new List<uint>();
            for (int iLId = 0; iLId < loopIds.Count; iLId++)
            {
                uint lId = loopIds[iLId];
                LoopEdgeItr lItr = cad.GetLoopEdgeItr(lId);
                while (true)
                {
                    for (; !lItr.IsEnd(); lItr.Next())
                    {
                        uint vId = lItr.GetVertexId();
                        if (vtxFlgs.Count <= vId)
                        {
                            int cnt = vtxFlgs.Count;
                            for (int iTmp = cnt; iTmp < vId + 1; iTmp++)
                            {
                                vtxFlgs.Add(0);
                            }
                        }
                        vtxFlgs[(int)vId] = 1;
                    }
                    if (!lItr.ShiftChildLoop())
                    {
                        break;
                    }
                }
            }
            // 辺に使われている頂点
            for (int iEId = 0; iEId < edgeIds.Count; iEId++)
            {
                uint eId = edgeIds[iEId];
                uint sVId;
                uint eVId;
                cad.GetEdgeVertexId(eId, out sVId, out eVId);
                uint[] vIds = { sVId, eVId };
                foreach (uint vId in vIds)
                {
                    if (vtxFlgs.Count <= vId)
                    {
                        int cnt = vtxFlgs.Count;
                        for (int iTmp = cnt; iTmp < vId + 1; iTmp++)
                        {
                            vtxFlgs.Add(0);
                        }
                    }
                    vtxFlgs[(int)vId] = 1;
                }
            }

            //-------------------------------
            // 分割済みのメッシュ情報を復元する
            Dictionary<uint, uint> initialVertexVecId = new Dictionary<uint, uint>();
            {
                // VEC
                foreach (OpenTK.Vector2d vec in initialMeshVec2Ds)
                {
                    Vecs.Add(vec);
                }

                // VERTEX
                foreach (MeshVertex2D vertex in initialMeshVertex2Ds)
                {
                    uint vId = vertex.VCadId;
                    if (vtxFlgs[(int)vId] == 0)
                    {
                        // 未使用
                        continue;
                    }
                    uint addId = GetFreeObjectId();
                    {
                        MeshVertex2D tmpVer = new MeshVertex2D();
                        tmpVer.Id = addId;
                        tmpVer.VCadId = vertex.VCadId;
                        tmpVer.V = vertex.V;
                        Vertexs.Add(tmpVer);

                        initialVertexVecId[tmpVer.Id] = tmpVer.V;
                    }
                    {
                        int typeLocCnt = TypeLocs.Count;
                        for (int i = typeLocCnt; i < addId + 1; i++)
                        {
                            var typeLoc = new MeshTypeLoc { Type = 0, Loc = -1 };
                            TypeLocs.Add(typeLoc);
                        }
                        TypeLocs[(int)addId].Loc = Vertexs.Count - 1;
                        TypeLocs[(int)addId].Type = 0; // VERTEX
                    }
                }

                // EDGE
                foreach (MeshBarArray2D barArray in initialMeshBarArray2Ds)
                {
                    uint sVId;
                    uint eVId;
                    cad.GetEdgeVertexId(barArray.ECadId, out sVId, out eVId);
                    if (vtxFlgs[(int)sVId] == 0 || vtxFlgs[(int)eVId] == 0)
                    {
                        // 未使用
                        continue;
                    }
                    //---------------------
                    // EDGEの始点、終点のVERTEXの要素IDを取得
                    uint sMeshId;
                    uint eMeshId;
                    {
                        uint loc;
                        uint type;
                        if (!FindElemLocTypeFromCadIdType(CadElementType.Vertex, sVId, out loc, out type))
                        {
                            System.Diagnostics.Debug.Assert(false);
                        }
                        System.Diagnostics.Debug.Assert(type == 0 && loc < Vertexs.Count);
                        MeshVertex2D sVer = Vertexs[(int)loc];
                        sMeshId = sVer.Id;
                        System.Diagnostics.Debug.Assert(sVer.V == initialVertexVecId[sMeshId]);
                        if (!FindElemLocTypeFromCadIdType(CadElementType.Vertex, eVId, out loc, out type))
                        {
                            System.Diagnostics.Debug.Assert(false);
                        }
                        System.Diagnostics.Debug.Assert(type == 0 && loc < Vertexs.Count);
                        MeshVertex2D eVer = Vertexs[(int)loc];
                        eMeshId = eVer.Id;
                        System.Diagnostics.Debug.Assert(eVer.V == initialVertexVecId[eMeshId]);
                    }
                    //-----------------------

                    uint addId = GetFreeObjectId();

                    IList<MeshBar> tmpBars = new List<MeshBar>();
                    foreach (MeshBar bar in barArray.Bars)
                    {
                        MeshBar tmpBar = new MeshBar();
                        tmpBar.V[0] = bar.V[0];
                        tmpBar.V[1] = bar.V[1];
                        tmpBar.S2[0] = 0;
                        tmpBar.S2[1] = 0;
                        tmpBar.R2[0] = 0;
                        tmpBar.R2[1] = 0;
                        tmpBars.Add(bar);
                    }
                    {
                        MeshBarArray2D tmpBarArray = new MeshBarArray2D();
                        tmpBarArray.Id = addId;
                        tmpBarArray.ECadId = barArray.ECadId;
                        tmpBarArray.SEId[0] = sMeshId;
                        tmpBarArray.SEId[1] = eMeshId;
                        tmpBarArray.LRId[0] = 0;
                        tmpBarArray.LRId[1] = 0;
                        tmpBarArray.Bars = tmpBars;
                        BarArrays.Add(tmpBarArray);
                    }
                    {
                        int typeLocCnt = TypeLocs.Count;
                        for (int i = typeLocCnt; i < addId + 1; i++)
                        {
                            var typeLoc = new MeshTypeLoc { Type = 0, Loc = -1 };
                            TypeLocs.Add(typeLoc);
                        }
                        TypeLocs[(int)addId].Loc = BarArrays.Count - 1;
                        TypeLocs[(int)addId].Type = 1; // BAR
                    }
                }
            }
            //-------------------------------

            //-------------------------------
            // 分割はここから
            // VERTEX
            // フラグを立てた頂点のメッシュオブジェクトを生成する
            for (uint vId = 0; vId < vtxFlgs.Count; vId++)
            {
                if (vtxFlgs[(int)vId] == 0)
                {
                    continue;
                }
                if (GetIdFromCadId(vId, CadElementType.Vertex) != 0)
                {
                    // 登録済み
                    continue;
                }
                uint addId = GetFreeObjectId();
                OpenTK.Vector2d vec = cad.GetVertexCoord(vId);
                Vecs.Add(vec);
                {
                    MeshVertex2D tmpVer = new MeshVertex2D();
                    tmpVer.Id = addId;
                    tmpVer.VCadId = vId;
                    tmpVer.Layer = cad.GetLayer(CadElementType.Vertex, vId);
                    tmpVer.V = (uint)(Vecs.Count - 1);

                    Vertexs.Add(tmpVer);

                }
                {
                    int typeLocCnt = TypeLocs.Count;
                    for (int iTmp = typeLocCnt; iTmp < addId + 1; iTmp++)
                    {
                        var typeLoc = new MeshTypeLoc { Type = 0, Loc = -1 };
                        TypeLocs.Add(typeLoc);
                    }
                    TypeLocs[(int)addId].Loc = Vertexs.Count - 1;
                    TypeLocs[(int)addId].Type = 0; // VERTEX
                }
            }
            System.Diagnostics.Debug.Assert(CheckMesh() == 0);

            // EDGE(1)
            for (int iEId = 0; iEId < edgeIds.Count; iEId++)
            {
                uint eId = edgeIds[iEId];
                double eLen = edgeELens[iEId];
                if (GetIdFromCadId(eId, CadElementType.Edge) != 0)
                {
                    // 登録済み
                    continue;
                }
                MakeMeshEdge(cad, eId, eLen);
                System.Diagnostics.Debug.Assert(CheckMesh() == 0);
            }

            // EDGE(2)
            for (int iLId = 0; iLId < loopIds.Count; iLId++)
            {
                // ループに必要な辺を作る
                uint lId = loopIds[iLId];
                double eLen = loopELens[iLId];
                for (LoopEdgeItr lItr = cad.GetLoopEdgeItr(lId); !lItr.IsChildEnd; lItr.ShiftChildLoop())
                {
                    for (lItr.Begin(); !lItr.IsEnd(); lItr.Next())
                    {
                        uint eId;
                        bool isSameDir;
                        if (!lItr.GetEdgeId(out eId, out isSameDir))
                        {
                            continue;
                        }
                        if (GetIdFromCadId(eId, CadElementType.Edge) != 0)
                        {
                            // 登録済み
                            continue;
                        }
                        MakeMeshEdge(cad, eId, eLen);
                        System.Diagnostics.Debug.Assert(CheckMesh() == 0);
                    }
                }
            }

            // LOOP
            for (int iLId = 0; iLId < loopIds.Count; iLId++)
            {
                // ループを作る
                uint lId = loopIds[iLId];
                double eLen = loopELens[iLId];
                System.Diagnostics.Debug.Assert(eLen > 0.0);
                MakeMeshLoop(cad, lId, eLen);

                System.Diagnostics.Debug.Assert(CheckMesh() == 0);
            }

            MakeIncludeRelation(cad);
            return true;
        }
    }
}
