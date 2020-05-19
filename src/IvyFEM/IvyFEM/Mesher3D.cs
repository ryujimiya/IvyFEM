using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Mesher3D : IMesher
    {
        public uint Dimension => 3;
        public Cad3D Cad { get; private set; } = null;
        private IList<CadIdELen> MeshingLoopCadIdELens = new List<CadIdELen>();
        private IList<CadIdELen> MeshingEdgeCadIdELens = new List<CadIdELen>();
        private uint MeshingMode = 0; // 0: for cad 1: for mesh

        private IList<MeshTypeLoc> TypeLocs = new List<MeshTypeLoc>();

        private IList<MeshVertex3D> Vertexs = new List<MeshVertex3D>();
        private IList<MeshBarArray3D> BarArrays = new List<MeshBarArray3D>();
        private IList<MeshTriArray3D> TriArrays = new List<MeshTriArray3D>();

        private IList<OpenTK.Vector3d> Vecs = new List<OpenTK.Vector3d>();

        public Mesher3D()
        {
            MeshingMode = 1;
        }

        public Mesher3D(Cad3D cad)
        {
            Cad = cad;
            MeshingMode = 0;
            double eLen = 1.0;
            IList<uint> lIds = cad.GetElementIds(CadElementType.Loop);
            for (uint i = 0; i < lIds.Count; i++)
            {
                var cadIdELen = new CadIdELen { CadId = lIds[(int)i], ELen = eLen };
                MeshingLoopCadIdELens.Add(cadIdELen);
            }

            MakeMesh(Cad);
        }

        public Mesher3D(Cad3D cad, double eLen)
        {
            Cad = cad;
            MeshingMode = 1;

            IList<uint> lIds = cad.GetElementIds(CadElementType.Loop);
            for (int i = 0; i < lIds.Count; i++)
            {
                var cadIdELen = new CadIdELen { CadId = lIds[(int)i], ELen = eLen };
                MeshingLoopCadIdELens.Add(cadIdELen);
            }

            MakeMesh(Cad);
        }

        public Mesher3D(Mesher3D src)
        {
            Clear();
            Cad = src.Cad;
            MeshingLoopCadIdELens = new List<CadIdELen>();
            foreach (var srcCadIdELen in src.MeshingLoopCadIdELens)
            {
                MeshingLoopCadIdELens.Add(new CadIdELen(srcCadIdELen));
            }
            MeshingEdgeCadIdELens = new List<CadIdELen>();
            foreach (var srcCadIdELen in src.MeshingEdgeCadIdELens)
            {
                MeshingEdgeCadIdELens.Add(new CadIdELen(srcCadIdELen));
            }
            MeshingMode = src.MeshingMode;

            TypeLocs = new List<MeshTypeLoc>();
            foreach (var srcTypeLoc in src.TypeLocs)
            {
                TypeLocs.Add(new MeshTypeLoc(srcTypeLoc));
            }

            Vertexs = new List<MeshVertex3D>(src.Vertexs);
            BarArrays = new List<MeshBarArray3D>(src.BarArrays);
            TriArrays = new List<MeshTriArray3D>(src.TriArrays);

            Vecs = new List<OpenTK.Vector3d>(src.Vecs);
        }

        public void Clear()
        {
            MeshingLoopCadIdELens.Clear();
            MeshingEdgeCadIdELens.Clear();
            MeshingMode = 1;
            ClearMeshData();
        }

        private void ClearMeshData()
        {
            TypeLocs.Clear();

            Vertexs.Clear();
            BarArrays.Clear();
            TriArrays.Clear();

            Vecs.Clear();
        }

        public void AddMeshingLoopCadId(uint lCadId, double eLen)
        {
            int index = IndexOfMeshingLoopCadId(lCadId);
            if (index != -1)
            {
                return;
            }
            MeshingLoopCadIdELens.Add(new CadIdELen { CadId = lCadId, ELen = eLen });
        }

        public bool IsMeshingLoopCadId(uint lCadId)
        {
            int index = IndexOfMeshingLoopCadId(lCadId);
            return index != -1;
        }

        private int IndexOfMeshingLoopCadId(uint lCadId)
        {
            int index = -1;
            for (int i = 0; i < MeshingLoopCadIdELens.Count; i++)
            {
                var cadIdELen = MeshingLoopCadIdELens[i];
                if (cadIdELen.CadId == lCadId)
                {
                    index = i;
                    break;
                }
            }
            return index;
        }

        public void RemoveMeshingLoopCadId(uint lCadId)
        {
            int index = IndexOfMeshingLoopCadId(lCadId);
            if (index == -1)
            {
                return;
            }
            MeshingLoopCadIdELens.RemoveAt(index);
        }

        public void ClearMeshingLoopCadIds()
        {
            MeshingLoopCadIdELens.Clear();
        }

        public IList<uint> GetMeshingLoopCadIds()
        {
            IList<uint> lIds = new List<uint>();
            for (int i = 0; i < MeshingLoopCadIdELens.Count; i++)
            {
                var cadIdELen = MeshingLoopCadIdELens[i];
                lIds.Add(cadIdELen.CadId);
            }
            return lIds;
        }

        public IList<double> GetMeshingLoopCadELens()
        {
            IList<double> eLens = new List<double>();
            for (int i = 0; i < MeshingLoopCadIdELens.Count; i++)
            {
                var cadIdELen = MeshingLoopCadIdELens[i];
                eLens.Add(cadIdELen.ELen);
            }
            return eLens;
        }

        public void AddMeshingEdgeCadId(uint eCadId, double eLen)
        {
            int index = IndexOfMeshingEdgeCadId(eCadId);
            if (index != -1)
            {
                return;
            }
            MeshingEdgeCadIdELens.Add(new CadIdELen { CadId = eCadId, ELen = eLen });
        }

        public bool IsMeshingEdgeCadId(uint eCadId)
        {
            int index = IndexOfMeshingEdgeCadId(eCadId);
            return index != -1;
        }

        private int IndexOfMeshingEdgeCadId(uint eCadId)
        {
            int index = -1;
            for (int i = 0; i < MeshingEdgeCadIdELens.Count; i++)
            {
                var cadIdELen = MeshingEdgeCadIdELens[i];
                if (cadIdELen.CadId == eCadId)
                {
                    index = i;
                    break;
                }
            }
            return index;
        }

        public void RemoveMeshingEdgeCadId(uint eCadId)
        {
            int index = IndexOfMeshingLoopCadId(eCadId);
            if (index == -1)
            {
                return;
            }
            MeshingEdgeCadIdELens.RemoveAt(index);
        }

        public void ClearMeshingEdgeCadIds()
        {
            MeshingEdgeCadIdELens.Clear();
        }

        public IList<uint> GetMeshingEdgeCadIds()
        {
            IList<uint> eIds = new List<uint>();
            for (int i = 0; i < MeshingEdgeCadIdELens.Count; i++)
            {
                var cadIdELen = MeshingEdgeCadIdELens[i];
                eIds.Add(cadIdELen.CadId);
            }
            return eIds;
        }

        public IList<double> GetMeshingEdgeCadELens()
        {
            IList<double> eLens = new List<double>();
            for (int i = 0; i < MeshingEdgeCadIdELens.Count; i++)
            {
                var cadIdELen = MeshingEdgeCadIdELens[i];
                eLens.Add(cadIdELen.ELen);
            }
            return eLens;
        }

        public IList<MeshTriArray3D> GetTriArrays()
        {
            return TriArrays;
        }

        public IList<MeshBarArray3D> GetBarArrays()
        {
            return BarArrays;
        }

        public IList<MeshVertex3D> GetVertexs()
        {
            return Vertexs;
        }

        public IList<OpenTK.Vector3d> GetVectors()
        {
            return Vecs;
        }

        public void GetCoords(out IList<double> coord)
        {
            coord = new List<double>();
            uint nodeCnt = (uint)Vecs.Count;
            uint dim = 3;
            for (int i = 0; i < nodeCnt * dim; i++)
            {
                coord.Add(0.0);
            }
            for (int iNode = 0; iNode < nodeCnt; iNode++)
            {
                coord[(int)(iNode * dim)] = Vecs[iNode].X;
                coord[(int)(iNode * dim + 1)] = Vecs[iNode].Y;
                coord[(int)(iNode * dim + 2)] = Vecs[iNode].Z;
            }
        }

        public IList<uint> GetIds()
        {
            IList<uint> ids = new List<uint>();
            for (uint id = 1; id < TypeLocs.Count; id++)
            {
                if (TypeLocs[(int)id].Loc == -1)
                {
                    continue;
                }
                ids.Add(id);
            }
            return ids;
        }

        public bool IsId(uint id)
        {
            if (TypeLocs.Count <= id)
            {
                return false;
            }
            int loc = TypeLocs[(int)id].Loc;
            if (loc == -1)
            {
                return false;
            }
            int type = TypeLocs[(int)id].Type;
            System.Diagnostics.Debug.Assert(type >= 0);
            System.Diagnostics.Debug.Assert(loc >= 0);
            if (type == 0)
            {
                System.Diagnostics.Debug.Assert(Vertexs.Count > loc);
                System.Diagnostics.Debug.Assert(Vertexs[loc].Id == id);
            }
            else if (type == 1)
            {
                System.Diagnostics.Debug.Assert(BarArrays.Count > loc);
                System.Diagnostics.Debug.Assert(BarArrays[loc].Id == id);
            }
            else if (type == 2)
            {
                System.Diagnostics.Debug.Assert(TriArrays.Count > loc);
                System.Diagnostics.Debug.Assert(TriArrays[loc].Id == id);
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return true;
        }

        public IList<uint> GetIncludeElementIds(uint id)
        {
            IList<uint> ids = new List<uint>();
            if (id >= TypeLocs.Count)
            {
                return ids;
            }
            if (TypeLocs[(int)id].Loc == -1)
            {
                return ids;
            }
            return TypeLocs[(int)id].IncludeRelations;
        }

        public void GetInfo(uint id, out uint cadId, out int layer)
        {
            cadId = 0;
            layer = 0;

            int type = TypeLocs[(int)id].Type;
            int loc = TypeLocs[(int)id].Loc;
            if (type == 0)
            {
                cadId = Vertexs[loc].VCadId;
            }
            else if (type == 1)
            {
                cadId = BarArrays[loc].ECadId;
            }
            else if (type == 2)
            {
                cadId = TriArrays[loc].LCadId;
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
        }

        public void GetConnectivity(uint id, out MeshType meshType, out int[] vertexs)
        {
            meshType = MeshType.NotSet;
            vertexs = null;

            System.Diagnostics.Debug.Assert(IsId(id));
            uint elemNodeCnt;
            uint elemCnt;
            int type = TypeLocs[(int)id].Type;
            int loc = TypeLocs[(int)id].Loc;
            System.Diagnostics.Debug.Assert(type != -1 && loc != -1);
            if (type == 0)
            {
                meshType = MeshType.Vertex;
                elemNodeCnt = 1;
                elemCnt = 1;
                vertexs = new int[elemNodeCnt * elemCnt];
                vertexs[0] = (int)Vertexs[loc].V;
            }
            else if (type == 1)
            {
                meshType = MeshType.Bar;
                elemNodeCnt = 2;
                IList<MeshBar> bars = BarArrays[loc].Bars;
                elemCnt = (uint)bars.Count;
                vertexs = new int[elemNodeCnt * elemCnt];
                for (int iElem = 0; iElem < elemCnt; iElem++)
                {
                    vertexs[iElem * 2] = (int)bars[iElem].V[0];
                    vertexs[iElem * 2 + 1] = (int)bars[iElem].V[1];
                }
            }
            else if (type == 2)
            {
                meshType = MeshType.Tri;
                elemNodeCnt = 3;
                IList<MeshTri3D> tris = TriArrays[loc].Tris;
                elemCnt = (uint)tris.Count;
                vertexs = new int[elemNodeCnt * elemCnt];
                for (int iElem = 0; iElem < elemCnt; iElem++)
                {
                    vertexs[iElem * 3] = (int)tris[iElem].V[0];
                    vertexs[iElem * 3 + 1] = (int)tris[iElem].V[1];
                    vertexs[iElem * 3 + 2] = (int)tris[iElem].V[2];
                }
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
        }

        public bool GetMeshInfo(uint id,
            out uint elemCount, out MeshType meshType, out int loc, out uint cadId)
        {
            elemCount = 0;
            meshType = MeshType.NotSet;
            loc = 0;
            cadId = 0;

            if (!IsId(id))
            {
                return false;
            }
            System.Diagnostics.Debug.Assert(id < TypeLocs.Count);
            int type = TypeLocs[(int)id].Type;
            loc = TypeLocs[(int)id].Loc;
            System.Diagnostics.Debug.Assert(type >= 0);
            if (type == 0)
            {
                System.Diagnostics.Debug.Assert(loc < Vertexs.Count);
                MeshVertex3D vertex = Vertexs[loc];
                cadId = vertex.VCadId;
                elemCount = 1;
                meshType = MeshType.Vertex;
            }
            else if (type == 1)
            {
                System.Diagnostics.Debug.Assert(loc < BarArrays.Count);
                MeshBarArray3D barArray = BarArrays[loc];
                cadId = barArray.ECadId;
                elemCount = (uint)barArray.Bars.Count;
                meshType = MeshType.Bar;
            }
            else if (type == 2)
            {
                System.Diagnostics.Debug.Assert(loc < TriArrays.Count);
                MeshTriArray3D triArray = TriArrays[loc];
                cadId = triArray.LCadId;
                elemCount = (uint)triArray.Tris.Count;
                meshType = MeshType.Tri;
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
                throw new NotImplementedException();
            }
            return true;
        }

        private uint FindMaxId()
        {
            uint maxId = 0;
            {
                for (uint iVer = 0; iVer < Vertexs.Count; iVer++)
                {
                    if (maxId < Vertexs[(int)iVer].Id)
                    {
                        maxId = Vertexs[(int)iVer].Id;
                    }
                }
                for (uint iBarArray = 0; iBarArray < BarArrays.Count; iBarArray++)
                {
                    if (maxId < BarArrays[(int)iBarArray].Id)
                    {
                        maxId = BarArrays[(int)iBarArray].Id;
                    }
                }
                for (uint iTriArray = 0; iTriArray < TriArrays.Count; iTriArray++)
                {
                    if (maxId < TriArrays[(int)iTriArray].Id)
                    {
                        maxId = TriArrays[(int)iTriArray].Id;
                    }
                }
            }
            return maxId;
        }

        private uint GetFreeObjectId()
        {
            uint maxId = FindMaxId();
            IList<uint> isUsedFlgs = new List<uint>();
            {
                for (uint iUse = 0; iUse < maxId + 1; iUse++)
                {
                    isUsedFlgs.Add(0);
                }
                for (uint iVer = 0; iVer < Vertexs.Count; iVer++)
                {
                    System.Diagnostics.Debug.Assert(isUsedFlgs[(int)Vertexs[(int)iVer].Id] == 0);
                    System.Diagnostics.Debug.Assert(Vertexs[(int)iVer].Id >= 1 && Vertexs[(int)iVer].Id <= maxId);
                    isUsedFlgs[(int)Vertexs[(int)iVer].Id] = 1;
                }
                for (uint iBarArray = 0; iBarArray < BarArrays.Count; iBarArray++)
                {
                    System.Diagnostics.Debug.Assert(isUsedFlgs[(int)BarArrays[(int)iBarArray].Id] == 0);
                    System.Diagnostics.Debug.Assert(BarArrays[(int)iBarArray].Id >= 1 && BarArrays[(int)iBarArray].Id <= maxId);
                    isUsedFlgs[(int)BarArrays[(int)iBarArray].Id] = 1;
                }
                for (uint iTriArray = 0; iTriArray < TriArrays.Count; iTriArray++)
                {
                    System.Diagnostics.Debug.Assert(isUsedFlgs[(int)TriArrays[(int)iTriArray].Id] == 0);
                    System.Diagnostics.Debug.Assert(TriArrays[(int)iTriArray].Id >= 1 && TriArrays[(int)iTriArray].Id <= maxId);
                    isUsedFlgs[(int)TriArrays[(int)iTriArray].Id] = 1;
                }
            }
            System.Diagnostics.Debug.Assert(isUsedFlgs[0] == 0);
            for (uint i = 1; i < isUsedFlgs.Count; i++)
            {
                if (isUsedFlgs[(int)i] == 0)
                {
                    return i;
                }
            }
            return maxId + 1;
        }

        public uint GetIdFromCadId(uint cadId, CadElementType cadType)
        {
            switch (cadType)
            {
                case CadElementType.Vertex:
                    for (uint iVer = 0; iVer < Vertexs.Count; iVer++)
                    {
                        if (Vertexs[(int)iVer].VCadId == cadId)
                        {
                            uint meshId = Vertexs[(int)iVer].Id;
                            System.Diagnostics.Debug.Assert(TypeLocs[(int)meshId].Loc == iVer);
                            System.Diagnostics.Debug.Assert(TypeLocs[(int)meshId].Type == 0); // VERTEX
                            return meshId;
                        }
                    }
                    break;

                case CadElementType.Edge:
                    for (uint iBar = 0; iBar < BarArrays.Count; iBar++)
                    {
                        if (BarArrays[(int)iBar].ECadId == cadId)
                        {
                            uint meshId = BarArrays[(int)iBar].Id;
                            System.Diagnostics.Debug.Assert(TypeLocs[(int)meshId].Loc == iBar);
                            System.Diagnostics.Debug.Assert(TypeLocs[(int)meshId].Type == 1); // BAR
                            return meshId;
                        }
                    }
                    break;

                case CadElementType.Loop:
                    for (uint iTri = 0; iTri < TriArrays.Count; iTri++)
                    {
                        if (TriArrays[(int)iTri].LCadId == cadId)
                        {
                            uint meshId = TriArrays[(int)iTri].Id;
                            System.Diagnostics.Debug.Assert(TypeLocs[(int)meshId].Loc == iTri);
                            System.Diagnostics.Debug.Assert(TypeLocs[(int)meshId].Type == 2); // TRI
                            return meshId;
                        }
                    }
                    break;
                default:
                    System.Diagnostics.Debug.Assert(false);
                    return 0;
            }
            return 0;
        }

        private bool FindElemLocTypeFromCadIdType(CadElementType cadType, uint cadId, out uint loc, out uint type)
        {
            loc = 0;
            type = 0;

            switch (cadType)
            {
                case CadElementType.Vertex:
                    for (uint iVer = 0; iVer < Vertexs.Count; iVer++)
                    {
                        if (Vertexs[(int)iVer].VCadId == cadId)
                        {
                            loc = iVer;
                            type = 0; // VERTEX
                            uint meshId = Vertexs[(int)iVer].Id;
                            System.Diagnostics.Debug.Assert(TypeLocs[(int)meshId].Loc == loc);
                            System.Diagnostics.Debug.Assert(TypeLocs[(int)meshId].Type == type);
                            return true;
                        }
                    }
                    break;

                case CadElementType.Edge:
                    for (uint iBar = 0; iBar < BarArrays.Count; iBar++)
                    {
                        if (BarArrays[(int)iBar].ECadId == cadId)
                        {
                            loc = iBar;
                            type = 1; // BAR
                            uint meshId = BarArrays[(int)iBar].Id;
                            System.Diagnostics.Debug.Assert(TypeLocs[(int)meshId].Loc == loc);
                            System.Diagnostics.Debug.Assert(TypeLocs[(int)meshId].Type == type);
                            return true;
                        }
                    }
                    break;

                case CadElementType.Loop:
                    for (uint iTri = 0; iTri < TriArrays.Count; iTri++)
                    {
                        if (TriArrays[(int)iTri].LCadId == cadId)
                        {
                            loc = iTri;
                            type = 2; // TRI
                            uint meshId = TriArrays[(int)iTri].Id;
                            System.Diagnostics.Debug.Assert(TypeLocs[(int)meshId].Loc == loc);
                            System.Diagnostics.Debug.Assert(TypeLocs[(int)meshId].Type == type);
                            return true;
                        }
                    }
                    break;

                default:
                    System.Diagnostics.Debug.Assert(false);
                    break;
            }
            return false;
        }

        public void SetMeshingModeZero()
        {
            MeshingMode = 0;
        }

        public void SetMeshingModeElemLength()
        {
            MeshingMode = 1;
        }

        public bool MakeMesh(Cad3D cad)
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
                return Tessellation(cad, meshingLoopIds);
            }
            else if (MeshingMode == 1)
            {
                return MakeMeshElemLength(cad, meshingLoopIds, meshingLoopELens, meshingEdgeIds, meshingEdgeELens);
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return false;
        }

        private bool Tessellation(Cad3D cad, IList<uint> loopIds)
        {
            // VERTEX
            {
                IList<uint> vIds = cad.GetElementIds(CadElementType.Vertex);
                for (uint iV = 0; iV < vIds.Count; iV++)
                {
                    uint vId = vIds[(int)iV];
                    System.Diagnostics.Debug.Assert(GetIdFromCadId(vId, CadElementType.Vertex) == 0);
                    uint addId = GetFreeObjectId();
                    OpenTK.Vector3d vec = cad.GetVertexCoord(vId);
                    Vecs.Add(vec);
                    {
                        MeshVertex3D tmpVer = new MeshVertex3D();
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

        private bool TessellateEdge(Cad3D cad, uint eId)
        {
            uint sVId;
            uint eVId;
            System.Diagnostics.Debug.Assert(cad.IsElementId(CadElementType.Edge, eId));
            if (!cad.GetEdgeVertexId(eId, out sVId, out eVId))
            {
                System.Diagnostics.Debug.WriteLine("error edge : " + eId);
                System.Diagnostics.Debug.Assert(false);
            }

            uint iSP;
            uint iEP;
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
                MeshVertex3D sVer = Vertexs[(int)loc];
                iSP = sVer.V;
                sMeshId = sVer.Id;
                if (!FindElemLocTypeFromCadIdType(CadElementType.Vertex, eVId, out loc, out type))
                {
                    System.Diagnostics.Debug.Assert(false);
                }
                System.Diagnostics.Debug.Assert(type == 0 && loc < Vertexs.Count);

                MeshVertex3D eVer = Vertexs[(int)loc];
                iEP = eVer.V;
                eMeshId = eVer.Id;
            }

            IList<MeshBar> tmpBars = new List<MeshBar>();
            {
                int iBar = 0;
                tmpBars.Add(new MeshBar());
                tmpBars[iBar].V[0] = iSP;
                tmpBars[iBar].V[1] = iEP;
                tmpBars[iBar].S2[0] = 0;
                tmpBars[iBar].S2[1] = 0;
                tmpBars[iBar].R2[0] = 0;
                tmpBars[iBar].R2[1] = 0;
            }
            uint newElemId = GetFreeObjectId();
            uint iBarArray = (uint)BarArrays.Count;
            BarArrays.Add(new MeshBarArray3D());
            BarArrays[(int)iBarArray].Bars = tmpBars;
            BarArrays[(int)iBarArray].ECadId = eId;
            BarArrays[(int)iBarArray].Id = newElemId;

            int typeLocCnt = TypeLocs.Count;
            for (int i = typeLocCnt; i < newElemId + 1; i++)
            {
                var typeLoc = new MeshTypeLoc { Type = 0, Loc = -1 };
                TypeLocs.Add(typeLoc);
            }
            TypeLocs[(int)newElemId].Loc = (int)iBarArray;
            TypeLocs[(int)newElemId].Type = 1; // BAR

            return true;
        }

        private bool TessellateLoop(Cad3D cad, uint lId)
        {
            Cad3DToXY cad2D = new Cad3DToXY(cad, CadElementType.Loop, lId);

            IList<OpenTK.Vector2d> initialVec2Ds;
            Dictionary<uint, uint> initialVecId3D2D;
            Dictionary<uint, uint> initialVecId2D3D;
            IList< MeshVertex2D> initialVertex2Ds;
            IList<MeshBarArray2D> initialBarArray2Ds;
            ConvertDoneMeshToXY(
                cad, lId,
                cad2D,
                out initialVec2Ds,
                out initialVecId3D2D,
                out initialVecId2D3D,
                out initialVertex2Ds,
                out initialBarArray2Ds);

            Mesher2D mesh2D = new Mesher2D(
                cad2D, initialVec2Ds, initialVertex2Ds, initialBarArray2Ds);

            // VEC
            var vec2Ds = mesh2D.GetVectors();
            Dictionary<uint, uint> vecId2D3D = new Dictionary<uint, uint>(initialVecId2D3D);

            // VERTEX
            // することはない

            // EDGE
            var barArray2Ds = mesh2D.GetBarArrays();
            foreach (MeshBarArray2D barArray2D in barArray2Ds)
            {
                uint eId2D = barArray2D.ECadId;
                uint eId = cad2D.GetEdgeId3DFrom2D(eId2D);
                if (GetIdFromCadId(eId, CadElementType.Edge) != 0)
                {
                    // 登録済み
                    continue;
                }

                IList<MeshBar> tmpBars = new List<MeshBar>();

                int elementCnt = barArray2D.Bars.Count;
                for (int iBar = 0; iBar < elementCnt; iBar++)
                {
                    MeshBar tmpBar = new MeshBar();

                    MeshBar bar2D = barArray2D.Bars[iBar];
                    for (int i = 0; i < 2; i++)
                    {
                        uint ivec2D = bar2D.V[i];
                        uint ivec;
                        if (vecId2D3D.ContainsKey(ivec2D))
                        {
                            ivec = vecId2D3D[ivec2D];
                        }
                        else
                        {
                            OpenTK.Vector2d vec2D = vec2Ds[(int)ivec2D];
                            OpenTK.Vector3d vec = cad2D.UnProject(vec2D);
                            Vecs.Add(vec);
                            ivec = (uint)(Vecs.Count - 1);
                            vecId2D3D[ivec2D] = ivec;
                        }
                        tmpBar.V[i] = ivec;
                    }
                    tmpBars.Add(tmpBar);
                }

                uint newElemId = GetFreeObjectId();
                uint iBarArray = (uint)BarArrays.Count;
                BarArrays.Add(new MeshBarArray3D());
                BarArrays[(int)iBarArray].Bars = tmpBars;
                BarArrays[(int)iBarArray].ECadId = eId;
                BarArrays[(int)iBarArray].Id = newElemId;

                int typeLocCnt = TypeLocs.Count;
                for (int i = typeLocCnt; i < newElemId + 1; i++)
                {
                    var typeLoc = new MeshTypeLoc { Type = 0, Loc = -1 };
                    TypeLocs.Add(typeLoc);
                }
                TypeLocs[(int)newElemId].Loc = (int)iBarArray;
                TypeLocs[(int)newElemId].Type = 1; // BAR
            }

            // LOOP
            System.Diagnostics.Debug.Assert(GetIdFromCadId(lId, CadElementType.Loop) == 0); // 未登録でなければならない
            var triArray2Ds = mesh2D.GetTriArrays();
            foreach (MeshTriArray2D triArray2D in triArray2Ds)
            {
                IList<MeshTri3D> tmpTris = new List<MeshTri3D>();

                int elementCnt = triArray2D.Tris.Count;
                for (int iTri = 0; iTri < elementCnt; iTri++)
                {
                    MeshTri3D tmpTri = new MeshTri3D();

                    MeshTri2D tri2D = triArray2D.Tris[iTri];
                    for (int iTriNo = 0; iTriNo < 3; iTriNo++)
                    {
                        uint ivec2D = tri2D.V[iTriNo];
                        uint ivec;
                        if (vecId2D3D.ContainsKey(ivec2D))
                        {
                            ivec = vecId2D3D[ivec2D];
                        }
                        else
                        {
                            OpenTK.Vector2d vec2D = vec2Ds[(int)ivec2D];
                            OpenTK.Vector3d vec = cad2D.UnProject(vec2D);
                            Vecs.Add(vec);
                            ivec = (uint)(Vecs.Count - 1);
                            vecId2D3D[ivec2D] = ivec;
                        }
                        tmpTri.V[iTriNo] = ivec;
                    }
                    tmpTris.Add(tmpTri);
                }

                uint newElemId = GetFreeObjectId();
                uint iTriArray = (uint)TriArrays.Count;
                TriArrays.Add(new MeshTriArray3D());
                TriArrays[(int)iTriArray].Tris = tmpTris;
                TriArrays[(int)iTriArray].LCadId = lId;
                TriArrays[(int)iTriArray].Id = newElemId;

                int typeLocCnt = TypeLocs.Count;
                for (int i = typeLocCnt; i < newElemId + 1; i++)
                {
                    var typeLoc = new MeshTypeLoc { Type = 0, Loc = -1 };
                    TypeLocs.Add(typeLoc);
                }
                TypeLocs[(int)newElemId].Type = 2;   // TRI
                TypeLocs[(int)newElemId].Loc = (int)iTriArray;

                System.Diagnostics.Debug.Assert(CheckMesh() == 0);
            }

            return true;
        }

        private void ConvertDoneMeshToXY(
            Cad3D cad, uint lId,
            Cad3DToXY cad2D,
            out IList<OpenTK.Vector2d> initialVec2Ds,
            out Dictionary<uint, uint> initialVecId3D2D,
            out Dictionary<uint, uint> initialVecId2D3D,
            out IList<MeshVertex2D> initialVertex2Ds,
            out IList<MeshBarArray2D> initialBarArray2Ds)
        {
            //------------------------------------------
            // 分割済みのメッシュ(VEC,VERTEX,EDGE)
            IList<uint> vIds = new List<uint>();
            IList<uint> eIds = new List<uint>();
            for (LoopEdgeItr lItr = cad.GetLoopEdgeItr(lId); !lItr.IsChildEnd; lItr.ShiftChildLoop())
            {
                for (lItr.Begin(); !lItr.IsEnd(); lItr.Next())
                {
                    uint vId = lItr.GetVertexId();
                    if (vIds.IndexOf(vId) == -1)
                    {
                        vIds.Add(vId);
                    }

                    // EDGE
                    uint eId;
                    bool isSameDir;
                    lItr.GetEdgeId(out eId, out isSameDir);
                    if (!cad.IsElementId(CadElementType.Edge, eId))
                    {
                        continue;
                    }
                    if (eIds.IndexOf(eId) == -1)
                    {
                        eIds.Add(eId);
                    }
                }
            }
            initialVec2Ds = new List<OpenTK.Vector2d>();
            initialVecId3D2D = new Dictionary<uint, uint>();
            initialVecId2D3D = new Dictionary<uint, uint>();
            initialVertex2Ds = new List<MeshVertex2D>();
            foreach (uint vId in vIds)
            {
                MeshVertex3D vertex;
                {
                    uint elemId = GetIdFromCadId(vId, CadElementType.Vertex);
                    int loc = TypeLocs[(int)elemId].Loc;
                    System.Diagnostics.Debug.Assert(TypeLocs[(int)elemId].Type == 0); // VERTEX
                    vertex = Vertexs[(int)loc];
                }
                MeshVertex2D vertex2D = new MeshVertex2D();
                vertex2D.Id = 0; // 未定
                vertex2D.VCadId = cad2D.GetVertexId2DFrom3D(vId);
                {
                    uint ivec = vertex.V;
                    OpenTK.Vector3d vec = Vecs[(int)ivec];
                    uint ivec2D;
                    if (initialVecId3D2D.ContainsKey(ivec))
                    {
                        ivec2D = initialVecId3D2D[ivec];
                    }
                    else
                    {
                        OpenTK.Vector2d vec2D = cad2D.Project(vec);
                        initialVec2Ds.Add(vec2D);
                        ivec2D = (uint)(initialVec2Ds.Count - 1);
                        initialVecId3D2D[ivec] = ivec2D;
                        initialVecId2D3D[ivec2D] = ivec;
                    }
                    vertex2D.V = ivec2D;
                }
                initialVertex2Ds.Add(vertex2D);
            }
            initialBarArray2Ds = new List<MeshBarArray2D>();
            foreach (uint eId in eIds)
            {
                MeshBarArray3D barArray;
                {
                    uint elemId = GetIdFromCadId(eId, CadElementType.Edge);
                    int loc = TypeLocs[(int)elemId].Loc;
                    System.Diagnostics.Debug.Assert(TypeLocs[(int)elemId].Type == 1); // BAR
                    barArray = BarArrays[(int)loc];
                }
                MeshBarArray2D barArray2D = new MeshBarArray2D();
                barArray2D.Id = 0; // 未定
                barArray2D.ECadId = cad2D.GetEdgeId2DFrom3D(eId);
                barArray2D.SEId[0] = 0; // 未定
                barArray2D.SEId[1] = 0; // 未定
                for (int iBar = 0; iBar < barArray.Bars.Count; iBar++)
                {
                    MeshBar bar = barArray.Bars[iBar];
                    MeshBar bar2D = new MeshBar();
                    for (int iBarNo = 0; iBarNo < 2; iBarNo++)
                    {
                        uint ivec = bar.V[iBarNo];
                        OpenTK.Vector3d vec = Vecs[(int)ivec];
                        uint ivec2D;
                        if (initialVecId3D2D.ContainsKey(ivec))
                        {
                            ivec2D = initialVecId3D2D[ivec];
                        }
                        else
                        {
                            OpenTK.Vector2d vec2D = cad2D.Project(vec);
                            initialVec2Ds.Add(vec2D);
                            ivec2D = (uint)(initialVec2Ds.Count - 1);
                            initialVecId3D2D[ivec] = ivec2D;
                            initialVecId2D3D[ivec2D] = ivec;
                        }
                        bar2D.V[iBarNo] = ivec2D;
                    }
                    barArray2D.Bars.Add(bar2D);
                }
                initialBarArray2Ds.Add(barArray2D);
            }
        }

        private void MakeIncludeRelation(Cad3D cad)
        {
            if (TypeLocs.Count == 0)
            {
                return;
            }

            System.Diagnostics.Debug.Assert(FindMaxId() + 1 == TypeLocs.Count);
            foreach (var typeLoc in TypeLocs)
            {
                typeLoc.IncludeRelations.Clear();
            }

            IList<uint> lIds = cad.GetElementIds(CadElementType.Loop);
            for (uint iLId = 0; iLId < lIds.Count; iLId++)
            {
                uint lId = lIds[(int)iLId];
                uint triId = GetIdFromCadId(lId, CadElementType.Loop);
                if (!IsId(triId))
                {
                    continue;
                }
                LoopEdgeItr lItr = cad.GetLoopEdgeItr(lId);
                while (true)
                {
                    for (; !lItr.IsEnd(); lItr.Next())
                    {
                        uint vCadId = lItr.GetVertexId();
                        uint meshVId = GetIdFromCadId(vCadId, CadElementType.Vertex);
                        TypeLocs[(int)triId].IncludeRelations.Add(meshVId);
                        System.Diagnostics.Debug.Assert(IsId(meshVId));

                        uint eId;
                        bool isSameDir;
                        if (!lItr.GetEdgeId(out eId, out isSameDir))
                        {
                            continue;
                        }
                        uint barId = GetIdFromCadId(eId, CadElementType.Edge);
                        System.Diagnostics.Debug.Assert(IsId(barId));
                        TypeLocs[(int)triId].IncludeRelations.Add(barId);
                    }
                    if (!lItr.ShiftChildLoop())
                    {
                        break;
                    }
                }
            }

            IList<uint> eIds = cad.GetElementIds(CadElementType.Edge);
            for (uint iEId = 0; iEId < eIds.Count; iEId++)
            {
                uint eId = eIds[(int)iEId];
                System.Diagnostics.Debug.Assert(cad.IsElementId(CadElementType.Edge, eId));
                uint barId = GetIdFromCadId(eId, CadElementType.Edge);
                if (!IsId(barId))
                {
                    // 浮いている辺があって，辺メッシュが切られなかった場合
                    continue;
                }
                uint sVId;
                uint eVId;
                if (!cad.GetEdgeVertexId(eId, out sVId, out eVId))
                {
                    System.Diagnostics.Debug.Assert(false);
                }
                uint meshSVId = GetIdFromCadId(sVId, CadElementType.Vertex);
                uint meshEVId = GetIdFromCadId(eVId, CadElementType.Vertex);
                System.Diagnostics.Debug.Assert(IsId(meshSVId));
                System.Diagnostics.Debug.Assert(IsId(meshEVId));
                TypeLocs[(int)barId].IncludeRelations.Add(meshSVId);
                TypeLocs[(int)barId].IncludeRelations.Add(meshEVId);
            }
        }

        private int CheckMesh()
        {
            {
                uint maxId = 0;
                for (uint iVer = 0; iVer < Vertexs.Count; iVer++)
                {
                    uint id0 = Vertexs[(int)iVer].Id;
                    if (maxId < id0)
                    {
                        maxId = id0;
                    }
                    System.Diagnostics.Debug.Assert(TypeLocs.Count > id0);
                    System.Diagnostics.Debug.Assert(TypeLocs[(int)id0].Type == 0); // VERTEX
                    int loc0 = TypeLocs[(int)id0].Loc;
                    System.Diagnostics.Debug.Assert(loc0 == (int)iVer);
                    System.Diagnostics.Debug.Assert(Vertexs.Count > loc0);
                    MeshVertex3D ver0 = Vertexs[(int)loc0];
                    System.Diagnostics.Debug.Assert(ver0.Id == id0);
                }
                for (uint iBarArray = 0; iBarArray < BarArrays.Count; iBarArray++)
                {
                    uint id0 = BarArrays[(int)iBarArray].Id;
                    if (maxId < id0)
                    {
                        maxId = id0;
                    }
                    System.Diagnostics.Debug.Assert(TypeLocs.Count > id0);
                    System.Diagnostics.Debug.Assert(TypeLocs[(int)id0].Type == 1); // BAR
                    int loc0 = TypeLocs[(int)id0].Loc;
                    System.Diagnostics.Debug.Assert(loc0 == (int)iBarArray);
                    System.Diagnostics.Debug.Assert(BarArrays.Count > loc0);
                    MeshBarArray3D bar0 = BarArrays[(int)loc0];
                    System.Diagnostics.Debug.Assert(bar0.Id == id0);
                }
                for (uint iTriArray = 0; iTriArray < TriArrays.Count; iTriArray++)
                {
                    uint id0 = TriArrays[(int)iTriArray].Id;
                    if (maxId < id0)
                    {
                        maxId = id0;
                    }
                    System.Diagnostics.Debug.Assert(TypeLocs.Count > id0);
                    System.Diagnostics.Debug.Assert(TypeLocs[(int)id0].Type == 2); // TRI
                    int loc0 = TypeLocs[(int)id0].Loc;
                    System.Diagnostics.Debug.Assert(loc0 == (int)iTriArray);
                    System.Diagnostics.Debug.Assert(TriArrays.Count > loc0);
                    MeshTriArray3D tri0 = TriArrays[(int)loc0];
                    System.Diagnostics.Debug.Assert(tri0.Id == id0);
                }
                System.Diagnostics.Debug.Assert(maxId == FindMaxId());
                uint[] isUsedFlgs = new uint[maxId + 1];
                for (uint iVer = 0; iVer < Vertexs.Count; iVer++)
                {
                    System.Diagnostics.Debug.Assert(isUsedFlgs[(int)Vertexs[(int)iVer].Id] == 0);
                    isUsedFlgs[(int)Vertexs[(int)iVer].Id] = 1;
                }
                for (uint iBarArray = 0; iBarArray < BarArrays.Count; iBarArray++)
                {
                    System.Diagnostics.Debug.Assert(isUsedFlgs[(int)BarArrays[(int)iBarArray].Id] == 0);
                    isUsedFlgs[(int)BarArrays[(int)iBarArray].Id] = 1;
                }
                for (uint iTriArray = 0; iTriArray < TriArrays.Count; iTriArray++)
                {
                    System.Diagnostics.Debug.Assert(isUsedFlgs[(int)TriArrays[(int)iTriArray].Id] == 0);
                    isUsedFlgs[(int)TriArrays[(int)iTriArray].Id] = 1;
                }
                System.Diagnostics.Debug.Assert(isUsedFlgs[0] == 0);
            }
            for (uint iVer = 0; iVer < Vertexs.Count; iVer++)
            {
                System.Diagnostics.Debug.Assert(IsId(Vertexs[(int)iVer].Id));
            }
            for (uint iBarArray = 0; iBarArray < BarArrays.Count; iBarArray++)
            {
                System.Diagnostics.Debug.Assert(IsId(BarArrays[(int)iBarArray].Id));
            }
            for (uint iTriArray = 0; iTriArray < TriArrays.Count; iTriArray++)
            {
                System.Diagnostics.Debug.Assert(IsId(TriArrays[(int)iTriArray].Id));
            }
            for (uint index = 0; index < TypeLocs.Count; index++)
            {
                if (TypeLocs[(int)index].Loc == -1)
                {
                    continue;
                }
                System.Diagnostics.Debug.Assert(TypeLocs[(int)index].Loc >= 0);
                System.Diagnostics.Debug.Assert(IsId(index));
            }
            return 0;
        }

        private bool MakeMeshElemLength(
            Cad3D cad, IList<uint> loopIds, IList<double> loopELens,
            IList<uint> edgeIds, IList<double> edgeELens)
        {
            System.Diagnostics.Debug.Assert(loopIds.Count == loopELens.Count);
            System.Diagnostics.Debug.Assert(edgeIds.Count == edgeELens.Count);

            ClearMeshData();

            {
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
                // フラグを立てた頂点のメッシュオブジェクトを生成する
                for (uint vId = 0; vId < vtxFlgs.Count; vId++)
                {
                    if (vtxFlgs[(int)vId] == 0)
                    {
                        continue;
                    }
                    uint addId = GetFreeObjectId();
                    OpenTK.Vector3d vec = cad.GetVertexCoord(vId);
                    Vecs.Add(vec);
                    {
                        MeshVertex3D tmpVer = new MeshVertex3D();
                        tmpVer.Id = addId;
                        tmpVer.VCadId = vId;
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
            }

            // 辺を作る
            for (int iEId = 0; iEId < edgeIds.Count; iEId++)
            {
                uint eId = edgeIds[iEId];
                double eLen = edgeELens[iEId];
                if (GetIdFromCadId(eId, CadElementType.Edge) != 0)
                {
                    // 既にこの辺はMeshに存在
                    continue;
                }
                MakeMeshEdge(cad, eId, eLen);
                System.Diagnostics.Debug.Assert(CheckMesh() == 0);
            }

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
                            // 既にこの辺はMeshに存在
                            continue;
                        }
                        MakeMeshEdge(cad, eId, eLen);
                        System.Diagnostics.Debug.Assert(CheckMesh() == 0);
                    }
                }
            }

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

        private bool MakeMeshEdge(Cad3D cad, uint eId, double len)
        {
            System.Diagnostics.Debug.Assert(GetIdFromCadId(eId, CadElementType.Edge) == 0);

            uint sVId;
            uint eVId;
            if (!cad.GetEdgeVertexId(eId, out sVId, out eVId))
            {
                System.Diagnostics.Debug.Assert(false);
            }

            // 始点、終点のメッシュ点番号をsPId,ePIdに代入
            uint sPId;
            uint ePId;
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
                MeshVertex3D sVP = Vertexs[(int)loc];
                sPId = sVP.V;
                sMeshId = sVP.Id;
                if (!FindElemLocTypeFromCadIdType(CadElementType.Vertex, eVId, out loc, out type))
                {
                    System.Diagnostics.Debug.Assert(false);
                }
                System.Diagnostics.Debug.Assert(type == 0 && loc < Vertexs.Count);
                MeshVertex3D eVP = Vertexs[(int)loc];
                ePId = eVP.V;
                eMeshId = eVP.Id;
            }

            uint newElemId = GetFreeObjectId();
            uint iBarArray0 = (uint)BarArrays.Count;
            BarArrays.Add(new MeshBarArray3D());
            {
                int typeLocCnt = TypeLocs.Count;
                for (int i = typeLocCnt; i < newElemId + 1; i++)
                {
                    var typeLoc = new MeshTypeLoc { Type = 0, Loc = -1 };
                    TypeLocs.Add(typeLoc);
                }
                TypeLocs[(int)newElemId].Loc = (int)iBarArray0;
                TypeLocs[(int)newElemId].Type = 1; // BAR
            }
            MeshBarArray3D barArray = BarArrays[(int)iBarArray0];
            IList<OpenTK.Vector3d> pts;
            cad.GetCurveAsPolyline(eId, out pts, len);
            ////////////////
            uint div = (uint)pts.Count + 1;
            IList<uint> ptIds = new List<uint>();
            {
                for (int i = 0; i < div + 1; i++)
                {
                    ptIds.Add(0);
                }
                ptIds[0] = sPId;
                for (int i = 1; i < div; i++)
                {
                    ptIds[i] = (uint)Vecs.Count;
                    Vecs.Add(pts[i - 1]);
                }
                ptIds[(int)div] = ePId;
            }
            {
                barArray.Id = newElemId;
                barArray.ECadId = eId;
                barArray.Bars.Clear();
                for (int ibar = 0; ibar < div; ibar++)
                {
                    MeshBar bar = new MeshBar();
                    bar.V[0] = ptIds[ibar];
                    bar.V[1] = ptIds[ibar + 1];
                    bar.S2[0] = 0;
                    bar.S2[1] = 0;
                    bar.R2[0] = 0;
                    bar.R2[1] = 0;
                    barArray.Bars.Add(bar);
                }
            }
            System.Diagnostics.Debug.Assert(CheckMesh() == 0);
            return true;
        }

        private bool MakeMeshLoop(Cad3D cad, uint lId, double len)
        {
            Cad3DToXY cad2D = new Cad3DToXY(cad, CadElementType.Loop, lId);

            IList<OpenTK.Vector2d> initialVec2Ds;
            Dictionary<uint, uint> initialVecId3D2D;
            Dictionary<uint, uint> initialVecId2D3D;
            IList<MeshVertex2D> initialVertex2Ds;
            IList<MeshBarArray2D> initialBarArray2Ds;
            ConvertDoneMeshToXY(
                cad, lId,
                cad2D,
                out initialVec2Ds,
                out initialVecId3D2D,
                out initialVecId2D3D,
                out initialVertex2Ds,
                out initialBarArray2Ds);

            Mesher2D mesh2D = new Mesher2D(
                cad2D, len, initialVec2Ds, initialVertex2Ds, initialBarArray2Ds);
            
            // VEC
            var vec2Ds = mesh2D.GetVectors();
            Dictionary<uint, uint> vecId2D3D = new Dictionary<uint, uint>(initialVecId2D3D);

            // VERTEX
            var vertex2Ds = mesh2D.GetVertexs();
            foreach (MeshVertex2D vertex2D in vertex2Ds)
            {
                uint vId2D = vertex2D.VCadId;
                uint vId = cad2D.GetVertexId3DFrom2D(vId2D);
                if (GetIdFromCadId(vId, CadElementType.Vertex) != 0)
                {
                    continue;
                }
                uint addId = GetFreeObjectId();
                OpenTK.Vector3d vec = cad.GetVertexCoord(vId);
                Vecs.Add(vec);
                {
                    MeshVertex3D tmpVer = new MeshVertex3D();
                    tmpVer.Id = addId;
                    tmpVer.VCadId = vId;
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

            // EDGE
            var barArray2Ds = mesh2D.GetBarArrays();
            foreach (MeshBarArray2D barArray2D in barArray2Ds)
            {
                uint eId2D = barArray2D.ECadId;
                uint eId = cad2D.GetEdgeId3DFrom2D(eId2D);
                if (GetIdFromCadId(eId, CadElementType.Edge) != 0)
                {
                    // 登録済み
                    continue;
                }

                IList<MeshBar> tmpBars = new List<MeshBar>();

                int elementCnt = barArray2D.Bars.Count;
                for (int iBar = 0; iBar < elementCnt; iBar++)
                {
                    MeshBar tmpBar = new MeshBar();

                    MeshBar bar2D = barArray2D.Bars[iBar];
                    for (int i = 0; i < 2; i++)
                    {
                        uint ivec2D = bar2D.V[i];
                        uint ivec;
                        if (vecId2D3D.ContainsKey(ivec2D))
                        {
                            ivec = vecId2D3D[ivec2D];
                        }
                        else
                        {
                            OpenTK.Vector2d vec2D = vec2Ds[(int)ivec2D];
                            OpenTK.Vector3d vec = cad2D.UnProject(vec2D);
                            Vecs.Add(vec);
                            ivec = (uint)(Vecs.Count - 1);
                            vecId2D3D[ivec2D] = ivec;
                        }
                        tmpBar.V[i] = ivec;
                    }
                    tmpBars.Add(tmpBar);
                }

                uint newElemId = GetFreeObjectId();
                uint iBarArray = (uint)BarArrays.Count;
                BarArrays.Add(new MeshBarArray3D());
                BarArrays[(int)iBarArray].Bars = tmpBars;
                BarArrays[(int)iBarArray].ECadId = eId;
                BarArrays[(int)iBarArray].Id = newElemId;

                int typeLocCnt = TypeLocs.Count;
                for (int i = typeLocCnt; i < newElemId + 1; i++)
                {
                    var typeLoc = new MeshTypeLoc { Type = 0, Loc = -1 };
                    TypeLocs.Add(typeLoc);
                }
                TypeLocs[(int)newElemId].Loc = (int)iBarArray;
                TypeLocs[(int)newElemId].Type = 1; // BAR
            }

            // LOOP
            System.Diagnostics.Debug.Assert(GetIdFromCadId(lId, CadElementType.Loop) == 0); // 未登録でなければならない
            var triArray2Ds = mesh2D.GetTriArrays();
            foreach (MeshTriArray2D triArray2D in triArray2Ds)
            {
                IList<MeshTri3D> tmpTris = new List<MeshTri3D>();

                int elementCnt = triArray2D.Tris.Count;
                for (int iTri = 0; iTri < elementCnt; iTri++)
                {
                    MeshTri3D tmpTri = new MeshTri3D();

                    MeshTri2D tri2D = triArray2D.Tris[iTri];
                    for (int iTriNo = 0; iTriNo < 3; iTriNo++)
                    {
                        uint ivec2D = tri2D.V[iTriNo];
                        uint ivec;
                        if (vecId2D3D.ContainsKey(ivec2D))
                        {
                            ivec = vecId2D3D[ivec2D];
                        }
                        else
                        {
                            OpenTK.Vector2d vec2D = vec2Ds[(int)ivec2D];
                            OpenTK.Vector3d vec = cad2D.UnProject(vec2D);
                            Vecs.Add(vec);
                            ivec = (uint)(Vecs.Count - 1);
                            vecId2D3D[ivec2D] = ivec;
                        }
                        tmpTri.V[iTriNo] = ivec;
                    }
                    tmpTris.Add(tmpTri);
                }

                uint newElemId = GetFreeObjectId();
                uint iTriArray = (uint)TriArrays.Count;
                TriArrays.Add(new MeshTriArray3D());
                TriArrays[(int)iTriArray].Tris = tmpTris;
                TriArrays[(int)iTriArray].LCadId = lId;
                TriArrays[(int)iTriArray].Id = newElemId;

                int typeLocCnt = TypeLocs.Count;
                for (int i = typeLocCnt; i < newElemId + 1; i++)
                {
                    var typeLoc = new MeshTypeLoc { Type = 0, Loc = -1 };
                    TypeLocs.Add(typeLoc);
                }
                TypeLocs[(int)newElemId].Type = 2;   // TRI
                TypeLocs[(int)newElemId].Loc = (int)iTriArray;

                System.Diagnostics.Debug.Assert(CheckMesh() == 0);
            }
            return true;
        }
    }
}
