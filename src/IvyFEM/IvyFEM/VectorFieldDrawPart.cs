using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK.Graphics.OpenGL;

namespace IvyFEM
{
    public class VectorFieldDrawPart
    {
        public uint MeshId { get; set; } = 0;
        public ElementType Type { get; private set; } = ElementType.NotSet;
        public int Layer { get; set; } = 0;
        public double[] Color { get; set; } = new double[3] { 0.0, 0.0, 0.0 };
        public uint ElemCount { get; set; } = 0;
        public uint ValueDof { get; private set; } = 0;
        public double[] Coords { get; set; } = null;
        public double[] Values { get; set; } = null;
        public VectorFieldDrawerType DrawerType { get; private set; } = VectorFieldDrawerType.NotSet;

        public uint Dimension { get; set; } = 0; // Worldからとってくる

        public VectorFieldDrawPart()
        {

        }

        public VectorFieldDrawPart(VectorFieldDrawPart src)
        {
            MeshId = src.MeshId;
            Type = src.Type;
            Layer = src.Layer;
            src.Color.CopyTo(Color, 0);
            ElemCount = src.ElemCount;
            ValueDof = src.ValueDof;
            Coords = null;
            if (src.Coords != null)
            {
                Coords = new double[src.Coords.Length];
                src.Coords.CopyTo(Coords, 0);
            }
            Values = null;
            if (src.Values != null)
            {
                Values = new double[src.Values.Length];
                src.Values.CopyTo(Values, 0);
            }
            DrawerType = src.DrawerType;
        }


        public VectorFieldDrawPart(uint meshId, FEWorld world)
        {
            var mesh = world.Mesh;
            if (!mesh.IsId(meshId))
            {
                return;
            }
            MeshId = meshId;
            //!!!!!!!!!!!
            Dimension = world.Dimension;

            uint cadId;
            int layer;
            uint elemCount;
            MeshType meshType;
            int loc;
            mesh.GetInfo(MeshId, out cadId, out layer);
            mesh.GetMeshInfo(MeshId, out elemCount, out meshType, out loc, out cadId);
            Layer = layer;
            ElemCount = elemCount;

            if (meshType == MeshType.Vertex)
            {
                Type = ElementType.Point;
            }
            else if (meshType == MeshType.Bar)
            {
                Type = ElementType.Line;
            }
            else if (meshType == MeshType.Tri)
            {
                Type = ElementType.Tri;
                if (Dimension == 2)
                {
                    SetTri(world);
                }
            }
            else if (meshType == MeshType.Tet)
            {
                Type = ElementType.Tet;
                if (Dimension == 3)
                {
                    SetTet(world);
                }
            }
            else
            {
                throw new NotImplementedException();
            }
        }

        private void SetTri(FEWorld world)
        {
            System.Diagnostics.Debug.Assert(Type == ElementType.Tri);
            if (Type != ElementType.Tri)
            {
                return;
            }
            var mesh = world.Mesh;

            int elemPtCount = 3;
            MeshType meshType;
            int[] vertexs;
            mesh.GetConnectivity(MeshId, out meshType, out vertexs);
            System.Diagnostics.Debug.Assert(elemPtCount * ElemCount == vertexs.Length);

            uint dim = Dimension;
            System.Diagnostics.Debug.Assert(dim == 2);
            Coords = new double[ElemCount * dim];
            for (int iTri = 0; iTri < ElemCount; iTri++)
            {
                double[] bubbleCoord = new double[dim];
                for (int iPt = 0; iPt < elemPtCount; iPt++)
                {
                    int coId = vertexs[iTri * elemPtCount + iPt];
                    double[] coord = world.GetVertexCoord(coId);
                    for (int iDimTmp = 0; iDimTmp < dim; iDimTmp++)
                    {
                        bubbleCoord[iDimTmp] += coord[iDimTmp];
                    }
                }
                for (int iDim = 0; iDim < dim; iDim++)
                {
                    bubbleCoord[iDim] /= elemPtCount;

                    Coords[iTri * dim + iDim] = bubbleCoord[iDim];
                }
            }
        }

        private void SetTet(FEWorld world)
        {
            System.Diagnostics.Debug.Assert(Type == ElementType.Tet);
            if (Type != ElementType.Tet)
            {
                return;
            }
            var mesh = world.Mesh;

            int elemPtCount = 4;
            MeshType meshType;
            int[] vertexs;
            mesh.GetConnectivity(MeshId, out meshType, out vertexs);
            System.Diagnostics.Debug.Assert(elemPtCount * ElemCount == vertexs.Length);

            uint dim = Dimension;
            System.Diagnostics.Debug.Assert(dim == 3);
            Coords = new double[ElemCount * dim];
            for (int iTri = 0; iTri < ElemCount; iTri++)
            {
                double[] bubbleCoord = new double[dim];
                for (int iPt = 0; iPt < elemPtCount; iPt++)
                {
                    int coId = vertexs[iTri * elemPtCount + iPt];
                    double[] coord = world.GetVertexCoord(coId);
                    for (int iDimTmp = 0; iDimTmp < dim; iDimTmp++)
                    {
                        bubbleCoord[iDimTmp] += coord[iDimTmp];
                    }
                }
                for (int iDim = 0; iDim < dim; iDim++)
                {
                    bubbleCoord[iDim] /= elemPtCount;

                    Coords[iTri * dim + iDim] = bubbleCoord[iDim];
                }
            }
        }

        public void Update(uint valueId, FieldDerivativeType dt, VectorFieldDrawerType drawerType, FEWorld world)
        {
            DrawerType = drawerType;

            if (DrawerType == VectorFieldDrawerType.Vector)
            {
                UpdateVector(valueId, dt, world);
            }
            else if (DrawerType == VectorFieldDrawerType.SymmetricTensor2)
            {
                UpdateSymmetricTensor2(valueId, dt, world);
            }
            else
            {
                throw new NotImplementedException();
            }
        }

        private void UpdateVector(uint valueId, FieldDerivativeType dt, FEWorld world)
        {
            //ValueDof = 2; // for 2D

            FieldValue fv = world.GetFieldValue(valueId);
            uint quantityId = fv.QuantityId;
            uint dof = fv.Dof;
            System.Diagnostics.Debug.Assert(fv.IsBubble == true);
            var mesh = world.Mesh;
            MeshType meshType;
            int[] vertexs;
            mesh.GetConnectivity(MeshId, out meshType, out vertexs);

            if (Dimension == 2 && Type == ElementType.Tri)
            {
                ValueDof = 2;
                Values = new double[ElemCount * ValueDof];
                for (int iTri = 0; iTri < ElemCount; iTri++)
                {
                    // Bubble
                    uint feId = world.GetTriangleFEIdFromMesh(quantityId, MeshId, (uint)iTri);
                    System.Diagnostics.Debug.Assert(feId != 0);
                    System.Diagnostics.Debug.Assert(dof >= ValueDof);
                    for (int iDof = 0; iDof < ValueDof; iDof++)
                    {
                        double u = fv.GetShowValue((int)(feId - 1), iDof, dt);
                        Values[iTri * ValueDof + iDof] = u;
                    }
                }
            }
            else if (Dimension == 3 && Type == ElementType.Tet)
            {
                ValueDof = 3;
                Values = new double[ElemCount * ValueDof];
                for (int iTet = 0; iTet < ElemCount; iTet++)
                {
                    // Bubble
                    uint feId = world.GetTetrahedronFEIdFromMesh(quantityId, MeshId, (uint)iTet);
                    System.Diagnostics.Debug.Assert(feId != 0);
                    System.Diagnostics.Debug.Assert(dof == ValueDof);
                    for (int iDof = 0; iDof < ValueDof; iDof++)
                    {
                        double u = fv.GetShowValue((int)(feId - 1), iDof, dt);
                        Values[iTet * ValueDof + iDof] = u;
                    }
                }
            }
        }

        private void UpdateSymmetricTensor2(uint valueId, FieldDerivativeType dt, FEWorld world)
        {
            ValueDof = 6;

            FieldValue fv = world.GetFieldValue(valueId);
            uint quantityId = fv.QuantityId;
            uint dof = fv.Dof;
            System.Diagnostics.Debug.Assert(fv.IsBubble == true);
            var mesh = world.Mesh;
            MeshType meshType;
            int[] vertexs;
            mesh.GetConnectivity(MeshId, out meshType, out vertexs);

            if (Type == ElementType.Tri)
            {
                Values = new double[ElemCount * ValueDof];
                for (int iTri = 0; iTri < ElemCount; iTri++)
                {
                    // Bubble
                    uint feId = world.GetTriangleFEIdFromMesh(quantityId, MeshId, (uint)iTri);
                    System.Diagnostics.Debug.Assert(feId != 0);
                    double[] sigma = new double[dof];
                    for (int iDof = 0; iDof < dof; iDof++)
                    {
                        sigma[iDof] = fv.GetShowValue((int)(feId - 1), iDof, dt);
                    }

                    double[] vecs;
                    double ls;
                    double[] vecl;
                    double ll;
                    GetPrincipalStressVectorForSymmetricTensor2(sigma,
                        out vecs, out ls,
                        out vecl, out ll);
                    Values[iTri * ValueDof + 0] = vecs[0];
                    Values[iTri * ValueDof + 1] = vecs[1];
                    Values[iTri * ValueDof + 2] = ls;
                    Values[iTri * ValueDof + 3] = vecl[0];
                    Values[iTri * ValueDof + 4] = vecl[1];
                    Values[iTri * ValueDof + 5] = ll;
                }
            }
        }

        // 主応力
        private void GetPrincipalStressVectorForSymmetricTensor2(double[] sigma, 
            out double[] vecs, out double ls,
            out double[] vecl, out double ll)
        {
            vecs = new double[2];
            ls = 0;
            vecl = new double[2];
            ll = 0;
            {
                double tmp1 = Math.Sqrt((sigma[0] - sigma[1]) * (sigma[0] - sigma[1]) + 4 * sigma[2] * sigma[2]);
                double tmp2 = sigma[0] + sigma[1];
                double l1 = 0.5 * (tmp2 - tmp1);
                double l2 = 0.5 * (tmp2 + tmp1);
                if (Math.Abs(l1) > Math.Abs(l2))
                {
                    ll = l1;
                    ls = l2;
                }
                else
                {
                    ll = l2;
                    ls = l1;
                }
            }
            {
                double[] a1 = { -sigma[2], sigma[0] - ls };
                double[] a2 = { sigma[1] - ls, -sigma[2] };
                double sqlen1 = a1[0] * a1[0] + a1[1] * a1[1];
                double sqlen2 = a2[0] * a2[0] + a2[1] * a2[1];
                if (sqlen1 > sqlen2)
                {
                    vecs[0] = a1[0];
                    vecs[1] = a1[1];
                }
                else
                {
                    vecs[0] = a2[0];
                    vecs[1] = a2[1];
                }
                double len = Math.Sqrt(vecs[0] * vecs[0] + vecs[1] * vecs[1]);
                if (len < 1.0e-10)
                {
                    vecs[0] = 0;
                    vecs[1] = 0;
                }
                else
                {
                    double normalizer = ls / len;
                    vecs[0] *= normalizer;
                    vecs[1] *= normalizer;
                }
            }
            {
                double[] a1 = { -sigma[2], sigma[0] - ll };
                double[] a2 = { sigma[1] - ll, -sigma[2] };
                double sqlen1 = a1[0] * a1[0] + a1[1] * a1[1];
                double sqlen2 = a2[0] * a2[0] + a2[1] * a2[1];
                if (sqlen1 > sqlen2)
                {
                    vecl[0] = a1[0];
                    vecl[1] = a1[1];
                }
                else
                {
                    vecl[0] = a2[0];
                    vecl[1] = a2[1];
                }
                double len = Math.Sqrt(vecl[0] * vecl[0] + vecl[1] * vecl[1]);
                if (len < 1.0e-10)
                {
                    vecl[0] = 0;
                    vecl[1] = 0;
                }
                else
                {
                    double normalizer = ll / len;
                    vecl[0] *= normalizer;
                    vecl[1] *= normalizer;
                }
            }
        }

        public void DrawElements()
        {
            uint dim = Dimension;
            if (dim == 2)
            {
                if (Type != ElementType.Tri)
                {
                    return;
                }

                if (DrawerType == VectorFieldDrawerType.Vector)
                {
                    System.Diagnostics.Debug.Assert(ValueDof == 2);
                    for (int iElem = 0; iElem < ElemCount; iElem++)
                    {
                        double[] co = { Coords[iElem * dim], Coords[iElem * dim + 1] };
                        double[] va = new double[ValueDof];
                        for (int iDof = 0; iDof < ValueDof; iDof++)
                        {
                            va[iDof] = Values[iElem * ValueDof + iDof];
                        }

                        GL.Color3(Color);
                        GL.LineWidth(1);
                        GL.Begin(PrimitiveType.Lines);
                        GL.Vertex2(co);
                        GL.Vertex2(co[0] + va[0], co[1] + va[1]);
                        GL.End();

                        double vaLen = Math.Sqrt(va[0] * va[0] + va[1] * va[1]);
                        double arrowTheta = Math.PI / 12.0;
                        double arrowLen = vaLen * 0.2;
                        if (Math.Sqrt(va[0] * va[0] + va[1] * va[1]) >= arrowLen *Math.Cos(arrowTheta) * 1.5)
                        {
                            // arrow
                            double vecTheta = Math.Atan2(va[1], va[0]);
                            double[] cPt = new double[]
                            {
                                co[0] + va[0] - arrowLen * Math.Cos(arrowTheta) * Math.Cos(vecTheta),
                                co[1] + va[1] - arrowLen * Math.Cos(arrowTheta) * Math.Sin(vecTheta)
                            };
                            double[] pt1 = new double[]
                            {
                                cPt[0] - arrowLen * Math.Sin(arrowTheta) * Math.Sin(vecTheta),
                                cPt[1] + arrowLen * Math.Sin(arrowTheta) * Math.Cos(vecTheta)
                            };
                            double[] pt2 = new double[]
                            {
                                cPt[0] + arrowLen * Math.Sin(arrowTheta) * Math.Sin(vecTheta),
                                cPt[1] - arrowLen * Math.Sin(arrowTheta) * Math.Cos(vecTheta)
                            };
                            GL.Begin(PrimitiveType.Lines);
                            GL.Vertex2(co[0] + va[0], co[1] + va[1]);
                            GL.Vertex2(pt1[0], pt1[1]);
                            GL.Vertex2(co[0] + va[0], co[1] + va[1]);
                            GL.Vertex2(pt2[0], pt2[1]);
                            GL.End();
                        }
                    }
                }
                else if (DrawerType == VectorFieldDrawerType.SymmetricTensor2)
                {
                    System.Diagnostics.Debug.Assert(ValueDof == 6);
                    for (int iElem = 0; iElem < ElemCount; iElem++)
                    {
                        double[] co = { Coords[iElem * dim], Coords[iElem * dim + 1] };
                        double[] va = new double[ValueDof];
                        for (int iDof = 0; iDof < ValueDof; iDof++)
                        {
                            va[iDof] = Values[iElem * ValueDof + iDof];
                        }

                        if (va[2] > 0)
                        {
                            GL.Color3(0.0, 0.0, 1.0);
                        }
                        else
                        {
                            GL.Color3(1.0, 0.0, 0.0);
                        }
                        GL.LineWidth(1);
                        GL.Begin(PrimitiveType.Lines);
                        GL.Vertex2(co);
                        GL.Vertex2(co[0] + va[0], co[1] + va[1]);

                        GL.Vertex2(co);
                        GL.Vertex2(co[0] - va[0], co[1] - va[1]);
                        GL.End();

                        if (va[5] > 0)
                        {
                            GL.Color3(0.0, 0.0, 1.0);
                        }
                        else
                        {
                            GL.Color3(1.0, 0.0, 0.0);
                        }
                        GL.Begin(PrimitiveType.Lines);
                        GL.Vertex2(co);
                        GL.Vertex2(co[0] + va[3], co[1] + va[4]);

                        GL.Vertex2(co);
                        GL.Vertex2(co[0] - va[3], co[1] - va[4]);
                        GL.End();
                    }
                }
            }
            else if (dim == 3)
            {
                if (Type != ElementType.Tet)
                {
                    return;
                }
                if (DrawerType == VectorFieldDrawerType.Vector)
                {
                    System.Diagnostics.Debug.Assert(ValueDof == 3);
                    for (int iElem = 0; iElem < ElemCount; iElem++)
                    {
                        double[] co = { Coords[iElem * dim], Coords[iElem * dim + 1], Coords[iElem * dim + 2] };
                        double[] va = new double[ValueDof];
                        for (int iDof = 0; iDof < ValueDof; iDof++)
                        {
                            va[iDof] = Values[iElem * ValueDof + iDof];
                        }

                        /*
                        GL.LineWidth(1);
                        GL.Begin(PrimitiveType.Lines);
                        GL.Color3(Color);
                        GL.Vertex3(co);
                        GL.Vertex3(co[0] + va[0], co[1] + va[1], co[2] + va[2]);
                        GL.End();

                        {
                            OpenTK.Vector3d p1 = new OpenTK.Vector3d(co[0], co[1], co[2]);
                            OpenTK.Vector3d p2 = new OpenTK.Vector3d(co[0] + va[0], co[1] + va[1], co[2] + va[2]);
                            OpenTK.Vector3d p3 = p1 + 0.8 * (p2 - p1);

                            // arrowの代わり
                            GL.LineWidth(4);
                            GL.Begin(PrimitiveType.Lines);
                            GL.Vertex3(p3);
                            GL.Vertex3(p2);
                            GL.End();
                        }
                        */
                        DrawArrow3D(
                            new OpenTK.Vector3d(co[0], co[1], co[2]),
                            new OpenTK.Vector3d(co[0] + va[0], co[1] + va[1], co[2] + va[2]));
                    }
                }
            }
        }

        private void rightAngleVector(OpenTK.Vector3d src3, out OpenTK.Vector3d dst13, out OpenTK.Vector3d dst23)
        {
            var tmp = new OpenTK.Vector3d(1, 0, 0);

            //tmpとsrc3の角度０またはそれに限りなく近いなら、別なベクトルを用意

            if (OpenTK.Vector3d.CalculateAngle(tmp, src3) < 0.1 * Math.PI / 180.0)
            {
                tmp.X = 0;
                tmp.Y = 1;
                tmp.Z = 0;
            }
            //外積を求める
            dst13 = OpenTK.Vector3d.Cross(tmp, src3);
            dst23 = OpenTK.Vector3d.Cross(src3, dst13);
        }

        private void DrawArrow3D(OpenTK.Vector3d from, OpenTK.Vector3d to)
        {
            GL.LineWidth(1);
            GL.Begin(PrimitiveType.Lines);
            GL.Color3(Color);
            GL.Vertex3(from);
            GL.Vertex3(to);
            GL.End();

            var v = to - from;
            double len = v.Length;

            v *= 0.8;

            OpenTK.Vector3d v1;
            OpenTK.Vector3d v2;
            rightAngleVector(v, out v1, out v2);

            v1.Normalize();
            v2.Normalize();
            v1 *= len * 0.03;
            v2 *= len * 0.03;

            var f = from;
            var t = to;

            GL.Begin(PrimitiveType.LineStrip);
            GL.Vertex3(f + v);
            GL.Vertex3(f + v + v1);
            GL.Vertex3(t);
            GL.End();
            GL.Begin(PrimitiveType.LineStrip);
            GL.Vertex3(f + v);
            GL.Vertex3(f + v + v2);
            GL.Vertex3(t);
            GL.End();
            GL.Begin(PrimitiveType.LineStrip);
            GL.Vertex3(f + v);
            GL.Vertex3(f + v - v1);
            GL.Vertex3(t);
            GL.End();
            GL.Begin(PrimitiveType.LineStrip);
            GL.Vertex3(f + v);
            GL.Vertex3(f + v - v2);
            GL.Vertex3(t);
            GL.End();
        }
    }
}
