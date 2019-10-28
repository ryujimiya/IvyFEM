using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshUtils
    {
        /// <summary>
        /// 線分にいくら点があるか
        /// </summary>
        public const uint EdNo = 2;

        /// <summary>
        /// 三角形にいくら頂点があるか
        /// </summary>
        public const uint TriNo = 3;
        /// <summary>
        /// 三角形にいくら辺があるか
        /// </summary>
        public const uint TriEdNo = 3;
        /// <summary>
        /// 三角形の各辺の頂点番号
        /// </summary>
        public static uint[][] TriElEdgeNo = new uint[(int)TriEdNo][]
        {
            new uint[(int)EdNo]{ 1, 2 },
            new uint[(int)EdNo]{ 2, 0 },
            new uint[(int)EdNo]{ 0, 1 }
        };
        /// <summary>
        /// 三角形の隣接関係
        /// </summary>
        public static uint[][] RelTriTri = new uint[3][]
        {
            new uint[3]{ 0, 2, 1 }, //  0
            new uint[3]{ 2, 1, 0 }, //  1 
	        new uint[3]{ 1, 0, 2 } //  2
        };

        /// <summary>
        /// 四角形にいくら頂点があるか
        /// </summary>
        public const uint QuadNo = 4;
        /// <summary>
        /// 四角形にいくら辺があるか
        /// </summary>
        public const uint QuadEdNo = 4;
        /// <summary>
        /// 四角形の各辺の頂点番号
        /// </summary>
        public static uint[][] QuadElEdgeNo = new uint[(int)QuadEdNo][]
        {
            new uint[(int)EdNo]{ 0, 1 },
            new uint[(int)EdNo]{ 1, 2 },
            new uint[(int)EdNo]{ 2, 3 },
            new uint[(int)EdNo]{ 3, 0 }
        };
        /// <summary>
        /// 四角形の隣接関係
        /// </summary>
        public static uint[][] RelQuadQuad = new uint[(int)QuadNo][]
        {
            new uint[(int)QuadNo]{ 0, 3, 2, 1 }, //  
            new uint[(int)QuadNo]{ 1, 0, 3, 2 }, //  1
            new uint[(int)QuadNo]{ 2, 1, 0, 3 }, //  2
            new uint[(int)QuadNo]{ 3, 2, 1, 0 } //  3
        };


        private static uint[] InvRelTriTri = new uint[3]
        {
            0, 1, 2
        };

        /// <summary>
        /// (0に相当するノード番号)*3+(1に相当するノード番号)  →　関係番号 変換
        /// </summary>
        private static int[] El2RelTriTri = new int[9]
        {
            -1, // 0 00
            -1, // 1 01
            0,  // 2 02
            2, // 3 10
            -1, // 4 11
            -1, // 5 12
            -1, // 6 20
            1,  // 7 21
            -1, // 8 22
        };

        /// <summary>
        /// (こちら側の辺番号)*3+(相手側の辺番号)　→　関係番号 変換
        /// </summary>
        private static uint[] Ed2RelTriTri = new uint[9]
        {
            0,  // 0 00
            2,  // 1 01
            1,  // 2 02
            2,  // 3 10
            1,  // 4 11
            0,  // 5 12
            1,  // 6 20
            0,  // 7 21
            2,  // 8 22
        };

        private static uint[][] IndexRot3 = new uint[3][]
        {
            new uint[3] { 0, 1, 2 },
            new uint[3] { 1, 2, 0 },
            new uint[3] { 2, 0, 1 },
        };

        /// <summary>
        /// ドロネー条件を満たすかどうか調べる
        /// </summary>
        /// <param name="p0"></param>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <returns></returns>
        public static int DetDelaunay(
            OpenTK.Vector2d p0, OpenTK.Vector2d p1, OpenTK.Vector2d p2, OpenTK.Vector2d p3)
        {
            double area = CadUtils.TriArea(p0, p1, p2);
            if (Math.Abs(area) < 1.0e-10)
            {
                return 3;
            }
            double tmpVal = 1.0 / (area * area * 16.0);

            double dtmp0 = CadUtils.SquareLength(p1, p2);
            double dtmp1 = CadUtils.SquareLength(p0, p2);
            double dtmp2 = CadUtils.SquareLength(p0, p1);

            double etmp0 = tmpVal * dtmp0 * (dtmp1 + dtmp2 - dtmp0);
            double etmp1 = tmpVal * dtmp1 * (dtmp0 + dtmp2 - dtmp1);
            double etmp2 = tmpVal * dtmp2 * (dtmp0 + dtmp1 - dtmp2);

            OpenTK.Vector2d outCenter = new OpenTK.Vector2d(
                etmp0 * p0.X + etmp1 * p1.X + etmp2 * p2.X,
                etmp0 * p0.Y + etmp1 * p1.Y + etmp2 * p2.Y);

            double qradius = CadUtils.SquareLength(outCenter, p0);
            double qdistance = CadUtils.SquareLength(outCenter, p3);

            //System.Diagnostics.Debug.Assert(Math.Abs(qradius - CadUtils.SquareLength(outCenter, p1)) < 1.0e-10 * qradius);
            //System.Diagnostics.Debug.Assert(Math.Abs(qradius - CadUtils.SquareLength(outCenter, p2)) < 1.0e-10 * qradius);

            const double tol = 1.0e-20;
            if (qdistance > qradius * (1.0 + tol))
            {
                // 外接円の外
                return 2;
            }
            else
            {
                if (qdistance < qradius * (1.0 - tol))
                {
                    // 外接円の中
                    return 0;
                }
                else
                {
                    // 外接円上
                    return 1;
                }
            }

            throw new InvalidOperationException();
            //return 0;
        }

        public static bool CheckTri(IList<MeshPoint2D> pts, IList<MeshTri2D> tris)
        {
            uint ptCnt = (uint)pts.Count;
            uint triCnt = (uint)tris.Count;

            ////////////////////////////////
            // 要素Indexのチェック

            for (uint itri = 0; itri < triCnt; itri++)
            {
               MeshTri2D tri = tris[(int)itri];
                for (uint iTriNo = 0; iTriNo < TriNo; iTriNo++)
                {
                    System.Diagnostics.Debug.Assert(tri.V[iTriNo] < ptCnt);
                }
                for (uint iTriEd = 0; iTriEd < TriEdNo; iTriEd++)
                {
                    if (tri.G2[iTriEd] == -2 || tri.G2[iTriEd] == -3)
                    {
                        uint iSTri = tri.S2[iTriEd];
                        uint iRel = tri.R2[iTriEd];
                        System.Diagnostics.Debug.Assert(iSTri < triCnt);
                        System.Diagnostics.Debug.Assert(iRel < 3);
                        // check sorounding
                        {
                            uint elDiaNo = RelTriTri[iRel][iTriEd];
                            System.Diagnostics.Debug.Assert(elDiaNo < 3);
                            if (tris[(int)iSTri].S2[elDiaNo] != itri)
                            {
                                System.Diagnostics.Debug.WriteLine(itri + " " + iTriEd);
                            }
                            System.Diagnostics.Debug.Assert(tris[(int)iSTri].S2[elDiaNo] == itri);
                        }
                        // check relation
                        for (uint iEdNo = 0; iEdNo < EdNo; iEdNo++)
                        {
                            uint iElNo = TriElEdgeNo[iTriEd][iEdNo];
                            if (tri.V[iElNo] != tris[(int)iSTri].V[(int)RelTriTri[iRel][iElNo]])
                            {
                                System.Diagnostics.Debug.WriteLine(itri + " " + iTriEd);
                            }
                            System.Diagnostics.Debug.Assert(tri.V[iElNo] ==
                                tris[(int)iSTri].V[(int)RelTriTri[iRel][iElNo]]);
                        }
                    }
                }
                {
                    if (tri.G2[0] == -1 && tri.G2[1] == -1 && tri.G2[2] == -1)
                    {
                        System.Diagnostics.Debug.WriteLine("Isolated Triangle " + itri);
                    }
                }
            }

            ////////////////////////////////
            // 頂点-要素間の一貫性のチェック

            for (uint iPt = 0; iPt < ptCnt; iPt++)
            {
                if (pts[(int)iPt].Elem >= 0)
                {
                    System.Diagnostics.Debug.Assert(pts[(int)iPt].Dir >= 0 && pts[(int)iPt].Dir < 3);
                    int itri0 = pts[(int)iPt].Elem;
                    uint inoel0 = pts[(int)iPt].Dir;
                    if (tris[itri0].V[inoel0] != iPt)
                    {
                        System.Diagnostics.Debug.WriteLine(itri0 + " " + inoel0 + "   " +
                            tris[itri0].V[inoel0] + " " + iPt);
                    }
                    System.Diagnostics.Debug.Assert(tris[itri0].V[inoel0] == iPt);
                }
            }

            ////////////////////////////////
            // Geometryのチェック

            for (uint iTri = 0; iTri < triCnt; iTri++)
            {
                MeshTri2D tri = tris[(int)iTri];
                {
                    double area = CadUtils.TriArea(
                        pts[(int)tri.V[0]].Point,
                        pts[(int)tri.V[1]].Point,
                        pts[(int)tri.V[2]].Point);
                    if (area < 1.0e-10)
                    {
                        System.Diagnostics.Debug.WriteLine("Negative Volume : " + iTri + " " + area);
                    }
                }

                // コメントアウトされてた部分　↓
                /*
                {
                    Vector2 v0 = pts[(int)tri.V[0]].Point;
                    Vector2 v1 = pts[(int)tri.V[1]].Point;
                    Vector2 v2 = pts[(int)tri.V[2]].Point;

                    double area = CadUtils.TriArea(v0, v1, v2);
                    double tmp1 = 0.5 / area;

                    double[] constTerm = new double[3];
                    constTerm[0] = tmp1 * (v1.X * v2.Y - v2.X * v1.Y);
                    constTerm[1] = tmp1 * (v2.X * v0.Y - v0.X * v2.Y);
                    constTerm[2] = tmp1 * (v0.X * v1.Y - v1.X * v0.Y);

                    double[][] dldx = new double[3][];
                    for (int i = 0; i < 3;  i++)
                    {
                        dldx[i] = new double[2];
                    }
                    dldx[0][0] = tmp1 * (v1.Y - v2.Y);
                    dldx[1][0] = tmp1 * (v2.Y - v0.Y);
                    dldx[2][0] = tmp1 * (v0.Y - v1.Y);

                    dldx[0][1] = tmp1 * (v2.X - v1.X);
                    dldx[1][1] = tmp1 * (v0.X - v2.X);
                    dldx[2][1] = tmp1 * (v1.X - v0.X);

                    System.Diagnostics.Debug.Assert(Math.Abs(dldx[0][0] + dldx[1][0] + dldx[2][0]) < 1.0e-15);
                    System.Diagnostics.Debug.Assert(Math.Abs(dldx[0][1] + dldx[1][1] + dldx[2][1]) < 1.0e-15);

                    System.Diagnostics.Debug.Assert(Math.Abs(constTerm[0] + dldx[0][0] * v0.X +
                        dldx[0][1] * v0.Y - 1.0) < 1.0e-10);
                    System.Diagnostics.Debug.Assert(Math.Abs(constTerm[0] + dldx[0][0] * v1.X +
                        dldx[0][1] * v1.Y) < 1.0e-10);
                    System.Diagnostics.Debug.Assert(Math.Abs(constTerm[0] + dldx[0][0] * v2.X +
                        dldx[0][1] * v2.Y) < 1.0e-10);

                    System.Diagnostics.Debug.Assert(Math.Abs(constTerm[1] + dldx[1][0] * v0.X +
                        dldx[1][1] * v0.Y) < 1.0e-10);
                    System.Diagnostics.Debug.Assert(Math.Abs(constTerm[1] + dldx[1][0] * v1.X +
                        dldx[1][1] * v1.Y - 1.0) < 1.0e-10);
                    System.Diagnostics.Debug.Assert(Math.Abs(constTerm[1] + dldx[1][0] * v2.X +
                        dldx[1][1] * v2.Y) < 1.0e-10);

                    System.Diagnostics.Debug.Assert(Math.Abs(constTerm[2] + dldx[2][0] * v0.X +
                        dldx[2][1] * v0.Y) < 1.0e-10);
                    System.Diagnostics.Debug.Assert(Math.Abs(constTerm[2] + dldx[2][0] * v1.X +
                        dldx[2][1] * v1.Y) < 1.0e-10);
                    System.Diagnostics.Debug.Assert(Math.Abs(constTerm[2] + dldx[2][0] * v2.X +
                        dldx[2][1] * v2.Y - 1.0) < 1.0e-10);
                }
                */
            }

            return true;
        }

        public static bool InsertPointElem(uint iInsPt, uint iInsTri,
                              IList<MeshPoint2D> points, IList<MeshTri2D> tris)
        {
            System.Diagnostics.Debug.Assert(iInsTri < tris.Count);
            System.Diagnostics.Debug.Assert(iInsPt < points.Count);

            int iTri0 = (int)iInsTri;
            int iTri1 = tris.Count;
            int iTri2 = tris.Count + 1;

            int triCnt = tris.Count;
            for (int i = triCnt; i < triCnt + 2; i++)
            {
                tris.Add(new MeshTri2D());
            }

            MeshTri2D oldTri = new MeshTri2D(tris[(int)iInsTri]);

            points[(int)iInsPt].Elem = iTri0;
            points[(int)iInsPt].Dir = 0;
            points[(int)oldTri.V[0]].Elem = iTri1;
            points[(int)oldTri.V[0]].Dir = 2;
            points[(int)oldTri.V[1]].Elem = iTri2;
            points[(int)oldTri.V[1]].Dir = 2;
            points[(int)oldTri.V[2]].Elem = iTri0;
            points[(int)oldTri.V[2]].Dir = 2;

            {
                MeshTri2D tri = tris[iTri0];

                tri.V[0] = iInsPt;
                tri.V[1] = oldTri.V[1];
                tri.V[2] = oldTri.V[2];
                tri.G2[0] = oldTri.G2[0];
                tri.G2[1] = -2;
                tri.G2[2] = -2;
                tri.S2[0] = oldTri.S2[0];
                tri.S2[1] = (uint)iTri1;
                tri.S2[2] = (uint)iTri2;

                if (oldTri.G2[0] == -2 || oldTri.G2[0] == -3)
                {
                    System.Diagnostics.Debug.Assert(oldTri.R2[0] < 3);
                    uint[] rel = RelTriTri[oldTri.R2[0]];
                    tri.R2[0] = (uint)El2RelTriTri[rel[0] * 3 + rel[1]];
                    System.Diagnostics.Debug.Assert(tri.R2[0] >= 0 && tri.R2[0] < 3);
                    System.Diagnostics.Debug.Assert(oldTri.S2[0] < tris.Count);
                    tris[(int)oldTri.S2[0]].S2[rel[0]] = (uint)iTri0;
                    tris[(int)oldTri.S2[0]].R2[rel[0]] = InvRelTriTri[tri.R2[0]];
                }
                tri.R2[1] = 0;
                tri.R2[2] = 0;
            }
            {
                MeshTri2D tri = tris[iTri1];

                tri.V[0] = iInsPt;
                tri.V[1] = oldTri.V[2];
                tri.V[2] = oldTri.V[0];
                tri.G2[0] = oldTri.G2[1];
                tri.G2[1] = -2;
                tri.G2[2] = -2;
                tri.S2[0] = oldTri.S2[1];
                tri.S2[1] = (uint)iTri2;
                tri.S2[2] = (uint)iTri0;

                if (oldTri.G2[1] == -2 || oldTri.G2[1] == -3)
                {
                    System.Diagnostics.Debug.Assert(oldTri.R2[1] < 3);
                    uint[] rel = RelTriTri[oldTri.R2[1]];
                    tri.R2[0] = (uint)El2RelTriTri[rel[1] * 3 + rel[2]];
                    System.Diagnostics.Debug.Assert(tri.R2[0] >= 0 && tri.R2[0] < 3);
                    System.Diagnostics.Debug.Assert(oldTri.S2[1] < tris.Count);
                    tris[(int)oldTri.S2[1]].S2[rel[1]] = (uint)iTri1;
                    tris[(int)oldTri.S2[1]].R2[rel[1]] = InvRelTriTri[tri.R2[0]];
                }
                tri.R2[1] = 0;
                tri.R2[2] = 0;
            }
            {
                MeshTri2D tri = tris[iTri2];

                tri.V[0] = iInsPt;
                tri.V[1] = oldTri.V[0];
                tri.V[2] = oldTri.V[1];
                tri.G2[0] = oldTri.G2[2];
                tri.G2[1] = -2;
                tri.G2[2] = -2;
                tri.S2[0] = oldTri.S2[2];
                tri.S2[1] = (uint)iTri0;
                tri.S2[2] = (uint)iTri1;

                if (oldTri.G2[2] == -2 || oldTri.G2[2] == -3)
                {
                    System.Diagnostics.Debug.Assert(oldTri.R2[2] < 3);
                    uint[] rel = RelTriTri[oldTri.R2[2]];
                    tri.R2[0] = (uint)El2RelTriTri[rel[2] * 3 + rel[0]];
                    System.Diagnostics.Debug.Assert(tri.R2[0] >= 0 && tri.R2[0] < 3);
                    System.Diagnostics.Debug.Assert(oldTri.S2[2] < tris.Count);
                    tris[(int)oldTri.S2[2]].S2[rel[2]] = (uint)iTri2;
                    tris[(int)oldTri.S2[2]].R2[rel[2]] = InvRelTriTri[tri.R2[0]];
                }
                tri.R2[1] = 0;
                tri.R2[2] = 0;
            }

            return true;
        }

        public static bool InsertPointElemEdge(uint iInsPt, uint iInsTri, uint iInsEd,
            IList<MeshPoint2D> points, IList<MeshTri2D> tris)
        {
            System.Diagnostics.Debug.Assert(iInsTri < tris.Count);
            System.Diagnostics.Debug.Assert(iInsPt < points.Count);

            if (tris[(int)iInsTri].G2[iInsEd] != -2)
            {
                // 未実装
                System.Diagnostics.Debug.Assert(false);
                new NotImplementedException();
            }

            uint iAdjTri = tris[(int)iInsTri].S2[iInsEd];
            uint iAdjEd = RelTriTri[(int)tris[(int)iInsTri].R2[iInsEd]][iInsEd];
            System.Diagnostics.Debug.Assert(iAdjTri < tris.Count);
            System.Diagnostics.Debug.Assert(iInsEd < 3);

            uint itri0 = iInsTri;
            uint itri1 = iAdjTri;
            uint itri2 = (uint)tris.Count;
            uint itri3 = (uint)(tris.Count + 1);

            int triCnt = tris.Count;
            for (int i = triCnt; i < triCnt + 2; i++)
            {
                tris.Add(new MeshTri2D());
            }

            MeshTri2D old0 = new MeshTri2D(tris[(int)iInsTri]);
            MeshTri2D old1 = new MeshTri2D(tris[(int)iAdjTri]);

            uint ino00 = iInsEd;
            uint ino10 = TriElEdgeNo[iInsEd][0];
            uint ino20 = TriElEdgeNo[iInsEd][1];

            uint ino01 = iAdjEd;
            uint ino11 = TriElEdgeNo[iAdjEd][0];
            uint ino21 = TriElEdgeNo[iAdjEd][1];

            System.Diagnostics.Debug.Assert(old0.V[ino10] == old1.V[ino21]);
            System.Diagnostics.Debug.Assert(old0.V[ino20] == old1.V[ino11]);
            System.Diagnostics.Debug.Assert(old0.S2[ino00] == itri1);
            System.Diagnostics.Debug.Assert(old1.S2[ino01] == itri0);

            points[(int)iInsPt].Elem = (int)itri0;
            points[(int)iInsPt].Dir = 0;
            points[(int)old0.V[ino20]].Elem = (int)itri0;
            points[(int)old0.V[ino20]].Dir = 1;
            points[(int)old0.V[ino00]].Elem = (int)itri1;
            points[(int)old0.V[ino00]].Dir = 1;
            points[(int)old1.V[ino21]].Elem = (int)itri2;
            points[(int)old1.V[ino21]].Dir = 1;
            points[(int)old1.V[ino01]].Elem = (int)itri3;
            points[(int)old1.V[ino01]].Dir = 1;

            {
                MeshTri2D  tri = tris[(int)itri0];
                tri.V[0] = iInsPt;
                tri.V[1] = old0.V[ino20];
                tri.V[2] = old0.V[ino00];
                tri.G2[0] = old0.G2[ino10];
                tri.G2[1] = -2;
                tri.G2[2] = -2;
                tri.S2[0] = old0.S2[ino10];
                tri.S2[1] = itri1;
                tri.S2[2] = itri3;
                if (old0.G2[ino10] == -2 || old0.G2[ino10] == -3)
                {
                    System.Diagnostics.Debug.Assert(old0.R2[ino10] < 3);
                    uint[] rel = RelTriTri[old0.R2[ino10]];
                    tri.R2[0] = (uint)El2RelTriTri[rel[ino10] * 3 + rel[ino20]];
                    System.Diagnostics.Debug.Assert(tri.R2[0] >= 0 && tri.R2[0] < 3);
                    System.Diagnostics.Debug.Assert(old0.S2[ino10] < tris.Count);
                    tris[(int)old0.S2[ino10]].S2[rel[ino10]] = itri0;
                    tris[(int)old0.S2[ino10]].R2[rel[ino10]] = InvRelTriTri[tri.R2[0]];
                }
                tri.R2[1] = 0;
                tri.R2[2] = 0;
            }
            {
                MeshTri2D tri = tris[(int)itri1];
                tri.V[0] = iInsPt;
                tri.V[1] = old0.V[ino00];
                tri.V[2] = old0.V[ino10];
                tri.G2[0] = old0.G2[ino20];
                tri.G2[1] = -2;
                tri.G2[2] = -2;
                tri.S2[0] = old0.S2[ino20];
                tri.S2[1] = itri2;
                tri.S2[2] = itri0;
                if (old0.G2[ino20] == -2 || old0.G2[ino20] == -3)
                {
                    System.Diagnostics.Debug.Assert(old0.R2[ino20] < 3);
                    uint[] rel = RelTriTri[old0.R2[ino20]];
                    tri.R2[0] = (uint)El2RelTriTri[rel[ino20] * 3 + rel[ino00]];
                    System.Diagnostics.Debug.Assert(tri.R2[0] >= 0 && tri.R2[0] < 3);
                    System.Diagnostics.Debug.Assert(old0.S2[ino20] < tris.Count);
                    tris[(int)old0.S2[ino20]].S2[rel[ino20]] = itri1;
                    tris[(int)old0.S2[ino20]].R2[rel[ino20]] = InvRelTriTri[tri.R2[0]];
                }
                tri.R2[1] = 0;
                tri.R2[2] = 0;
            }
            {
                MeshTri2D tri = tris[(int)itri2];
                tri.V[0] = iInsPt;
                tri.V[1] = old1.V[ino21];
                tri.V[2] = old1.V[ino01];
                tri.G2[0] = old1.G2[ino11];
                tri.G2[1] = -2;
                tri.G2[2] = -2;
                tri.S2[0] = old1.S2[ino11];
                tri.S2[1] = itri3;
                tri.S2[2] = itri1;
                if (old1.G2[ino11] == -2 || old0.G2[ino20] == -3)
                {
                    System.Diagnostics.Debug.Assert(old1.R2[ino11] < 3);
                    uint[] rel = RelTriTri[old1.R2[ino11]];
                    tri.R2[0] = (uint)El2RelTriTri[rel[ino11] * 3 + rel[ino21]];
                    System.Diagnostics.Debug.Assert(tri.R2[0] >= 0 && tri.R2[0] < 3);
                    System.Diagnostics.Debug.Assert(old1.S2[ino11] < tris.Count);
                    tris[(int)old1.S2[ino11]].S2[rel[ino11]] = itri2;
                    tris[(int)old1.S2[ino11]].R2[rel[ino11]] = InvRelTriTri[tri.R2[0]];
                }
                tri.R2[1] = 0;
                tri.R2[2] = 0;
            }
            {
                MeshTri2D tri = tris[(int)itri3];
                tri.V[0] = iInsPt;
                tri.V[1] = old1.V[ino01];
                tri.V[2] = old1.V[ino11];
                tri.G2[0] = old1.G2[ino21];
                tri.G2[1] = -2;
                tri.G2[2] = -2;
                tri.S2[0] = old1.S2[ino21];
                tri.S2[1] = itri0;
                tri.S2[2] = itri2;
                if (old1.G2[ino21] == -2 || old1.G2[ino21] == -3)
                {
                    System.Diagnostics.Debug.Assert(old1.R2[ino21] < 3);
                    uint[] rel = RelTriTri[old1.R2[ino21]];
                    tri.R2[0] = (uint)El2RelTriTri[rel[ino21] * 3 + rel[ino01]];
                    System.Diagnostics.Debug.Assert(tri.R2[0] >= 0 && tri.R2[0] < 3);
                    System.Diagnostics.Debug.Assert(old1.S2[ino21] < tris.Count);
                    tris[(int)old1.S2[ino21]].S2[rel[ino21]] = itri3;
                    tris[(int)old1.S2[ino21]].R2[rel[ino21]] = InvRelTriTri[tri.R2[0]];
                }
                tri.R2[1] = 0;
                tri.R2[2] = 0;
            }
            return true;
        }

        public static bool DelaunayAroundPoint(uint iPt0, IList<MeshPoint2D> points, IList<MeshTri2D> tris)
        {
            System.Diagnostics.Debug.Assert(iPt0 < points.Count);
            if (points[(int)iPt0].Elem == -1)
            {
                return true;
            }

            System.Diagnostics.Debug.Assert(points[(int)iPt0].Elem >= 0 &&
                points[(int)iPt0].Elem < tris.Count );
            System.Diagnostics.Debug.Assert(tris[points[(int)iPt0].Elem].V[points[(int)iPt0].Dir] == iPt0);

            uint iTri0 = (uint)points[(int)iPt0].Elem;
            uint iTriNo0 = points[(int)iPt0].Dir;

            uint iCurTri = iTri0;
            uint iCurTriNo = points[(int)iPt0].Dir;
            bool isWall = false;
            while (true)
            {
                System.Diagnostics.Debug.Assert(tris[(int)iCurTri].V[iCurTriNo] == iPt0);

                if (tris[(int)iCurTri].G2[iCurTriNo] == -2)
                {
                    // 向かいの要素を調べる
                    uint iTriDia = tris[(int)iCurTri].S2[iCurTriNo];
                    uint[] diaRel = RelTriTri[tris[(int)iCurTri].R2[iCurTriNo]];
                    uint iTriDiaNo = diaRel[iCurTriNo];
                    System.Diagnostics.Debug.Assert(tris[(int)iTriDia].G2[iTriDiaNo] == -2);
                    System.Diagnostics.Debug.Assert(tris[(int)iTriDia].S2[iTriDiaNo] == iCurTri);
                    uint iDiaPt = tris[(int)iTriDia].V[iTriDiaNo];
                    if (DetDelaunay(
                        points[(int)tris[(int)iCurTri].V[0]].Point,
                        points[(int)tris[(int)iCurTri].V[1]].Point,
                        points[(int)tris[(int)iCurTri].V[2]].Point,
                        points[(int)iDiaPt].Point) == 0)
                    {
                        // Delaunay条件が満たされない場合

                        // 辺を切り替える
                        // FlipEdgeによってiCurTriは時計回り側の３角形に切り替わる
                        FlipEdge(iCurTri, iCurTriNo, points, tris);

                        iCurTriNo = 2;
                        System.Diagnostics.Debug.Assert(tris[(int)iCurTri].V[iCurTriNo] == iPt0);
                        // Flipによってtris[itri0].V[inotri0] != ipo0 でなくなってしまうのを防ぐため
                        if (iCurTri == iTri0)
                        {
                            iTriNo0 = iCurTriNo;
                        }
                        continue; // ループの始めに戻る
                    }
                }

                {
                    // 次の要素へ進める
                    uint iTriNo1 = IndexRot3[1][iCurTriNo];
                    if (tris[(int)iCurTri].G2[iTriNo1] != -2 && tris[(int)iCurTri].G2[iTriNo1] != -3)
                    {
                        isWall = true;
                        break;
                    }
                    uint iNexTri = tris[(int)iCurTri].S2[iTriNo1];
                    uint[] nexRel = RelTriTri[tris[(int)iCurTri].R2[iTriNo1]];
                    uint iNexTriNo = nexRel[iCurTriNo];
                    System.Diagnostics.Debug.Assert(tris[(int)iNexTri].V[iNexTriNo] == iPt0);
                    if (iNexTri == iTri0)
                    {
                        break;   // 一周したら終わり
                    }
                    iCurTri = iNexTri;
                    iCurTriNo = iNexTriNo;
                }
            }
            if (!isWall)
            {
                return true;
            }

            ////////////////////////////////
            // 逆向きへの回転

            iCurTri = iTri0;
            iCurTriNo = iTriNo0;
            while (true)
            {
                System.Diagnostics.Debug.Assert(tris[(int)iCurTri].V[iCurTriNo] == iPt0);

                if (tris[(int)iCurTri].G2[iCurTriNo] == -2)
                {
                    // 向かいの要素を調べる
                    uint iDiaTri = tris[(int)iCurTri].S2[iCurTriNo];
                    uint[] diaRel = RelTriTri[tris[(int)iCurTri].R2[iCurTriNo]];
                    uint iDiaNoTri = diaRel[iCurTriNo];
                    System.Diagnostics.Debug.Assert(tris[(int)iDiaTri].G2[iDiaNoTri] == -2);
                    System.Diagnostics.Debug.Assert(tris[(int)iDiaTri].S2[iDiaNoTri] == iCurTri);
                    uint iDiaPt = tris[(int)iDiaTri].V[iDiaNoTri];
                    if (DetDelaunay(
                        points[(int)tris[(int)iCurTri].V[0]].Point,
                        points[(int)tris[(int)iCurTri].V[1]].Point,
                        points[(int)tris[(int)iCurTri].V[2]].Point,
                        points[(int)iDiaPt].Point) == 0)
                    {
                        // Delaunay条件が満たされない場合

                        // 辺を切り替える
                        FlipEdge(iCurTri, iCurTriNo, points, tris);

                        iCurTri = iDiaTri;
                        iCurTriNo = 1;
                        System.Diagnostics.Debug.Assert(tris[(int)iCurTri].V[iCurTriNo] == iPt0);
                        continue;   // ループの始めに戻る
                    }
                }

                { 
                    // 次の要素へ進める
                    uint iTriNo2 = IndexRot3[2][iCurTriNo];
                    if (tris[(int)iCurTri].G2[iTriNo2] != -2 && tris[(int)iCurTri].G2[iTriNo2] != -3)
                    {
                        return true;
                    }
                    uint iNexTri = tris[(int)iCurTri].S2[iTriNo2];
                    uint[] nexRel = RelTriTri[tris[(int)iCurTri].R2[iTriNo2]];
                    uint iNexTriNo = nexRel[iCurTriNo];
                    System.Diagnostics.Debug.Assert(tris[(int)iNexTri].V[iNexTriNo] == iPt0);
                    System.Diagnostics.Debug.Assert(iNexTri != iTri0);  // 一周したら終わり
                    iCurTri = iNexTri;
                    iCurTriNo = iNexTriNo;
                }
            }

            throw new InvalidOperationException();
            //return true;
        }

        public static bool FlipEdge(uint iTri0, uint iEd0, IList<MeshPoint2D> points, IList<MeshTri2D> tris)
        {
            System.Diagnostics.Debug.Assert(iTri0 < tris.Count);
            System.Diagnostics.Debug.Assert(iEd0 < 3);
            System.Diagnostics.Debug.Assert(tris[(int)iTri0].G2[iEd0] == -2);

            uint iTri1 = tris[(int)iTri0].S2[iEd0];
            uint iEd1 = RelTriTri[tris[(int)iTri0].R2[iEd0]][iEd0];
            System.Diagnostics.Debug.Assert(iTri1 < tris.Count);
            System.Diagnostics.Debug.Assert(iEd1 < 3);
            System.Diagnostics.Debug.Assert(tris[(int)iTri1].G2[iEd1] == -2);

            //	std::cout << itri0 << "-" << ied0 << "    " << itri1 << "-" << ied1 << std::endl;

            MeshTri2D old0 = new MeshTri2D(tris[(int)iTri0]);
            MeshTri2D old1 = new MeshTri2D(tris[(int)iTri1]);

            uint no00 = iEd0;
            uint no10 = TriElEdgeNo[iEd0][0];
            uint no20 = TriElEdgeNo[iEd0][1];

            uint no01 = iEd1;
            uint no11 = TriElEdgeNo[iEd1][0];
            uint no21 = TriElEdgeNo[iEd1][1];

            System.Diagnostics.Debug.Assert(old0.V[no10] == old1.V[no21]);
            System.Diagnostics.Debug.Assert(old0.V[no20] == old1.V[no11]);

            points[(int)old0.V[no10]].Elem = (int)iTri0;
            points[(int)old0.V[no10]].Dir = 0;
            points[(int)old0.V[no00]].Elem = (int)iTri0;
            points[(int)old0.V[no00]].Dir = 2;
            points[(int)old1.V[no11]].Elem = (int)iTri1;
            points[(int)old1.V[no11]].Dir = 0;
            points[(int)old1.V[no01]].Elem = (int)iTri1;
            points[(int)old1.V[no01]].Dir = 2;

            {
                MeshTri2D tri = tris[(int)iTri0];
                tri.V[0] = old0.V[no10];
                tri.V[1] = old1.V[no01];
                tri.V[2] = old0.V[no00];
                tri.G2[0] = -2;
                tri.G2[1] = old0.G2[no20];
                tri.G2[2] = old1.G2[no11];
                tri.S2[0] = iTri1;
                tri.S2[1] = old0.S2[no20];
                tri.S2[2] = old1.S2[no11];

                tri.R2[0] = 0;
                if (old0.G2[no20] == -2 || old0.G2[no20] == -3)
                {
                    System.Diagnostics.Debug.Assert(old0.R2[no20] < 3);
                    uint[] rel = RelTriTri[old0.R2[no20]];
                    System.Diagnostics.Debug.Assert(old0.S2[no20] < tris.Count);
                    System.Diagnostics.Debug.Assert(old0.S2[no20] != iTri0);
                    System.Diagnostics.Debug.Assert(old0.S2[no20] != iTri1);
                    tri.R2[1] = (uint)El2RelTriTri[rel[no10] * 3 + rel[no20]];
                    System.Diagnostics.Debug.Assert(tri.R2[1] >= 0 && tri.R2[1] < 3);
                    tris[(int)old0.S2[no20]].S2[rel[no20]] = iTri0;
                    tris[(int)old0.S2[no20]].R2[rel[no20]] = InvRelTriTri[tri.R2[1]];
                }
                if (old1.G2[no11] == -2 || old1.G2[no11] == -3)
                {
                    System.Diagnostics.Debug.Assert(old1.R2[no11] < 3);
                    uint[] rel = RelTriTri[old1.R2[no11]];
                    System.Diagnostics.Debug.Assert(old1.S2[no11] < tris.Count);
                    tri.R2[2] = (uint)El2RelTriTri[rel[no21] * 3 + rel[no01]];
                    System.Diagnostics.Debug.Assert(tri.R2[2] >= 0 && tri.R2[2] < 3);
                    tris[(int)old1.S2[no11]].S2[rel[no11]] = iTri0;
                    tris[(int)old1.S2[no11]].R2[rel[no11]] = InvRelTriTri[tri.R2[2]];
                }
            }

            {
                MeshTri2D tri = tris[(int)iTri1];
                tri.V[0] = old1.V[no11];
                tri.V[1] = old0.V[no00];
                tri.V[2] = old1.V[no01];
                tri.G2[0] = -2;
                tri.G2[1] = old1.G2[no21];
                tri.G2[2] = old0.G2[no10];
                tri.S2[0] = iTri0; tri.S2[1] = old1.S2[no21];
                tri.S2[2] = old0.S2[no10];

                tri.R2[0] = 0;
                if (old1.G2[no21] == -2 || old1.G2[no21] == -3)
                {
                    System.Diagnostics.Debug.Assert(old1.R2[no21] < 3);
                    uint[] rel = RelTriTri[old1.R2[no21]];
                    System.Diagnostics.Debug.Assert(old1.S2[no21] < tris.Count);
                    tri.R2[1] = (uint)El2RelTriTri[rel[no11] * 3 + rel[no21]];
                    System.Diagnostics.Debug.Assert(tri.R2[1] >= 0 && tri.R2[1] < 3);
                    tris[(int)old1.S2[no21]].S2[rel[no21]] = iTri1;
                    tris[(int)old1.S2[no21]].R2[rel[no21]] = InvRelTriTri[tri.R2[1]];
                }
                if (old0.G2[no10] == -2 || old0.G2[no10] == -3)
                {
                    System.Diagnostics.Debug.Assert(old0.R2[no10] < 3);
                    uint[] rel = RelTriTri[old0.R2[no10]];
                    System.Diagnostics.Debug.Assert(old0.S2[no10] < tris.Count);
                    tri.R2[2] = (uint)El2RelTriTri[rel[no20] * 3 + rel[no00]];
                    System.Diagnostics.Debug.Assert(tri.R2[2] >= 0 && tri.R2[2] < 3);
                    tris[(int)old0.S2[no10]].S2[rel[no10]] = iTri1;
                    tris[(int)old0.S2[no10]].R2[rel[no10]] = InvRelTriTri[tri.R2[2]];
                }
            }
            return true;
        }

        // 辺[iPt0-iPt1]の左側の３角形itri0を探索する
        // 三角形がなければ->falseを返す。
        // 三角形があれば  ->true を返す。
        // 但し、その場合
        // tri[iTri0].v[iTriNo0]==iPt0
        // tri[iTri0].v[iTriNo1]==iPt1
        // を満たす
        public static bool FindEdge(uint iPt0, uint iPt1, out uint iTri0, out uint iTriNo0, out uint iTriNo1,
            IList<MeshPoint2D> points, IList<MeshTri2D> tris)
        {
            iTri0 = 0;
            iTriNo0 = 0;
            iTriNo1 = 0;

            uint iIniTri = (uint)points[(int)iPt0].Elem;
            uint iIniTriNo = points[(int)iPt0].Dir;
            uint iCurTri = iIniTri;
            uint iCurTriNo = iIniTriNo;
            while (true)
            {
                //　時計周りに検索する。
                System.Diagnostics.Debug.Assert(tris[(int)iCurTri].V[iCurTriNo] == iPt0);
                { 
                    // この要素がOKか調べる
                    uint iTriNo2 = IndexRot3[1][iCurTriNo];
                    if (tris[(int)iCurTri].V[iTriNo2] == iPt1)
                    {
                        iTri0 = iCurTri;
                        iTriNo0 = iCurTriNo;
                        iTriNo1 = iTriNo2;
                        System.Diagnostics.Debug.Assert(tris[(int)iTri0].V[iTriNo0] == iPt0);
                        System.Diagnostics.Debug.Assert(tris[(int)iTri0].V[iTriNo1] == iPt1);
                        return true;
                    }
                }
                {  
                    // 次の要素へ進める
                    uint iTriNo2 = IndexRot3[2][iCurTriNo];
                    if (tris[(int)iCurTri].G2[iTriNo2] != -2 && tris[(int)iCurTri].G2[iTriNo2] != -3)
                    {
                        break;
                    }
                    uint iNexTri = tris[(int)iCurTri].S2[iTriNo2];
                    uint[] rel = RelTriTri[tris[(int)iCurTri].R2[iTriNo2]];
                    uint iTriNo3 = rel[iCurTriNo];
                    System.Diagnostics.Debug.Assert(tris[(int)iNexTri].V[iTriNo3] == iPt0);
                    if (iNexTri == iIniTri)
                    {
                        return false;
                    }
                    iCurTri = iNexTri;
                    iCurTriNo = iTriNo3;
                }
            }

            iCurTriNo = iIniTriNo;
            iCurTri = iIniTri;
            while (true)
            {   
                //　反時計周りの検索
                System.Diagnostics.Debug.Assert(tris[(int)iCurTri].V[iCurTriNo] == iPt0);
                {  
                    // 次の要素へ進める
                    uint iTriNo2 = IndexRot3[1][iCurTriNo];
                    if (tris[(int)iCurTri].G2[iTriNo2] != -2 || tris[(int)iCurTri].G2[iTriNo2] != -3)
                    {
                        break;
                    }
                    uint iNexTri = tris[(int)iCurTri].S2[iTriNo2];
                    uint[] rel = RelTriTri[tris[(int)iCurTri].R2[iTriNo2]];
                    uint iTriNo3 = rel[iCurTriNo];
                    System.Diagnostics.Debug.Assert(tris[(int)iNexTri].V[iTriNo3] == iPt0);
                    if (iNexTri == iIniTri)
                    {   
                        // 一周したら終わり
                        iTri0 = 0;
                        iTriNo0 = 0; iTriNo1 = 0;
                        return false;
                    }
                    iCurTri = iNexTri;
                    iCurTriNo = iTriNo3;
                }
                {
                    // 要素の向きを調べる
                    uint iTriNo2 = IndexRot3[1][iCurTriNo];
                    if (tris[(int)iCurTri].V[iTriNo2] == iPt1)
                    {
                        iTri0 = iCurTri;
                        iTriNo0 = iCurTriNo;
                        iTriNo1 = iTriNo2;
                        System.Diagnostics.Debug.Assert(tris[(int)iTri0].V[iTriNo0] == iPt0);
                        System.Diagnostics.Debug.Assert(tris[(int)iTri0].V[iTriNo1] == iPt1);
                        return true;
                    }
                }
            }

            return false;
        }

        public static bool FindEdgePointAcrossEdge(uint iPt0, uint iPt1, 
            out uint iTri0, out uint iTriNo0, out uint iTriNo1, out double ratio,
            IList<MeshPoint2D> points, IList<MeshTri2D> tris)
        {
            uint iIniTri = (uint)points[(int)iPt0].Elem;
            uint iIniTriNo = points[(int)iPt0].Dir;
            uint iCurTri = iIniTri;
            uint iCurTriNo = iIniTriNo;
            while (true)
            {
                //　反時計周りの検索
                System.Diagnostics.Debug.Assert(tris[(int)iCurTri].V[iCurTriNo] == iPt0);
                {
                    uint iTriNo2 = IndexRot3[1][iCurTriNo];
                    uint iTriNo3 = IndexRot3[2][iCurTriNo];
                    double area0 = CadUtils.TriArea(
                        points[(int)iPt0].Point,
                        points[(int)tris[(int)iCurTri].V[iTriNo2]].Point,
                        points[(int)iPt1].Point);
                    if (area0 > -1.0e-20)
                    {
                        double area1 = CadUtils.TriArea(
                            points[(int)iPt0].Point,
                            points[(int)iPt1].Point,
                            points[(int)tris[(int)iCurTri].V[iTriNo3]].Point);
                        if (area1 > -1.0e-20)
                        {
                            System.Diagnostics.Debug.Assert(area0 + area1 > 1.0e-20);
                            ratio = area0 / (area0 + area1);
                            iTri0 = iCurTri;
                            iTriNo0 = iTriNo2;
                            iTriNo1 = iTriNo3;
                            return true;
                        }
                    }
                }
                { 
                    // 次の要素へ進める
                    uint iTriNo2 = IndexRot3[1][iCurTriNo];
                    if (tris[(int)iCurTri].G2[iTriNo2] != -2 && tris[(int)iCurTri].G2[iTriNo2] != -3)
                    {
                        break;
                    }
                    uint iNexTri = tris[(int)iCurTri].S2[iTriNo2];
                    uint[] rel = RelTriTri[tris[(int)iCurTri].R2[iTriNo2]];
                    uint iTriNo3 = rel[iCurTriNo];
                    System.Diagnostics.Debug.Assert(tris[(int)iNexTri].V[iTriNo3] == iPt0);
                    if (iNexTri == iIniTri)
                    { 
                        // 一周したら終わり
                        iTri0 = 0;
                        iTriNo0 = 0; iTriNo1 = 0;
                        ratio = 0.0;
                        return false;
                    }
                    iCurTri = iNexTri;
                    iCurTriNo = iTriNo3;
                }
            }

            iCurTriNo = iIniTriNo;
            iCurTri = iIniTri;
            while (true)
            {  
                //　時計周りに検索する。
                System.Diagnostics.Debug.Assert(tris[(int)iCurTri].V[iCurTriNo] == iPt0);
                {
                    uint iTriNo2 = IndexRot3[1][iCurTriNo];
                    uint iTriNo3 = IndexRot3[2][iCurTriNo];
                    double area0 = CadUtils.TriArea(
                        points[(int)iPt0].Point,
                        points[(int)tris[(int)iCurTri].V[iTriNo2]].Point,
                        points[(int)iPt1].Point);
                    if (area0 > -1.0e-20)
                    {
                        double area1 = CadUtils.TriArea(
                            points[(int)iPt0].Point,
                            points[(int)iPt1].Point,
                            points[(int)tris[(int)iCurTri].V[iTriNo3]].Point);
                        if (area1 > -1.0e-20)
                        {
                            System.Diagnostics.Debug.Assert(area0 + area1 > 1.0e-20);
                            ratio = area0 / (area0 + area1);
                            iTri0 = iCurTri;
                            iTriNo0 = iTriNo2;
                            iTriNo1 = iTriNo3;
                            return true;
                        }
                    }
                }
                {
                    // 次の要素へ進める
                    uint iTriNo2 = IndexRot3[2][iCurTriNo];
                    if (tris[(int)iCurTri].G2[iTriNo2] != -2 && tris[(int)iCurTri].G2[iTriNo2] != -3)
                    {
                        break;
                    }
                    uint iNexTri = tris[(int)iCurTri].S2[iTriNo2];
                    uint[] rel = RelTriTri[tris[(int)iCurTri].R2[iTriNo2]];
                    uint iTriNo3 = rel[iCurTriNo];
                    System.Diagnostics.Debug.Assert(tris[(int)iNexTri].V[iTriNo3] == iPt0);
                    if (iNexTri == iIniTri)
                    {
                        System.Diagnostics.Debug.Assert(false);  // 一周しないはず
                    }
                    iCurTri = iNexTri;
                    iCurTriNo = iTriNo3;
                }
            }

            // 失敗したときの値を入れる
            iTri0 = 0;
            iTriNo0 = 0;
            iTriNo1 = 0;
            ratio = 0.0;

            return false;
        }


        public static bool MakePointSurTri(IList<MeshTri2D> tris, uint npoin, uint[] elsupInd,
            out uint nelsup, out uint[] elsup )
        {

            uint nnotri = 3;

            for (int ipoin = 0; ipoin < npoin + 1; ipoin++)
            {
                elsupInd[ipoin] = 0;
            }
            for (int itri = 0; itri < tris.Count; itri++)
            {
                for (int inotri = 0; inotri < nnotri; inotri++)
                {
                    elsupInd[tris[itri].V[inotri] + 1]++;
                }
            }
            for (int ipoin = 0; ipoin < npoin; ipoin++)
            {
                elsupInd[ipoin + 1] += elsupInd[ipoin];
            }
            nelsup = elsupInd[npoin];
            elsup = new uint[nelsup];
            for (int itri = 0; itri < tris.Count; itri++)
            {
                for (int inotri = 0; inotri < nnotri; inotri++)
                {
                    uint ipoin0 = tris[itri].V[inotri];
                    uint ielsup = elsupInd[ipoin0];
                    elsup[ielsup] = (uint)itri;
                    elsupInd[ipoin0]++;
                }
            }
            for (int ipoin = (int)npoin; ipoin > 0; ipoin--)
            {
                elsupInd[ipoin] = elsupInd[ipoin - 1];
            }
            elsupInd[0] = 0;

            return true;
        }

        public static bool MakeInnerRelationTri(IList<MeshTri2D> tris, uint npoin,
            uint[] elsupInd, uint nelsup, uint[] elsup)
        {
            uint[][] EdEd2Rel = {
                new uint[(int)TriEdNo]{ 0, 2, 1 },
                new uint[(int)TriEdNo]{ 2, 1, 0 },
                new uint[(int)TriEdNo]{ 1, 0, 2 }
            };
            System.Diagnostics.Debug.Assert(EdEd2Rel.Length == TriEdNo);

            uint[] tmpPoin = new uint[npoin];
            for (int ipoin = 0; ipoin < npoin; ipoin++)
            {
                tmpPoin[ipoin] = 0;
            }
            uint[] inpofa = new uint[2];

            int nTri = tris.Count;
            for (int itri = 0; itri < nTri; itri++)
            {
                for (int iedtri = 0; iedtri < TriEdNo; iedtri++)
                {
                    for (int ipoed = 0; ipoed < EdNo; ipoed++)
                    {
                        inpofa[ipoed] = tris[itri].V[TriElEdgeNo[iedtri][ipoed]];
                        tmpPoin[inpofa[ipoed]] = 1;
                    }
                    uint ipoin0 = inpofa[0];
                    bool iflg = false;
                    for (uint ielsup = elsupInd[ipoin0]; ielsup < elsupInd[ipoin0 + 1]; ielsup++)
                    {
                        uint jtri0 = elsup[ielsup];
                        if (jtri0 == itri)
                        {
                            continue;
                        }
                        for (int jedtri = 0; jedtri < TriEdNo; jedtri++)
                        {
                            iflg = true;
                            for (int jpoed = 0; jpoed < EdNo; jpoed++)
                            {
                                uint jpoin0 = tris[(int)jtri0].V[TriElEdgeNo[jedtri][jpoed]];
                                if (tmpPoin[jpoin0] == 0)
                                {
                                    iflg = false;
                                    break;
                                }
                            }
                            if (iflg)
                            {
                                tris[itri].G2[iedtri] = -2;
                                tris[itri].S2[iedtri] = jtri0;
                                tris[itri].R2[iedtri] = EdEd2Rel[iedtri][jedtri];
                                break;
                            }
                        }
                        if (iflg)
                        {
                            break;
                        }
                    }
                    if (!iflg)
                    {
                        tris[itri].G2[iedtri] = -1;
                    }
                    for (int ipofa = 0; ipofa < EdNo; ipofa++)
                    {
                        tmpPoin[inpofa[ipofa]] = 0;
                    }
                }
            }
            return true;
        }

        public static void LaplacianSmoothing(IList<MeshPoint2D> points, IList<MeshTri2D> tris, IList<uint> isntMoves)
        {
            for (uint iPt = 0; iPt < points.Count; iPt++)
            {  
                // 点周りの点を探索して調べる。
                if (iPt < isntMoves.Count)
                {
                    if (isntMoves[(int)iPt] == 1)
                    {
                        continue;
                    }
                }
                uint iTriIni = (uint)points[(int)iPt].Elem;
                uint iNoElCIni = points[(int)iPt].Dir;
                System.Diagnostics.Debug.Assert(iTriIni < tris.Count);
                System.Diagnostics.Debug.Assert(iNoElCIni < 3);
                System.Diagnostics.Debug.Assert(tris[(int)iTriIni].V[iNoElCIni] == iPt);
                uint iTri0 = iTriIni;
                uint iNoElC0 = iNoElCIni;
                uint iNoElB0 = TriElEdgeNo[iNoElC0][0];
                bool isBound = false;
                OpenTK.Vector2d vecDelta = points[(int)iPt].Point;
                uint ntriAround = 1;
                while (true)
                {
                    System.Diagnostics.Debug.Assert(iTri0 < tris.Count);
                    System.Diagnostics.Debug.Assert(iNoElC0 < 3);
                    System.Diagnostics.Debug.Assert(tris[(int)iTri0].V[iNoElC0] == iPt);
                    {
                        vecDelta.X += points[(int)tris[(int)iTri0].V[iNoElB0]].Point.X;
                        vecDelta.Y += points[(int)tris[(int)iTri0].V[iNoElB0]].Point.Y;
                        ntriAround++;
                    }
                    if (tris[(int)iTri0].G2[iNoElB0] == -2)
                    {
                        uint iTri1 = tris[(int)iTri0].S2[iNoElB0];
                        uint rel01 = tris[(int)iTri0].R2[iNoElB0];
                        uint iNoElC1 = RelTriTri[rel01][iNoElC0];
                        uint iNoElB1 = RelTriTri[rel01][TriElEdgeNo[iNoElC0][1]];
                        System.Diagnostics.Debug.Assert(iTri1 < tris.Count);
                        System.Diagnostics.Debug.Assert(
                            tris[(int)iTri1].S2[RelTriTri[rel01][iNoElB0]] == iTri0);
                        System.Diagnostics.Debug.Assert(tris[(int)iTri1].V[iNoElC1] == iPt);
                        if (iTri1 == iTriIni)
                        {
                            break;
                        }
                        iTri0 = iTri1;
                        iNoElC0 = iNoElC1;
                        iNoElB0 = iNoElB1;
                    }
                    else
                    {   
                        // この点は境界上の点だから動かしてはならない。
                        isBound = true;
                        break;
                    }
                }
                if (isBound)
                {
                    continue;
                }
                points[(int)iPt].Point = new OpenTK.Vector2d(
                    vecDelta.X / ntriAround,
                    vecDelta.Y / ntriAround);
            }
        }

        public static void LaplaceDelaunaySmoothing(
            IList<MeshPoint2D> points, IList<MeshTri2D> tris, IList<uint> isntMoves)
        {
            for (uint iPt = 0; iPt < points.Count; iPt++)
            {
                // 点周りの点を探索して調べる。
                if (iPt < isntMoves.Count)
                {
                    if (isntMoves[(int)iPt] == 1)
                    {
                        continue;
                    }
                }
                uint iTriIni = (uint)points[(int)iPt].Elem;
                uint iNoElCIni = points[(int)iPt].Dir;
                System.Diagnostics.Debug.Assert(iTriIni < tris.Count);
                System.Diagnostics.Debug.Assert(iNoElCIni < 3);
                System.Diagnostics.Debug.Assert(tris[(int)iTriIni].V[iNoElCIni] == iPt);
                uint iTri0 = iTriIni;
                uint iNoElC0 = iNoElCIni;
                uint iNoElB0 = TriElEdgeNo[iNoElC0][0];
                bool isBound = false;
                OpenTK.Vector2d vecDelta = points[(int)iPt].Point;
                uint ntriAround = 1;
                while (true)
                { 
                    // 点の周りの要素を一回りする
                    System.Diagnostics.Debug.Assert(iTri0 < tris.Count);
                    System.Diagnostics.Debug.Assert(iNoElC0 < 3);
                    System.Diagnostics.Debug.Assert(tris[(int)iTri0].V[iNoElC0] == iPt);
                    {
                        vecDelta.X += points[(int)tris[(int)iTri0].V[iNoElB0]].Point.X;
                        vecDelta.Y += points[(int)tris[(int)iTri0].V[iNoElB0]].Point.Y;
                        ntriAround++;
                    }
                    if (tris[(int)iTri0].G2[iNoElB0] == -2)
                    {
                        uint iTri1 = tris[(int)iTri0].S2[iNoElB0];
                        uint rel01 = tris[(int)iTri0].R2[iNoElB0];
                        uint iNoElC1 = RelTriTri[rel01][iNoElC0];
                        uint iNoElB1 = RelTriTri[rel01][TriElEdgeNo[iNoElC0][1]];
                        System.Diagnostics.Debug.Assert(iTri1 < tris.Count);
                        System.Diagnostics.Debug.Assert(
                            tris[(int)iTri1].S2[RelTriTri[rel01][iNoElB0]] == iTri0);
                        System.Diagnostics.Debug.Assert(tris[(int)iTri1].V[iNoElC1] == iPt);
                        if (iTri1 == iTriIni)
                        {
                            break;
                        }
                        iTri0 = iTri1;
                        iNoElC0 = iNoElC1;
                        iNoElB0 = iNoElB1;
                    }
                    else
                    {
                        // この点は境界上の点だから動かしてはならない。
                        isBound = true;
                        break;
                    }
                }
                if (isBound)
                {
                    continue;
                }
                points[(int)iPt].Point = new OpenTK.Vector2d(
                    vecDelta.X / ntriAround,
                    vecDelta.Y / ntriAround);
                DelaunayAroundPoint(iPt, points, tris);
            }
        }


    }
}
