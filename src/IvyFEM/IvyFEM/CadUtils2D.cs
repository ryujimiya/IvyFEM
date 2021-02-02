using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class CadUtils2D
    {
        /// <summary>
        /// 三角形最小面積
        /// </summary>
        public const double DefMinTriArea = 1.0e-10;

        public static OpenTK.Vector2d GetDirection(OpenTK.Vector2d v1, OpenTK.Vector2d v2)
        {
            OpenTK.Vector2d dir = v2 - v1;
            dir = OpenTK.Vector2d.Normalize(dir);
            return dir;
        }

        public static OpenTK.Vector2d GetNormal(OpenTK.Vector2d v1, OpenTK.Vector2d v2)
        {
            // Note: tはdirの逆方向
            var t = v1 - v2;
            t = OpenTK.Vector2d.Normalize(t);
            // n = e3 x t
            OpenTK.Vector2d normal = new OpenTK.Vector2d(-t.Y, t.X);
            return normal;
        }

        public static double Dot2D(double[] v1, double[] v2)
        {
            double ret = v1[0] * v2[0] + v1[1] * v2[1];
            return ret;
        }

        public static double[] Normalize2D(double[] v)
        {
            double len = Math.Sqrt(v[0] * v[0] + v[1] * v[1]);
            double[] ret = { v[0] / len, v[1] / len };
            return ret;
        }

        public static double GetDistance2D(double[] v1, double[] v2)
        {
            double[] vec = { v2[0] - v1[0], v2[1] - v1[1] };
            double d = Math.Sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
            return d;
        }

        public static double[] GetDirection2D(double[] v1, double[] v2)
        {
            double[] dir = { v2[0] - v1[0], v2[1] - v1[1] };
            dir = Normalize2D(dir);
            return dir;
        }

        public static double[] GetNormal2D(double[] v1, double[] v2)
        {
            double[] tan = { v1[0] - v2[0], v1[1] - v2[1] };
            tan = Normalize2D(tan);
            // n = e3 x t
            double[] normal = new double[] { -tan[1], tan[0] };
            return normal;
        }

        public static double Cross(OpenTK.Vector2d v1, OpenTK.Vector2d v2)
        {
            double cross = v1.X * v2.Y - v1.Y * v2.X;
            return cross;
        }

        public static OpenTK.Vector2d GetProjectedPointOnCircle(OpenTK.Vector2d c, double r, OpenTK.Vector2d v)
        {
            OpenTK.Vector2d cv = v - c;
            double k = (r / cv.Length);
            return cv * k + c;
        }

        public static double SquareDistance(OpenTK.Vector2d pt0, OpenTK.Vector2d pt1)
        {
            OpenTK.Vector2d v = pt1 - pt0;
            double len = v.Length;
            return len * len;
        }

        public static double SquareLength(OpenTK.Vector2d point)
        {
            double len = point.Length;
            return len * len;
        }

        public static double TriHeight(OpenTK.Vector2d v1, OpenTK.Vector2d v2, OpenTK.Vector2d v3)
        {
            double area = TriArea(v1, v2, v3);
            double len = Math.Sqrt(SquareDistance(v2, v3));
            return area * 2.0 / len;
        }

        public static double TriArea(OpenTK.Vector2d v1, OpenTK.Vector2d v2, OpenTK.Vector2d v3)
        {
            return 0.5 * ((v2.X - v1.X) * (v3.Y - v1.Y) - (v3.X - v1.X) * (v2.Y - v1.Y));
        }

        public static void UnitNormalAreaTri3D(out double[] n, out double a,
            OpenTK.Vector3d v1, OpenTK.Vector3d v2, OpenTK.Vector3d v3)
        {
            n = new double[3];
            n[0] = (v2.Y - v1.Y) * (v3.Z - v1.Z) - (v3.Y - v1.Y) * (v2.Z - v1.Z);
            n[1] = (v2.Z - v1.Z) * (v3.X - v1.X) - (v3.Z - v1.Z) * (v2.X - v1.X);
            n[2] = (v2.X - v1.X) * (v3.Y - v1.Y) - (v3.X - v1.X) * (v2.Y - v1.Y);
            a = Math.Sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) * 0.5;
            double invA = 0.5 / a;
            n[0] *= invA;
            n[1] *= invA;
            n[2] *= invA;
        }

        public static double SquareCircumradius(OpenTK.Vector2d p0, OpenTK.Vector2d p1, OpenTK.Vector2d p2)
        {
            double area = TriArea(p0, p1, p2);

            double dtmp0 = SquareDistance(p1, p2);
            double dtmp1 = SquareDistance(p0, p2);
            double dtmp2 = SquareDistance(p0, p1);

            return dtmp0 * dtmp1 * dtmp2 / (16.0 * area * area);
        }

        public static bool CenterCircumcircle(OpenTK.Vector2d p0, OpenTK.Vector2d p1, OpenTK.Vector2d p2,
            double minTriArea,
            out OpenTK.Vector2d center)
        {
            center = new OpenTK.Vector2d();

            double area = TriArea(p0, p1, p2);
            if (Math.Abs(area) < minTriArea)
            {
                return false;
            }
            double tmpVal = 1.0 / (area * area * 16.0);

            double dtmp0 = SquareDistance(p1, p2);
            double dtmp1 = SquareDistance(p0, p2);
            double dtmp2 = SquareDistance(p0, p1);

            double etmp0 = tmpVal * dtmp0 * (dtmp1 + dtmp2 - dtmp0);
            double etmp1 = tmpVal * dtmp1 * (dtmp0 + dtmp2 - dtmp1);
            double etmp2 = tmpVal * dtmp2 * (dtmp0 + dtmp1 - dtmp2);

            center = new OpenTK.Vector2d(
                etmp0 * p0.X + etmp1 * p1.X + etmp2 * p2.X,
                etmp0 * p0.Y + etmp1 * p1.Y + etmp2 * p2.Y);
            return true;
        }

        public static bool IsCrossLineSegLineSeg(
            OpenTK.Vector2d sPt0, OpenTK.Vector2d ePt0, OpenTK.Vector2d sPt1, OpenTK.Vector2d ePt1)
        {
            double minX0 = (sPt0.X < ePt0.X) ? sPt0.X : ePt0.X;
            double maxX0 = (sPt0.X > ePt0.X) ? sPt0.X : ePt0.X;
            double maxX1 = (sPt1.X > ePt1.X) ? sPt1.X : ePt1.X;
            double minX1 = (sPt1.X < ePt1.X) ? sPt1.X : ePt1.X;
            double minY0 = (sPt0.Y < ePt0.Y) ? sPt0.Y : ePt0.Y;
            double maxY0 = (sPt0.Y > ePt0.Y) ? sPt0.Y : ePt0.Y;
            double maxY1 = (sPt1.Y > ePt1.Y) ? sPt1.Y : ePt1.Y;
            double minY1 = (sPt1.Y < ePt1.Y) ? sPt1.Y : ePt1.Y;
            double len = ((maxX0 - minX0) + (maxY0 - minY0) + (maxX1 - minX1) + (maxY1 - minY1)) * 0.0001;
            if (maxX1 + len < minX0)
            {
                return false;
            }
            if (maxX0 + len < minX1)
            {
                return false;
            }
            if (maxY1 + len < minY0)
            {
                return false;
            }
            if (maxY0 + len < minY1)
            {
                return false;
            }

            double area1 = TriArea(sPt0, ePt0, sPt1);
            double area2 = TriArea(sPt0, ePt0, ePt1);
            double area3 = TriArea(sPt1, ePt1, sPt0);
            double area4 = TriArea(sPt1, ePt1, ePt0);
            double a12 = area1 * area2;
            if (a12 > 0)
            {
                return false;
            }
            double a34 = area3 * area4;
            if (a34 > 0)
            {
                return false;
            }
            return true;
        }

        public static bool IsCrossCircleCircle(
            OpenTK.Vector2d cPt0, double radius0,
            OpenTK.Vector2d cPt1, double radius1,
            out OpenTK.Vector2d pt0, out OpenTK.Vector2d pt1)
        {
            pt0 = new OpenTK.Vector2d();
            pt1 = new OpenTK.Vector2d();

            double sqDist = SquareDistance(cPt0, cPt1);
            double dist = Math.Sqrt(sqDist);
            if (radius0 + radius1 < dist)
            {
                return false;
            }
            if (Math.Abs(radius0 - radius1) > dist)
            {
                return false;
            }
            if (dist < 1.0e-30)
            {
                return false;
            }
            double ct = 0.5 * (sqDist + radius0 * radius0 - radius1 * radius1) / (radius0 * dist);
            System.Diagnostics.Debug.Assert(ct >= -1 && ct <= 1);
            double st = Math.Sqrt(1 - ct * ct);
            OpenTK.Vector2d e0 = (cPt1 - cPt0) * (1 / dist);
            OpenTK.Vector2d e1 = new OpenTK.Vector2d(e0.Y, -e0.X);
            pt0 = cPt0 + e0 * (radius0 * ct) + e1 * (radius0 * st);
            pt1 = cPt0 + e0 * (radius0 * ct) - e1 * (radius0 * st);
            return true;
        }

        public static bool IsCrossLineCircle(
            OpenTK.Vector2d cPt, double radius, OpenTK.Vector2d sPt, OpenTK.Vector2d ePt,
            out double t0, out double t1)
        {
            t0 = 0.0f;
            t1 = 0.0f;

            double minX = (sPt.X < ePt.X) ? sPt.X : ePt.X;
            if (cPt.X + radius < minX)
            {
                return false;
            }
            double maxX = (sPt.X > ePt.X) ? sPt.X : ePt.X;
            if (cPt.X - radius > maxX)
            {
                return false;
            }
            double minY = (sPt.Y < ePt.Y) ? sPt.Y : ePt.Y;
            if (cPt.Y + radius < minY)
            {
                return false;
            }
            double maxY = (sPt.Y > ePt.Y) ? sPt.Y : ePt.Y;
            if (cPt.Y - radius > maxY)
            {
                return false;
            }

            OpenTK.Vector2d es = ePt - sPt;
            OpenTK.Vector2d cs = cPt - sPt;
            double a = SquareLength(es);
            double b = OpenTK.Vector2d.Dot(es, cs);
            double c = SquareLength(cs) - radius * radius;
            double det = b * b - a * c;
            if (det < 0)
            {
                return false;
            }
            t0 = (b - Math.Sqrt(det)) / a;
            t1 = (b + Math.Sqrt(det)) / a;
            return true;
        }

        public static double FindNearestPointParameterLinePoint(
            OpenTK.Vector2d cPt, OpenTK.Vector2d sPt, OpenTK.Vector2d ePt)
        {
            OpenTK.Vector2d es = ePt - sPt;
            OpenTK.Vector2d sc = sPt - cPt;
            double a = SquareLength(es);
            double b = OpenTK.Vector2d.Dot(es, sc);
            return -b / a;
        }

        public static bool IsDirectionArc(
            OpenTK.Vector2d pt, OpenTK.Vector2d sPt, OpenTK.Vector2d ePt, OpenTK.Vector2d cPt, bool isLeftSide)
        {
            if (isLeftSide)
            {
                if (TriArea(sPt, cPt, ePt) > 0.0)
                {
                    if (TriArea(sPt, cPt, pt) > 0.0 && TriArea(pt, cPt, ePt) > 0.0)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                else
                {
                    if (TriArea(sPt, cPt, pt) > 0.0 || TriArea(pt, cPt, ePt) > 0.0)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
            }
            else
            {
                if (TriArea(ePt, cPt, sPt) > 0.0)
                {
                    if (TriArea(ePt, cPt, pt) > 0.0 && TriArea(pt, cPt, sPt) > 0.0)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                else
                {
                    if (TriArea(ePt, cPt, pt) > 0.0 || TriArea(pt, cPt, sPt) > 0.0)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
            }

            throw new InvalidOperationException();
            //return true;
        }


        public static double GetDistancePointArc(OpenTK.Vector2d pt, OpenTK.Vector2d sPt1, OpenTK.Vector2d ePt1,
            OpenTK.Vector2d cPt1, double radius1, bool isLeftSide1)
        {
            double minDist = OpenTK.Vector2d.Distance(pt, sPt1);
            double d0 = OpenTK.Vector2d.Distance(pt, ePt1);
            minDist = (minDist < d0) ? minDist : d0;
            if (IsDirectionArc(GetProjectedPointOnCircle(cPt1, radius1, pt), sPt1, ePt1, cPt1, isLeftSide1))
            {
                d0 = Math.Abs(OpenTK.Vector2d.Distance(pt, cPt1) - radius1);
                minDist = (d0 < minDist) ? d0 : minDist;
            }
            return minDist;
        }

        public static double GetDistanceLineSegLineSeg(
            OpenTK.Vector2d sPt0, OpenTK.Vector2d ePt0, OpenTK.Vector2d sPt1, OpenTK.Vector2d ePt1)
        {
            if (IsCrossLineSegLineSeg(sPt0, ePt0, sPt1, ePt1))
            {
                return -1;
            }
            double sD1 = GetDistanceLineSegPoint(sPt0, sPt1, ePt1);
            double eD1 = GetDistanceLineSegPoint(ePt0, sPt1, ePt1);
            double sD0 = GetDistanceLineSegPoint(sPt1, sPt0, ePt0);
            double eD0 = GetDistanceLineSegPoint(ePt1, sPt0, ePt0);
            double minDist = sD1;
            minDist = (eD1 < minDist) ? eD1 : minDist;
            minDist = (sD0 < minDist) ? sD0 : minDist;
            minDist = (eD0 < minDist) ? eD0 : minDist;
            return minDist;
        }

        public static double GetDistanceLineSegPoint(OpenTK.Vector2d cPt, OpenTK.Vector2d sPt, OpenTK.Vector2d ePt)
        {
            OpenTK.Vector2d es = ePt - sPt;
            OpenTK.Vector2d sc = sPt - cPt;
            double a = SquareLength(es);
            double b = OpenTK.Vector2d.Dot(es, sc);
            double t = -b / a;
            if (t < 0)
            {
                return OpenTK.Vector2d.Distance(sPt, cPt);
            }
            if (t > 1)
            {
                return OpenTK.Vector2d.Distance(ePt, cPt);
            }
            OpenTK.Vector2d p = sPt + t * (ePt - sPt);
            return OpenTK.Vector2d.Distance(p, cPt);
        }

        public static double GetDistanceLineSegArc(
            OpenTK.Vector2d sPt0, OpenTK.Vector2d ePt0,
            OpenTK.Vector2d sPt1, OpenTK.Vector2d ePt1,
            OpenTK.Vector2d cPt1, double radius1, bool isLeftSide1)
        {
            double t0 = 0;
            double t1 = 0;
            if (IsCrossLineCircle(cPt1, radius1, sPt0, ePt0, out t0, out t1))
            {
                if (0 < t0 && t0 < 1 &&
                    IsDirectionArc(sPt0 + (ePt0 - sPt0) * t0, sPt1, ePt1, cPt1, isLeftSide1))
                {
                    return -1;
                }
                if (0 < t1 && t1 < 1 &&
                    IsDirectionArc(sPt0 + (ePt0 - sPt0) * t1, sPt1, ePt1, cPt1, isLeftSide1))
                {
                    return -1;
                }
            }
            double sminDist0 = GetDistancePointArc(sPt0, sPt1, ePt1, cPt1, radius1, isLeftSide1);
            double eminDist0 = GetDistancePointArc(ePt0, sPt1, ePt1, cPt1, radius1, isLeftSide1);
            double minDist = (sminDist0 < eminDist0) ? sminDist0 : eminDist0;
            double t = FindNearestPointParameterLinePoint(cPt1, sPt0, ePt0);
            if (t > 0 && t < 1)
            {
                OpenTK.Vector2d v = sPt0 + (ePt0 - sPt0) * t;
                double d0 = OpenTK.Vector2d.Distance(v, cPt1) - radius1;
                if (d0 > 0)
                {
                    if (IsDirectionArc(GetProjectedPointOnCircle(cPt1, radius1, v), sPt1, ePt1, cPt1, isLeftSide1))
                    {
                        minDist = (d0 < minDist) ? d0 : minDist;
                    }
                }
            }
            return minDist;
        }

        public static int CheckEdgeIntersection(IList<Edge2D> edges)
        {
            uint edgeCnt = (uint)edges.Count;
            for (int iedge = 0; iedge < edgeCnt; iedge++)
            {
                Edge2D iE = edges[iedge];
                if (iE.IsCrossEdgeSelf())
                {
                    return 1;
                }
                uint iPt0 = iE.GetVertexId(true);
                uint iPt1 = iE.GetVertexId(false);
                BoundingBox2D iBB = iE.GetBoundingBox();
                for (int jedge = iedge + 1; jedge < edgeCnt; jedge++)
                {
                    Edge2D jE = edges[jedge];
                    uint jPt0 = jE.GetVertexId(true);
                    uint jPt1 = jE.GetVertexId(false);
                    if (iPt0 != jPt0 && iPt0 != jPt1 &&
                        iPt1 != jPt0 &&  iPt1 != jPt1)
                    {
                        BoundingBox2D jBB = jE.GetBoundingBox();
                        if (jBB.MinX > iBB.MaxX || jBB.MaxX < iBB.MinX)
                        {
                            continue;
                        }
                        if (jBB.MinY > iBB.MaxY || jBB.MaxY < iBB.MinY) continue;
                        if (!iE.IsCrossEdge(jE))
                        {
                            continue;
                        }
                        return 1;
                    }
                    if (iPt0 == jPt0 && iPt1 == jPt1)
                    {
                        if (iE.IsCrossEdgeShareBothPoints(jE, true))
                        {
                            return 1;
                        }
                    }
                    else if (iPt0 == jPt1 && iPt1 == jPt0)
                    {
                        if (iE.IsCrossEdgeShareBothPoints(jE, false))
                        {
                            return 1;
                        }
                    }
                    else if (iPt0 == jPt0)
                    {
                        if (iE.IsCrossEdgeShareOnePoint(jE, true, true))
                        {
                            return 1;
                        }
                    }
                    else if (iPt0 == jPt1)
                    {
                        if (iE.IsCrossEdgeShareOnePoint(jE, true, false))
                        {
                            return 1;
                        }
                    }
                    else if (iPt1 == jPt0)
                    {
                        if (iE.IsCrossEdgeShareOnePoint(jE, false, true))
                        {
                            return 1;
                        }
                    }
                    else if (iPt1 == jPt1)
                    {
                        if (iE.IsCrossEdgeShareOnePoint(jE, false, false))
                        {
                            return 1;
                        }
                    }
                    continue;
                }
            }
            return 0;
        }

        public static void RotMatrix33(OpenTK.Quaterniond quat, double[] m)
        {
            double real = quat.W;
            double x = quat.X;
            double y = quat.Y;
            double z = quat.Z;

            m[0] = 1.0 - 2.0 * (y * y + z * z);
            m[1] = 2.0 * (x * y - z * real);
            m[2] = 2.0 * (z * x + y * real);

            m[3] = 2.0 * (x * y + z * real);
            m[4] = 1.0 - 2.0 * (z * z + x * x);
            m[5] = 2.0 * (y * z - x * real);

            m[6] = 2.0 * (z * x - y * real);
            m[7] = 2.0 * (y * z + x * real);
            m[8] = 1.0 - 2.0 * (y * y + x * x);
        }

        /// <summary>
        /// 回転移動する
        /// </summary>
        /// <param name="srcPt"></param>
        /// <param name="rotAngle"></param>
        /// <param name="rotOrigin"></param>
        /// <returns></returns>
        public static double[] GetRotCoord2D(double[] srcPt, double rotAngle, double[] rotOrigin = null)
        {
            System.Diagnostics.Debug.Assert(srcPt.Length == 2);
            if (rotOrigin != null)
            {
                System.Diagnostics.Debug.Assert(rotOrigin.Length == 2);
            }
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

        public static bool IsPointInsideTriangle(
            OpenTK.Vector2d p, OpenTK.Vector2d v1, OpenTK.Vector2d v2, OpenTK.Vector2d v3)
        {
            var vec1 = v3 - v1;
            var vec2 = v2 - v1;
            var vec3 = p - v1;

            double dot11 = OpenTK.Vector2d.Dot(vec1, vec1);
            double dot12 = OpenTK.Vector2d.Dot(vec1, vec2);
            double dot13 = OpenTK.Vector2d.Dot(vec1, vec3);
            double dot22 = OpenTK.Vector2d.Dot(vec2, vec2);
            double dot23 = OpenTK.Vector2d.Dot(vec2, vec3);

            // Compute barycentric coordinates
            double invDenom = 1 / (dot11 * dot22 - dot12 * dot12);
            double u = (dot22 * dot13 - dot12 * dot23) * invDenom;
            double v = (dot11 * dot23 - dot12 * dot13) * invDenom;

            // Check if point is in triangle
            return (u >= 0) && (v >= 0) && (u + v < 1);
        }

        public static bool GetTwoLineSegIntersectPoint(
            OpenTK.Vector2d line1P1, OpenTK.Vector2d line1P2,
            OpenTK.Vector2d line2P1, OpenTK.Vector2d line2P2,
            out OpenTK.Vector2d intersectP)
        {
            intersectP = new OpenTK.Vector2d();

            double a1 = line1P2.Y - line1P1.Y;
            double b1 = line1P1.X - line1P2.X;
            double c1 = a1 * line1P1.X + b1 * line1P1.Y;

            double a2 = line2P2.Y - line2P1.Y;
            double b2 = line2P1.X - line2P2.X;
            double c2 = a2 * line2P1.X + b2 * line2P1.Y;
            double eps = Constants.PrecisionLowerLimit;

            //lines are parallel
            double det = a1 * b2 - a2 * b1;
            if (Math.Abs(det) < eps)
            {
                //parallel lines
                return false;
            }
            else
            {
                double x = (b2 * c1 - b1 * c2) / det;
                double y = (a1 * c2 - a2 * c1) / det;
                bool online1 = (
                    (Math.Min(line1P1.X, line1P2.X) < x || Math.Abs(Math.Min(line1P1.X, line1P2.X) - x) < eps) &&
                    (Math.Max(line1P1.X, line1P2.X) > x || Math.Abs(Math.Max(line1P1.X, line1P2.X) - x) < eps) &&
                    (Math.Min(line1P1.Y, line1P2.Y) < y || Math.Abs(Math.Min(line1P1.Y, line1P2.Y) - y) < eps) &&
                    (Math.Max(line1P1.Y, line1P2.Y) > y || Math.Abs(Math.Max(line1P1.Y, line1P2.Y) - y) < eps));
                bool online2 = (
                    (Math.Min(line2P1.X, line2P2.X) < x || Math.Abs(Math.Min(line2P1.X, line2P2.X) - x) < eps) &&
                    (Math.Max(line2P1.X, line2P2.X) > x || Math.Abs(Math.Max(line2P1.X, line2P2.X) - x) < eps) &&
                    (Math.Min(line2P1.Y, line2P2.Y) < y || Math.Abs(Math.Min(line2P1.Y, line2P2.Y) - y) < eps) &&
                    (Math.Max(line2P1.Y, line2P2.Y) > y || Math.Abs(Math.Max(line2P1.Y, line2P2.Y) - y) < eps)
                    );

                if (online1 && online2)
                {
                    intersectP = new OpenTK.Vector2d(x, y);
                    return true;
                }
            }

            //intersection is at out of at least one segment.
            return false; 
        }

        public static bool IsPointInsidePolygon(OpenTK.Vector2d test, OpenTK.Vector2d[] polyPoints)
        {
            bool result = false;
            for (int i = 0; i < polyPoints.Length; i++)
            {
                int j = i == 0 ? (polyPoints.Length - 1) : (i - 1);
                if ((polyPoints[i].Y > test.Y) != (polyPoints[j].Y > test.Y) &&
                    (test.X < (polyPoints[j].X - polyPoints[i].X) * 
                    (test.Y - polyPoints[i].Y) / (polyPoints[j].Y - polyPoints[i].Y) + polyPoints[i].X))
                {
                    result = !result;
                }
            }
            return result;
        }

        public static OpenTK.Vector2d[] GetLineSegPolyIntersectPoints(
            OpenTK.Vector2d line1P1, OpenTK.Vector2d line1P2, OpenTK.Vector2d[] polyPoints)
        {
            IList<OpenTK.Vector2d> intersectionPoints = new List<OpenTK.Vector2d>();
            for (int i = 0; i < polyPoints.Length; i++)
            {
                int next = (i + 1 == polyPoints.Length) ? 0 : i + 1;
                OpenTK.Vector2d intersectP;
                bool isIntersect = GetTwoLineSegIntersectPoint(
                    line1P1, line1P2, polyPoints[i], polyPoints[next], out intersectP);
                if (isIntersect)
                {
                    intersectionPoints.Add(intersectP);
                }

            }
            return intersectionPoints.ToArray();
        }

        private static void _AddPoint(IList<OpenTK.Vector2d> points, OpenTK.Vector2d newPoint)
        {
            bool found = false;
            foreach (OpenTK.Vector2d p in points)
            {
                if ((p - newPoint).Length < Constants.PrecisionLowerLimit)
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                points.Add(newPoint);
            }
        }

        private static IList<OpenTK.Vector2d> _OrderCounterClockwise(IList<OpenTK.Vector2d> points)
        {
            double mX = 0;
            double mY = 0;
            foreach (OpenTK.Vector2d p in points)
            {
                mX += p.X;
                mY += p.Y;
            }
            mX /= points.Count;
            mY /= points.Count;

            return points.OrderBy(v => Math.Atan2(v.Y - mY, v.X - mX)).ToList();
        }

        public static OpenTK.Vector2d[] GetTwoPolygonsIntersectPoints(
            OpenTK.Vector2d[] polyPoints1, OpenTK.Vector2d[] polyPoints2)
        {
            IList<OpenTK.Vector2d> clippedPoints = new List<OpenTK.Vector2d>();

            // poly1の点でpoly2内部にある頂点を追加       
            for (int i = 0; i < polyPoints1.Length; i++)
            {
                if (IsPointInsidePolygon(polyPoints1[i], polyPoints2))
                {
                    _AddPoint(clippedPoints, polyPoints1[i]);
                }
            }

            // poly2の点でpoly1内部にある頂点を追加
            for (int i = 0; i < polyPoints2.Length; i++)
            {
                if (IsPointInsidePolygon(polyPoints2[i], polyPoints1))
                {
                    _AddPoint(clippedPoints, polyPoints2[i]);
                }
            }

            // poly1の線分とpoly2の交点を追加
            for (int i = 0, next = 1; i < polyPoints1.Length; i++, next = (i + 1 == polyPoints1.Length) ? 0 : i + 1)
            {
                OpenTK.Vector2d[] intesectPoints =
                    GetLineSegPolyIntersectPoints(polyPoints1[i], polyPoints1[next], polyPoints2);
                foreach (OpenTK.Vector2d intersectP in intesectPoints)
                {
                    _AddPoint(clippedPoints, intersectP);
                }
            }

            clippedPoints = _OrderCounterClockwise(clippedPoints);
            return clippedPoints.ToArray();
        }
    }
}
