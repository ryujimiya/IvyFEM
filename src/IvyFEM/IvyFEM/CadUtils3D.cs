using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class CadUtils3D
    {
        public static OpenTK.Vector3d GetDirection(OpenTK.Vector3d v1, OpenTK.Vector3d v2)
        {
            OpenTK.Vector3d dir = v2 - v1;
            dir = OpenTK.Vector3d.Normalize(dir);
            return dir;
        }

        public static double Dot3D(double[] v1, double[] v2)
        {
            double ret = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
            return ret;
        }

        public static double[] Normalize3D(double[] v)
        {
            double len = Math.Sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            double[] ret = { v[0] / len, v[1] / len, v[2] / len };
            return ret;
        }

        public static double GetDistance3D(double[] v1, double[] v2)
        {
            double[] vec = { v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2] };
            double d = Math.Sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
            return d;
        }

        public static double[] GetDirection3D(double[] v1, double[] v2)
        {
            double[] dir = { v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2] };
            dir = Normalize3D(dir);
            return dir;
        }

		public static double[] TriNormal3D(double[] v1, double[] v2, double[] v3)
		{
			// Area vector
			double Ax = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]);
			double Ay = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]);
			double Az = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]);
			double[] n = { Ax, Ay, Az };
			n = Normalize3D(n);
			return n;
		}

		public static OpenTK.Vector2d ProjectToPlane(
			OpenTK.Vector3d p,
			OpenTK.Vector3d planeOrigin, OpenTK.Vector3d planeNormal, OpenTK.Vector3d planeXDir)
		{
			double x = OpenTK.Vector3d.Dot((p - planeOrigin), planeXDir);
			double y = OpenTK.Vector3d.Dot((p - planeOrigin), OpenTK.Vector3d.Cross(planeNormal, planeXDir));
			return new OpenTK.Vector2d(x, y);
		}

		public static OpenTK.Vector3d UnProjectFromPlane(
			OpenTK.Vector2d p,
			OpenTK.Vector3d planeOrigin, OpenTK.Vector3d planeNormal, OpenTK.Vector3d planeXDir)
		{
			return planeOrigin + planeXDir * p.X + OpenTK.Vector3d.Cross(planeNormal, planeXDir) * p.Y;
		}

		public static OpenTK.Vector3d GetVerticalUnitVector(OpenTK.Vector3d v)
        {
			double eps = Constants.PrecisionLowerLimit;
			double x;
			double y;
			double z;
			if (Math.Abs(v.X) >= eps)
			{
				if (Math.Abs(v.Y) >= eps)
                {
					y = 1.0;
					z = 0.0;
				}
				else
                {
					y = 0.0;
					z = 1.0;
				}
				x = (-v.Y * y - v.Z * z) / v.X;
			}
			else if (Math.Abs(v.Y) >= eps)
			{
				if (Math.Abs(v.X) >= eps)
				{
					x = 1.0;
					z = 0.0;
				}
				else
				{
					x = 0.0;
					z = 1.0;
				}
				y = (-v.X * x - v.Z * z) / v.Y;
			}
			else if (Math.Abs(v.Z) >= eps)
			{
				if (Math.Abs(v.X) >= eps)
				{
					x = 1.0;
					y = 0.0;
				}
				else
				{
					x = 0.0;
					y = 1.0;
				}
				z = (-v.X * x - v.Y * y) / v.Z;
			}
			else
            {
				// fail
				x = 0;
				y = 0;
				z = 0;
            }
			OpenTK.Vector3d normal = new OpenTK.Vector3d(x, y, z);
			normal.Normalize();
			return normal;
		}

		public static double TriArea(OpenTK.Vector3d v1, OpenTK.Vector3d v2, OpenTK.Vector3d v3)
        {
            // Area vector
            double Ax = (v2.Y - v1.Y) * (v3.Z - v1.Z) - (v3.Y - v1.Y) * (v2.Z - v1.Z);
            double Ay = (v2.Z - v1.Z) * (v3.X - v1.X) - (v3.Z - v1.Z) * (v2.X - v1.X);
            double Az = (v2.X - v1.X) * (v3.Y - v1.Y) - (v3.X - v1.X) * (v2.Y - v1.Y);
            return 0.5 * Math.Sqrt(Ax * Ax + Ay * Ay + Az * Az);
        }

		public static OpenTK.Vector3d TriNormal(OpenTK.Vector3d v1, OpenTK.Vector3d v2, OpenTK.Vector3d v3)
        {
            // Area vector
            double Ax = (v2.Y - v1.Y) * (v3.Z - v1.Z) - (v3.Y - v1.Y) * (v2.Z - v1.Z);
            double Ay = (v2.Z - v1.Z) * (v3.X - v1.X) - (v3.Z - v1.Z) * (v2.X - v1.X);
            double Az = (v2.X - v1.X) * (v3.Y - v1.Y) - (v3.X - v1.X) * (v2.Y - v1.Y);
            OpenTK.Vector3d n = new OpenTK.Vector3d(Ax, Ay, Az);
            n = OpenTK.Vector3d.Normalize(n);
            return n;
        }

		// 四面体の高さ
		public static double TetHeight(OpenTK.Vector3d v1, OpenTK.Vector3d v2, OpenTK.Vector3d v3, OpenTK.Vector3d v4)
        {
            // get normal vector
            double Ax = (v2.Y - v1.Y) * (v3.Z - v1.Z) - (v2.Z - v1.Z) * (v3.Y - v1.Y);
            double Ay = (v2.Z - v1.Z) * (v3.X - v1.X) - (v2.X - v1.X) * (v3.Z - v1.Z);
            double Az = (v2.X - v1.X) * (v3.Y - v1.Y) - (v2.Y - v1.Y) * (v3.X - v1.X);

            // normalize normal vector
            double invNorm = 1.0 / Math.Sqrt(Ax * Ax + Ay * Ay + Az * Az);
            Ax *= invNorm;
            Ay *= invNorm;
            Az *= invNorm;

            return (v4.X - v1.X) * Ax + (v4.Y - v1.Y) * Ay + (v4.Z - v1.Z) * Az;
        }

		// 四面体の体積
        public static double TetVolume(OpenTK.Vector3d v1, OpenTK.Vector3d v2, OpenTK.Vector3d v3, OpenTK.Vector3d v4)
        {
            double vol = (
				(v2.X - v1.X) * ((v3.Y - v1.Y) * (v4.Z - v1.Z) - (v4.Y - v1.Y) * (v3.Z - v1.Z)) -
				(v2.Y - v1.Y) * ((v3.X - v1.X) * (v4.Z - v1.Z) - (v4.X - v1.X) * (v3.Z - v1.Z)) +
				(v2.Z - v1.Z) * ((v3.X - v1.X) * (v4.Y - v1.Y) - (v4.X - v1.X) * (v3.Y - v1.Y))) / 6.0;
            return vol;
        }

		public static double SquareDistance(OpenTK.Vector3d pt0, OpenTK.Vector3d pt1)
		{
			var v = pt1 - pt0;
			double len = v.Length;
			return len * len;
		}

		public static double SquareLength(OpenTK.Vector3d point)
		{
			double len = point.Length;
			return len * len;
		}

		// 四面体の外接球の半径の２乗
		public static double SquareCircumradius(
			OpenTK.Vector3d pt0,
			OpenTK.Vector3d pt1,
			OpenTK.Vector3d pt2,
			OpenTK.Vector3d pt3)
		{
			double[][] _base = new double[3][] {
				new double[3] { pt1.X - pt0.X, pt1.Y - pt0.Y, pt1.Z - pt0.Z },
				new double[3] { pt2.X - pt0.X, pt2.Y - pt0.Y, pt2.Z - pt0.Z },
				new double[3] { pt3.X - pt0.X, pt3.Y - pt0.Y, pt3.Z - pt0.Z }
			};
			double[] s = new double[6] {
				_base[0][0] * _base[0][0] + _base[0][1] * _base[0][1] + _base[0][2] * _base[0][2],
				_base[1][0] * _base[1][0] + _base[1][1] * _base[1][1] + _base[1][2] * _base[1][2],
				_base[2][0] * _base[2][0] + _base[2][1] * _base[2][1] + _base[2][2] * _base[2][2],
				_base[1][0] * _base[2][0] + _base[1][1] * _base[2][1] + _base[1][2] * _base[2][2],
				_base[2][0] * _base[0][0] + _base[2][1] * _base[0][1] + _base[2][2] * _base[0][2],
				_base[0][0] * _base[1][0] + _base[0][1] * _base[1][1] + _base[0][2] * _base[1][2]
			};
			double vol = TetVolume(pt0, pt1, pt2, pt3) * 6.0;
			if (vol < 1.0e-20)
			{
				System.Diagnostics.Debug.Assert(false);
			}
			double invDet = 1.0 / (vol * vol);
			double[] t = new double[6] {
				(s[1] * s[2] - s[3] * s[3]) * 0.5 * invDet,
				(s[2] * s[0] - s[4] * s[4]) * 0.5 * invDet,
				(s[0] * s[1] - s[5] * s[5]) * 0.5 * invDet,
				(s[4] * s[5] - s[0] * s[3]) * 0.5 * invDet,
				(s[5] * s[3] - s[1] * s[4]) * 0.5 * invDet,
				(s[3] * s[4] - s[2] * s[5]) * 0.5 * invDet,
			};
			double[] u = new double[3] {
				t[0] * s[0] + t[5] * s[1] + t[4] * s[2],
				t[5] * s[0] + t[1] * s[1] + t[3] * s[2],
				t[4] * s[0] + t[3] * s[1] + t[2] * s[2],
			};
			return 0.5 * (u[0] * s[0] + u[1] * s[1] + u[2] * s[2]);
		}

		// 四面体の外接球の半径
		public static double Circumradius(
			OpenTK.Vector3d pt0,
			OpenTK.Vector3d pt1,
			OpenTK.Vector3d pt2,
			OpenTK.Vector3d pt3)
		{
			return Math.Sqrt(SquareCircumradius(pt0, pt1, pt2, pt3));
		}

		public static bool IsPointInsideTriangle(
			OpenTK.Vector3d p, OpenTK.Vector3d v1, OpenTK.Vector3d v2, OpenTK.Vector3d v3)
		{
			var vec1 = v3 - v1;
			var vec2 = v2 - v1;
			var vec3 = p - v1;

			double dot11 = OpenTK.Vector3d.Dot(vec1, vec1);
			double dot12 = OpenTK.Vector3d.Dot(vec1, vec2);
			double dot13 = OpenTK.Vector3d.Dot(vec1, vec3);
			double dot22 = OpenTK.Vector3d.Dot(vec2, vec2);
			double dot23 = OpenTK.Vector3d.Dot(vec2, vec3);

			// Compute barycentric coordinates
			double invDenom = 1 / (dot11 * dot22 - dot12 * dot12);
			double u = (dot22 * dot13 - dot12 * dot23) * invDenom;
			double v = (dot11 * dot23 - dot12 * dot13) * invDenom;

			// Check if point is in triangle
			return (u >= 0) && (v >= 0) && (u + v < 1);
		}

		private static bool IsPointSameSideOfFace(
			OpenTK.Vector3d p,
			OpenTK.Vector3d v1, OpenTK.Vector3d v2, OpenTK.Vector3d v3, OpenTK.Vector3d v4)
		{
			double eps = 1.0e-12;
			OpenTK.Vector3d normal = OpenTK.Vector3d.Cross(v2 - v1, v3 - v1);
			double dotV4 = OpenTK.Vector3d.Dot(normal, v4 - v1);
			double dotP = OpenTK.Vector3d.Dot(normal, p - v1);
			return (Math.Sign(dotV4) == Math.Sign(dotP) ||
					Math.Abs(dotP) < eps);
		}

		public static bool IsPointInsideTetrahedron(
			OpenTK.Vector3d p,
			OpenTK.Vector3d v1, OpenTK.Vector3d v2, OpenTK.Vector3d v3, OpenTK.Vector3d v4)
		{
			return IsPointSameSideOfFace(p, v1, v2, v3, v4) &&
				   IsPointSameSideOfFace(p, v2, v3, v4, v1) &&
				   IsPointSameSideOfFace(p, v3, v4, v1, v2) &&
				   IsPointSameSideOfFace(p, v4, v1, v2, v3);
		}

		public static bool GetLineTriangleIntersectPoint(
			OpenTK.Vector3d lineOrigin, OpenTK.Vector3d lineDir,
			OpenTK.Vector3d v1, OpenTK.Vector3d v2, OpenTK.Vector3d v3,
			out OpenTK.Vector3d intersectPoint)
		{
			intersectPoint = new OpenTK.Vector3d();
			double eps1 = 1.0e-12;
			double eps2 = 1.0e-12;
			OpenTK.Vector3d edge1;
			OpenTK.Vector3d edge2;
			OpenTK.Vector3d h;
			OpenTK.Vector3d s;
			OpenTK.Vector3d q;
			double a;
			double f;
			double u;
			double v;
			edge1 = v2 - v1;
			edge2 = v3 - v1;
			h = OpenTK.Vector3d.Cross(lineDir, edge2);
			a = OpenTK.Vector3d.Dot(edge1, h);
			if (a > -eps1 && a < eps1)
            {
				// 直線は三角形と平行
				return false;
			}
			f = 1.0 / a;
			s = lineOrigin - v1;
			u = f * OpenTK.Vector3d.Dot(s, h);
			if (u < 0.0 -  eps2 || u > 1.0 + eps2)
            {
				return false;
			}
			q = OpenTK.Vector3d.Cross(s, edge1);
			v = f * OpenTK.Vector3d.Dot(lineDir, q);
			if (v < 0.0 - eps2 || u + v > 1.0 + eps2)
            {
				return false;
			}

			double t = f * OpenTK.Vector3d.Dot(edge2, q);
			intersectPoint = lineOrigin + lineDir * t;
			return true;
		}

		public static bool GetLinePlaneIntersectPoint(
			OpenTK.Vector3d lineOrigin, OpenTK.Vector3d lineDir,
			OpenTK.Vector3d planeOrigin, OpenTK.Vector3d planeNormal,
			out OpenTK.Vector3d intersectPoint)
		{
			intersectPoint = new OpenTK.Vector3d();
			double eps = 1.0e-12;
			double a = OpenTK.Vector3d.Dot(planeNormal, lineDir);
			if (Math.Abs(a) < eps)
            {
				// 直線と平面は平行
				return false;
            }
			double b = OpenTK.Vector3d.Dot(planeNormal, planeOrigin - lineOrigin);
			double t = b / a;
			intersectPoint = lineOrigin + lineDir * t;
			return true;
		}
	}
}
