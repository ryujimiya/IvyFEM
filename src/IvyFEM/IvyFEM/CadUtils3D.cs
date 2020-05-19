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
    }
}
