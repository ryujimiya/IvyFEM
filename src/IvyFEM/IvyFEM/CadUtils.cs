using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class CadUtils
    {
        public static double Dot(double[] v1, double[] v2)
        {
            double ret;
            if (v1.Length == 2)
            {
                ret = CadUtils2D.Dot2D(v1, v2);
            }
            else if (v1.Length == 3)
            {
                ret = CadUtils3D.Dot3D(v1, v2);
            }
            else
            {
                throw new NotImplementedException();
            }
            return ret;
        }

        public static double[] Normalize(double[] v)
        {
            double[] ret;
            if (v.Length == 2)
            {
                ret = CadUtils2D.Normalize2D(v);
            }
            else if (v.Length == 3)
            {
                ret = CadUtils3D.Normalize3D(v);
            }
            else
            {
                throw new NotImplementedException();
            }
            return ret;
        }

        public static double GetDistance(double[] v1, double[] v2)
        {
            double ret;
            if (v1.Length == 2)
            {
                ret = CadUtils2D.GetDistance2D(v1, v2);
            }
            else if (v1.Length == 3)
            {
                ret = CadUtils3D.GetDistance3D(v1, v2);
            }
            else
            {
                throw new NotImplementedException();
            }
            return ret;
        }

        public static double[] GetDirection(double[] v1, double[] v2)
        {
            double[] ret;
            if (v1.Length == 2)
            {
                ret = CadUtils2D.GetDirection2D(v1, v2);
            }
            else if (v1.Length == 3)
            {
                ret = CadUtils3D.GetDirection3D(v1, v2);
            }
            else
            {
                throw new NotImplementedException();
            }
            return ret;
        }

        public static double TriArea(double[] v1, double[] v2, double[] v3)
        {
            double ret;
            if (v1.Length == 2)
            {
                ret = CadUtils2D.TriArea(
                    new OpenTK.Vector2d(v1[0], v1[1]),
                    new OpenTK.Vector2d(v2[0], v2[1]),
                    new OpenTK.Vector2d(v3[0], v3[1]));
            }
            else if (v1.Length == 3)
            {
                ret = CadUtils3D.TriArea(
                    new OpenTK.Vector3d(v1[0], v1[1], v1[2]),
                    new OpenTK.Vector3d(v2[0], v2[1], v2[2]),
                    new OpenTK.Vector3d(v3[0], v3[1], v3[2]));
            }
            else
            {
                throw new NotImplementedException();
            }
            return ret;
        }
    }
}
