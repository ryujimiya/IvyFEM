using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Camera3D : Camera
    {
        public OpenTK.Quaterniond RotQuat { get; set; } = new OpenTK.Quaterniond(0.0, 0.0, 0.0, 1.0);

        public override void MouseRotation(double movBeginX, double movBeginY, double movEndX, double movEndY)
        {
            if (Math.Abs(movBeginX - movEndX) < 1.0e-8 &&
                Math.Abs(movBeginY - movEndY) < 1.0e-8)
            {
                return;
            }
            OpenTK.Vector3d p1;
            {
                double norm2 = movBeginX * movBeginX + movBeginY * movBeginY;
                double z;
                if (norm2 < 0.5)
                {
                    z = Math.Sqrt(1.0 - norm2);
                }
                else
                {
                    z = 0.343443391 / norm2;
                }
                p1 = new OpenTK.Vector3d(movBeginX, movBeginY, z);
            }
            OpenTK.Vector3d p2;
            {
                double norm2 = movEndX * movEndX + movEndY * movEndY;
                double z;
                if (norm2 < 0.5)
                { 
                    z = Math.Sqrt(1.0 - norm2);
                }
                else
                {
                    z = 0.343443391 / norm2;
                }
                p2 = new OpenTK.Vector3d(movEndX, movEndY, z);
            }
            OpenTK.Vector3d norm = OpenTK.Vector3d.Cross(p1, p2);
            norm = OpenTK.Vector3d.Normalize(norm);
            p2 -= p1;
            double t = p2.Length / 2.0;
            if (t > 1.0)
            { 
                t = 1.0;
            }
            else if (t < -1.0)
            {
                t = -1.0;
            }
            double phi = 5.0 * Math.Asin(t);
            RotQuat = OpenTK.Quaterniond.FromAxisAngle(norm, phi) * RotQuat;
            RotQuat = OpenTK.Quaterniond.Normalize(RotQuat);
        }

        public override OpenTK.Matrix3d RotMatrix33()
        {
            OpenTK.Matrix3d m = new OpenTK.Matrix3d();
            double vx = RotQuat.X;
            double vy = RotQuat.Y;
            double vz = RotQuat.Z;
            double real = RotQuat.W;

            m[0, 0] = 1.0 - 2.0 * (vy * vy + vz * vz);
            m[0, 1] = 2.0 * (vx * vy - vz * real);
            m[0, 2] = 2.0 * (vz * vx + vy * real);

            m[1, 0] = 2.0 * (vx * vy + vz * real);
            m[1, 1] = 1.0 - 2.0 * (vz * vz + vx * vx);
            m[1, 2] = 2.0 * (vy * vz - vx * real);

            m[2, 0] = 2.0 * (vz * vx - vy * real);
            m[2, 1] = 2.0 * (vy * vz + vx * real);
            m[2, 2] = 1.0 - 2.0 * (vy * vy + vx * vx);
            return m;
        }
    }
}
