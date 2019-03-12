using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Camera2D : Camera
    {
        public double Theta { get; set; } = 0;

        public override void MouseRotation(double movBeginX, double movBeginY, double movEndX, double movEndY)
        {
            {
                double diffX = movEndX - movBeginX;
                double diffY = movEndY - movBeginY;
                double qLenDiff = diffX * diffX + diffY * diffY;
                if (Math.Abs(qLenDiff) < 1.0e-16)
                {
                    return;
                }
            }
            double centerX = WindowCenter.X;
            double centerY = WindowCenter.Y;
            double qLenEnd = (movEndX - centerX) * (movEndX - centerX) +
                (movEndY - centerY) * (movEndY - centerY);
            double qLenBegin = (movBeginX - centerX) * (movBeginX - centerX) + 
                (movBeginY - centerY) * (movBeginY - centerY);
            double dSin = ((movEndX - centerX) * (movBeginY - centerY) -
                (movEndY - centerY) * (movBeginX - centerX)) / Math.Sqrt(qLenEnd * qLenBegin);
            double dTheta = Math.Asin(dSin);
            Theta -= dTheta;
        }

        public override OpenTK.Matrix3d RotMatrix33()
        {
            var rot = new OpenTK.Matrix3d();
            double ct = Math.Cos(Theta);
            double st = Math.Sin(Theta);

            rot[0, 0] = ct;
            rot[1, 0] = st;
            rot[2, 0] = 0.0;

            rot[0, 1] = -st;
            rot[1, 1] = ct;
            rot[2, 1] = 0.0;

            rot[0, 2] = 0.0;
            rot[1, 2] = 0.0;
            rot[2, 2] = 1.0;
            return rot;
        }

        public OpenTK.Vector3d ProjectionOnPlane(double posX,   double posY,
            double planeX = 0, double planeY = 0, double planeZ = 0,
            double normX = 0,  double normY = 0,  double normZ = 1)
      	{
            System.Diagnostics.Debug.Assert(Math.Abs(planeX) < 1.0e-10);
            System.Diagnostics.Debug.Assert(Math.Abs(planeY) < 1.0e-10);
            System.Diagnostics.Debug.Assert(Math.Abs(planeZ) < 1.0e-10);

            System.Diagnostics.Debug.Assert(Math.Abs(normX) < 1.0e-10);
            System.Diagnostics.Debug.Assert(Math.Abs(normY) < 1.0e-10);
            System.Diagnostics.Debug.Assert(Math.Abs(normZ - 1.0) < 1.0e-10);

            double hw = HalfViewHeight * WindowAspect;
            double hh = HalfViewHeight;

            double objCX = ObjectCenter.X;
            double objCY = ObjectCenter.Y;
            double objCZ = ObjectCenter.Z;

            double ox = hw * Math.Cos(Theta) * (posX * InvScale - WindowCenter.X) +
                hh * Math.Sin(Theta) * (posY * InvScale - WindowCenter.Y) + objCX;
            double oy = -hw * Math.Sin(Theta) * (posX * InvScale - WindowCenter.X) +
                hh * Math.Cos(Theta) * (posY * InvScale - WindowCenter.Y) + objCY;
            double oz = +objCZ;

            return new OpenTK.Vector3d(ox, oy, oz);
        }

    }
}
