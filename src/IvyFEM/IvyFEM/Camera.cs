using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract class Camera
    {
        public double WindowAspect { get; set; }
        public OpenTK.Vector2d WindowCenter { get; set; } = new OpenTK.Vector2d();
        protected double InvScale;
        public double Scale { get => GetScale(); set => SetScale(value); }
        public double HalfViewHeight { get; private set; }

        public double MinLen { get; set; } = 1.0e-4;
        private double ObjectW;
        private double ObjectH;
        private double ObjectD;
        public OpenTK.Vector3d ObjectCenter { get; set; } = new OpenTK.Vector3d();

        private bool isPers;
        public bool IsPers
        {
            get
            {
                return isPers;
            }
            set
            {
                SetIsPers(value);
            }
        }
        private double fovY;
        public double FovY
        {
            get
            {
                return GetFovY();
            }
            set
            {
                SetFovY(value);
            }
        }
        private double Dist;

        public Camera()
        {
            WindowAspect = 1.0;

            HalfViewHeight = 1.0;
            InvScale = 1.0;
            WindowCenter = new OpenTK.Vector2d(0.0, 0.0);

            ObjectW = 1.0;
            ObjectH = 1.0;
            ObjectD = 1.0;
            ObjectCenter = OpenTK.Vector3d.Zero;

            IsPers = false;
            fovY = 30.0 * Math.PI / 180.0;
            Dist = HalfViewHeight / Math.Tan(fovY * 0.5) + ObjectD * 0.5;
        }

        private void SetIsPers(bool value)
        {
            isPers = value;
            if (IsPers)
            {
                Dist = HalfViewHeight * InvScale / Math.Tan(fovY * 0.5) + ObjectD * 0.5;
            }
        }

        private double GetScale()
        {
            return 1.0 / InvScale;
        }

        private void SetScale(double value)
        {
            if (value < 0.01)
            {
                InvScale = 100;
            }
            else if (value > 100.0)
            {
                InvScale = 0.01;
            }
            else
            {
                InvScale = 1.0 / value;
            }
            Dist = HalfViewHeight * InvScale / Math.Tan(fovY * 0.5) + ObjectD * 0.5;
        }

        private double GetFovY()
        {
            return fovY * 180.0 / Math.PI;
        }

        private void SetFovY(double fovY)
        {
            if (!IsPers)
            {
                return;
            }
            double _fovY = fovY;
            if (fovY < 15.0)
            {
                _fovY = 15.0;
            }
            else if (fovY > 90.0)
            {
                _fovY = 90.0;
            }
            else
            {
                _fovY = fovY;
            }
            this.fovY = _fovY * Math.PI / 180.0;
            Dist = HalfViewHeight * InvScale / Math.Tan(this.fovY * 0.5) + ObjectD * 0.5;
        }

        public void GetPerspective(out double fovY, out double aspect,
            out double clipNear, out double clipFar)
        {
            fovY = this.fovY * 180.0 / Math.PI;
            aspect = WindowAspect;
            clipNear = 0.001;
            clipFar = Dist * 2.0 + ObjectD * 20.0 + 100.0;
        }

        public void GetOrtho(out double hw, out double hh, out double hd)
        {
            hw = HalfViewHeight * InvScale * WindowAspect;
            hh = HalfViewHeight * InvScale;
            hd = (ObjectD * 20.0 > MinLen) ? ObjectD * 20.0 : MinLen;
            hd = (ObjectW * 20.0 > hd) ? ObjectW * 20.0 : hd;
            hd = (ObjectH * 20.0 > hd) ? ObjectH * 20.0 : hd;
        }

        public void GetCenterPosition(out double x, out double y, out double z)
        {
            if (IsPers)
            {
                System.Diagnostics.Debug.Assert(false);
                x = 0.0;
                y = 0.0;
                z = -Dist;
            }
            else
            {
                x = HalfViewHeight * WindowCenter.X * WindowAspect;
                y = HalfViewHeight * WindowCenter.Y;
                z = 0.0;
            }
        }

        public OpenTK.Vector3d GetCenterPosition()
        {
            double x;
            double y;
            double z;
            GetCenterPosition(out x, out y, out z);
            return new OpenTK.Vector3d(x, y, z);
        }

        public void GetObjectCenter(out double x, out double y, out double z)
    	{
            x = ObjectCenter.X;
            y = ObjectCenter.Y;
            z = ObjectCenter.Z;
        }

        public OpenTK.Vector3d GetObjectCenter()
        {
            return ObjectCenter;
        }

        public void SetObjectCenter(double x, double y, double z)
        {
            ObjectCenter = new OpenTK.Vector3d(x, y, z);
        }

        public void GetObjectSize(out double w, out double h, out double d)
        {
            w = ObjectW;
            h = ObjectH;
            d = ObjectD;
        }

        public void SetObjectSize(double w, double h, double d)
        {
            ObjectW = w;
            ObjectH = h;
            ObjectD = d;
        }

        public void SetObjectBoundingBox(BoundingBox3D bb)
        {
            ObjectCenter = new OpenTK.Vector3d(
                ((bb.MinX + bb.MaxX) * 0.5),
                ((bb.MinY + bb.MaxY) * 0.5),
                ((bb.MinZ + bb.MaxZ) * 0.5)
                );
            ObjectW = bb.MaxX - bb.MinX;
            ObjectH = bb.MaxY - bb.MinY;
            ObjectD = bb.MaxZ - bb.MinZ;
        }

        public void Fit()
        {
            double margin = 1.5;
            double objAspect = ObjectW / ObjectH;
            InvScale = 1.0;
            if (objAspect < WindowAspect)
            {
                HalfViewHeight = ObjectH * 0.5 * margin;
            }
            else
            {
                double tmpH = ObjectW / WindowAspect;
                HalfViewHeight = tmpH * 0.5 * margin;
            }
            Dist = HalfViewHeight * InvScale / Math.Tan(fovY * 0.5) + ObjectD * 0.5;
            WindowCenter = new OpenTK.Vector2d(0.0, 0.0);
        }

        public void Fit(BoundingBox3D bb)
        {
            SetObjectBoundingBox(bb);
            Fit();
        }

        public void MousePan(double movBeginX, double movBeginY, double movEndX, double movEndY)
        {
            double x = WindowCenter.X;
            double y = WindowCenter.Y;
            
            x += (movEndX - movBeginX) * InvScale;
            y += (movEndY - movBeginY) * InvScale;

            WindowCenter = new OpenTK.Vector2d(x, y);
        }

        public abstract void MouseRotation(double movBeginX, double movBeginY, double movEndX, double movEndY);

        public abstract OpenTK.Matrix3d RotMatrix33();

        public OpenTK.Matrix4d RotMatrix44()
        {
            OpenTK.Matrix4d rot = new OpenTK.Matrix4d();

            var rot1 = RotMatrix33();

            rot[0, 0] = rot1[0, 0];
            rot[1, 0] = rot1[1, 0];
            rot[2, 0] = rot1[2, 0];
            rot[3, 0] = 0.0;

            rot[0, 1] = rot1[0, 1];
            rot[1, 1] = rot1[1, 1];
            rot[2, 1] = rot1[2, 1];
            rot[3, 1] = 0.0;

            rot[0, 2] = rot1[0, 2];
            rot[1, 2] = rot1[1, 2];
            rot[2, 2] = rot1[2, 2];
            rot[3, 2] = 0.0;

            rot[0, 3] = 0.0;
            rot[1, 3] = 0.0;
            rot[2, 3] = 0.0;
            rot[3, 3] = 1.0;

            return rot;
        }

    }
}
