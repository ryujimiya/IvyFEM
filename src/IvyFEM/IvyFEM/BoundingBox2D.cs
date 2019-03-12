using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class BoundingBox2D
    {
        public double MinX { get; set; } = 0;
        public double MaxX { get; set; } = 0;
        public double MinY { get; set; } = 0;
        public double MaxY { get; set; } = 0;
        public bool IsntEmpty { get; set; } = false;

        public BoundingBox2D()
        {

        }

        public BoundingBox2D(double xMin, double xMax, double yMin, double yMax)
        {
            MinX = xMin;
            MaxX = xMax;
            MinY = yMin;
            MaxY = yMax;
            System.Diagnostics.Debug.Assert(MinX <= MaxX);
            System.Diagnostics.Debug.Assert(MinY <= MaxY);
            IsntEmpty = true;
        }

        public BoundingBox2D(BoundingBox2D src)
        {
            MinX = src.MinX;
            MaxX = src.MaxX;
            MinY = src.MinY;
            MaxY = src.MaxY;
            IsntEmpty = src.IsntEmpty;
        }

        public void Copy(BoundingBox2D src)
        {
            MinX = src.MinX;
            MaxX = src.MaxX;
            MinY = src.MinY;
            MaxY = src.MaxY;
            IsntEmpty = src.IsntEmpty;
        }

        public string Dump()
        {
            string ret = "";
            string CRLF = System.Environment.NewLine;

            ret += "■BoundingBox2D" + CRLF;
            ret += "MinX = " + MinX + CRLF;
            ret += "MaxX = " + MaxX + CRLF;
            ret += "MinY = " + MinY + CRLF;
            ret += "MaxY = " + MaxY + CRLF;
            ret += "IsntEmpty = " + IsntEmpty + CRLF;
            return ret;
        }

        public static BoundingBox2D operator +(BoundingBox2D bb1, BoundingBox2D bb2)
        {
            BoundingBox2D result;

            if (!bb2.IsntEmpty)
            {
                result = new BoundingBox2D(bb1);
                return result;
            }
            if (!bb1.IsntEmpty)
            {
                result = new BoundingBox2D(bb2);
                return result;
            }

            double xMax = (bb1.MaxX > bb2.MaxX) ? bb1.MaxX : bb2.MaxX;
            double xMin = (bb1.MinX < bb2.MinX) ? bb1.MinX : bb2.MinX;
            double yMax = (bb1.MaxY > bb2.MaxY) ? bb1.MaxY : bb2.MaxY;
            double yMin = (bb1.MinY < bb2.MinY) ? bb1.MinY : bb2.MinY;
            result = new BoundingBox2D(xMin, xMax, yMin, yMax);
            return result;
        }

        public bool IsInside(OpenTK.Vector2d vec)
        {
            if (!IsntEmpty)
            {
                return false;
            }
            if (vec.X >= MinX && vec.X <= MaxX &&
                vec.Y >= MinY && vec.Y <= MaxY)
            {
                return true;
            }

            return false;
        }

        public bool IsIntersectSphere(OpenTK.Vector2d vec, double radius)
        {
            if (!IsntEmpty)
            {
                return false;
            }
            if (vec.X < MinX - radius || vec.X > MaxX + radius ||
               vec.Y < MinY - radius || vec.Y > MaxY + radius)
            {
                return false;
            }
            return true;
        }

        public bool IsIntersect(BoundingBox2D bb, double clearance)
        {
            if (bb.MinX > MaxX + clearance || bb.MaxX < MinX - clearance)
            {
                return false;
            }
            if (bb.MinY > MaxY + clearance || bb.MaxY < MinY - clearance)
            {
                return false;
            }
            return true;
        }
    }
}
