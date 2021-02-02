using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class BoundingBox3D
    {
        public double MinX { get; set; } = 0;
        public double MaxX { get; set; } = 0;
        public double MinY { get; set; } = 0;
        public double MaxY { get; set; } = 0;
        public double MinZ { get; set; } = 0;
        public double MaxZ { get; set; } = 0;
        public bool IsntEmpty { get; set; } = false;

        public BoundingBox3D()
        {

        }

        public BoundingBox3D(
            double minX, double maxX, 
            double minY, double maxY, 
            double minZ, double maxZ)
        {
            MinX = minX;
            MaxX = maxX;
            MinY = minY;
            MaxY = maxY;
            MinZ = minZ;
            MaxZ = maxZ;

            System.Diagnostics.Debug.Assert(MinX <= MaxX);
            System.Diagnostics.Debug.Assert(MinY <= MaxY);
            System.Diagnostics.Debug.Assert(MinZ <= MaxZ);
            IsntEmpty = true;
        }

        public BoundingBox3D(BoundingBox3D src)
        {
            MinX = src.MinX;
            MaxX = src.MaxX;
            MinY = src.MinY;
            MaxY = src.MaxY;
            MinZ = src.MinZ;
            MaxZ = src.MaxZ;
            IsntEmpty = src.IsntEmpty;
        }

        public static BoundingBox3D operator +(BoundingBox3D bb1, BoundingBox3D bb2)
        {
            BoundingBox3D result = new BoundingBox3D(bb1);
            if (!bb2.IsntEmpty)
            {
                return result;
            }
            if (!bb1.IsntEmpty)
            {
                result.MaxX = bb2.MaxX;
                result.MinX = bb2.MinX;
                result.MaxY = bb2.MaxY;
                result.MinY = bb2.MinY;
                result.MaxZ = bb2.MaxZ;
                result.MinZ = bb2.MinZ;
                result.IsntEmpty = bb2.IsntEmpty;
                return result;
            }
            result.MaxX = (bb1.MaxX > bb2.MaxX) ? bb1.MaxX : bb2.MaxX;
            result.MinX = (bb1.MinX < bb2.MinX) ? bb1.MinX : bb2.MinX;
            result.MaxY = (bb1.MaxY > bb2.MaxY) ? bb1.MaxY : bb2.MaxY;
            result.MinY = (bb1.MinY < bb2.MinY) ? bb1.MinY : bb2.MinY;
            result.MaxZ = (bb1.MaxZ > bb2.MaxZ) ? bb1.MaxZ : bb2.MaxZ;
            result.MinZ = (bb1.MinZ < bb2.MinZ) ? bb1.MinZ : bb2.MinZ;
            return result;

        }

        public bool IsInside(OpenTK.Vector3d vec)
        {
            if (!IsntEmpty)
            {
                return false;
            }
            if (vec.X >= MinX && vec.X <= MaxX &&
                vec.Y >= MinY && vec.Y <= MaxY &&
                vec.Z >= MinZ && vec.Z <= MaxZ)
            {
                return true;
            }
            return false;
        }

        public bool IsPossibilityIntersectSphere(OpenTK.Vector3d vec, double radius)
        {
            if (!IsntEmpty)
            {
                return false;
            }
            if (vec.X < MinX - radius || vec.X > MaxX + radius ||
                vec.Y < MinY - radius || vec.Y > MaxY + radius ||
                vec.Z < MinZ - radius || vec.Z > MaxZ + radius)
            {
                return false;
            }
            // 干渉しないやつが含まれているが、アバウトなままでよい
            return true;
        }

        public bool AddPoint(OpenTK.Vector3d vec, double eps)
        {
            if (eps <= 0)
            {
                return false;
            }
            if (IsntEmpty)
            {
                MinX = (MinX < vec.X - eps) ? MinX : vec.X - eps;
                MinY = (MinY < vec.Y - eps) ? MinY : vec.Y - eps;
                MinZ = (MinZ < vec.Z - eps) ? MinZ : vec.Z - eps;
                MaxX = (MaxX > vec.X + eps) ? MaxX : vec.X + eps;
                MaxY = (MaxY > vec.Y + eps) ? MaxY : vec.Y + eps;
                MaxZ = (MaxZ > vec.Z + eps) ? MaxZ : vec.Z + eps;
            }
            else
            {
                IsntEmpty = true;
                MinX = vec.X - eps;
                MinY = vec.Y - eps;
                MinZ = vec.Z - eps;
                MaxX = vec.X + eps;
                MaxY = vec.Y + eps;
                MaxZ = vec.Z + eps;
            }
            return true;
        }

        public void SetValueToArray(double[] bb)
        {
            bb[0] = MinX;
            bb[2] = MinY;
            bb[4] = MinZ;
            bb[1] = MaxX;
            bb[3] = MaxY;
            bb[5] = MaxZ;
        }

        public void ProjectOnLine(out double minR, out double maxR, OpenTK.Vector3d org, OpenTK.Vector3d dir)
        {
            double[] d = new double[8]
            {
                OpenTK.Vector3d.Dot(dir, new OpenTK.Vector3d(MinX, MinY, MinZ) - org),
                OpenTK.Vector3d.Dot(dir, new OpenTK.Vector3d(MaxX, MinY, MinZ) - org),
                OpenTK.Vector3d.Dot(dir, new OpenTK.Vector3d(MinX, MaxY, MinZ) - org),
                OpenTK.Vector3d.Dot(dir, new OpenTK.Vector3d(MaxX, MaxY, MinZ) - org),
                OpenTK.Vector3d.Dot(dir, new OpenTK.Vector3d(MinX, MinY, MaxZ) - org),
                OpenTK.Vector3d.Dot(dir, new OpenTK.Vector3d(MaxX, MinY, MaxZ) - org),
                OpenTK.Vector3d.Dot(dir, new OpenTK.Vector3d(MinX, MaxY, MaxZ) - org),
                OpenTK.Vector3d.Dot(dir, new OpenTK.Vector3d(MaxX, MaxY, MaxZ) - org)
            };
            minR = d[0];
            maxR = d[0];
            for (uint i = 1; i < 8; i++)
            {
                minR = (d[i] < minR) ? d[i] : minR;
                maxR = (d[i] > maxR) ? d[i] : maxR;
            }
        }

    }
}
