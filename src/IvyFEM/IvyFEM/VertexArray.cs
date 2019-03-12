using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class VertexArray
    {
        public double[] VertexCoordArray { get; set; } = null;
        public double[] UVCoordArray { get; set; } = null;
        public uint PointCount { get; private set; } = 0;
        public uint Dimension { get; private set; } = 0;

        public VertexArray()
        {

        }

        public VertexArray(uint ptCnt, uint dim)
        {
            PointCount = ptCnt;
            Dimension = dim;
            VertexCoordArray = new double[PointCount * Dimension];
            UVCoordArray = null;
        }

        public void SetSize(uint ptCnt, uint dim)
        {
            if (PointCount == ptCnt && Dimension == dim)
            {
                return;
            }
            PointCount = ptCnt;
            Dimension = dim;
            VertexCoordArray = new double[ptCnt * dim];
            UVCoordArray = new double[ptCnt * 2];
        }

        public BoundingBox3D GetBoundingBox(OpenTK.Matrix3d rot)
        {
            if (VertexCoordArray == null || VertexCoordArray.Length == 0)
            {
                return new BoundingBox3D();
            }
            if (rot == null)
            {
                if (Dimension == 2)
                {
                    BoundingBox3D bb;
                    {
                        double x1 = VertexCoordArray[0];
                        double y1 = VertexCoordArray[1];
                        double z1 = 0.0;
                        bb = new BoundingBox3D(x1, x1, y1, y1, z1, z1);
                    }
                    for (uint iPt = 1; iPt < PointCount; iPt++)
                    {
                        double x1 = VertexCoordArray[iPt * 2];
                        double y1 = VertexCoordArray[iPt * 2 + 1];
                        double z1 = 0.0;
                        bb.MaxX = (x1 > bb.MaxX) ? x1 : bb.MaxX; bb.MinX = (x1 < bb.MinX) ? x1 : bb.MinX;
                        bb.MaxY = (y1 > bb.MaxY) ? y1 : bb.MaxY; bb.MinY = (y1 < bb.MinY) ? y1 : bb.MinY;
                        bb.MaxZ = (z1 > bb.MaxZ) ? z1 : bb.MaxZ; bb.MinZ = (z1 < bb.MinZ) ? z1 : bb.MinZ;
                    }
                    return bb;
                }
                if (Dimension == 3)
                {
                    BoundingBox3D bb;
                    {
                        double x1 = VertexCoordArray[0];
                        double y1 = VertexCoordArray[1];
                        double z1 = 0.0;
                        bb = new BoundingBox3D(x1, x1, y1, y1, z1, z1);
                    }
                    for (uint iPt = 1; iPt < PointCount; iPt++)
                    {
                        double x1 = VertexCoordArray[iPt * 3];
                        double y1 = VertexCoordArray[iPt * 3 + 1];
                        double z1 = VertexCoordArray[iPt * 3 + 2];
                        bb.MaxX = (x1 > bb.MaxX) ? x1 : bb.MaxX; bb.MinX = (x1 < bb.MinX) ? x1 : bb.MinX;
                        bb.MaxY = (y1 > bb.MaxY) ? y1 : bb.MaxY; bb.MinY = (y1 < bb.MinY) ? y1 : bb.MinY;
                        bb.MaxZ = (z1 > bb.MaxZ) ? z1 : bb.MaxZ; bb.MinZ = (z1 < bb.MinZ) ? z1 : bb.MinZ;
                    }
                    return bb;
                }
            }
            if (Dimension == 2)
            {
                double minX;
                double maxX;
                double minY;
                double maxY;
                double minZ;
                double maxZ;
                {
                    double x1 = VertexCoordArray[0];
                    double y1 = VertexCoordArray[1];
                    double z1 = 0.0;
                    minX = maxX = x1 * rot[0, 0] + y1 * rot[0, 1] + z1 * rot[0, 2];
                    minY = maxY = x1 * rot[1, 0] + y1 * rot[1, 1] + z1 * rot[1, 2];
                    minZ = maxZ = x1 * rot[2, 0] + y1 * rot[2, 1] + z1 * rot[2, 2];
                }
                for (uint iPt = 1; iPt < PointCount; iPt++)
                {
                    double x1 = VertexCoordArray[iPt * 2];
                    double y1 = VertexCoordArray[iPt * 2 + 1];
                    double z1 = 0.0;
                    double x2 = x1 * rot[0, 0] + y1 * rot[0, 1] + z1 * rot[0, 2];
                    double y2 = x1 * rot[1, 0] + y1 * rot[1, 1] + z1 * rot[1, 2];
                    double z2 = x1 * rot[2, 0] + y1 * rot[2, 1] + z1 * rot[2, 2];
                    maxX = (x2 > maxX) ? x2 : maxX; minX = (x2 < minX) ? x2 : minX;
                    maxY = (y2 > maxY) ? y2 : maxY; minY = (y2 < minY) ? y2 : minY;
                    maxZ = (z2 > maxZ) ? z2 : maxZ; minZ = (z2 < minZ) ? z2 : minZ;
                }

                double c1X = (minX + maxX) * 0.5;
                double c1Y = (minY + maxY) * 0.5;
                double c1Z = (minZ + maxZ) * 0.5;
                double c2X = c1X * rot[0, 0] + c1Y * rot[1, 0] + c1Z * rot[2, 0];
                double c2Y = c1X * rot[0, 1] + c1Y * rot[1, 1] + c1Z * rot[2, 1];
                double c2Z = c1X * rot[0, 2] + c1Y * rot[1, 2] + c1Z * rot[2, 2];
                double hX = (maxX - minX) * 0.5;
                double hY = (maxY - minY) * 0.5;
                double hZ = (maxZ - minZ) * 0.5;
                BoundingBox3D bb = new BoundingBox3D(
                    c2X - hX, c2X + hX, 
                    c2Y - hY, c2Y + hY, 
                    c2Z - hZ, c2Z + hZ);
                return bb;
            }
            if (Dimension == 3) // view axis alligned bounding box
            {
                double minX;
                double maxX;
                double minY;
                double maxY;
                double minZ;
                double maxZ;
                {
                    double x1 = VertexCoordArray[0];
                    double y1 = VertexCoordArray[1];
                    double z1 = VertexCoordArray[2];
                    minX = maxX = x1 * rot[0, 0] + y1 * rot[0, 1] + z1 * rot[0, 2];
                    minY = maxY = x1 * rot[1, 0] + y1 * rot[1, 1] + z1 * rot[1, 2];
                    minZ = maxZ = x1 * rot[2, 0] + y1 * rot[2, 1] + z1 * rot[2, 2];
                }
                for (uint iPt = 1; iPt < PointCount; iPt++)
                {
                    double x1 = VertexCoordArray[iPt * 3];
                    double y1 = VertexCoordArray[iPt * 3 + 1];
                    double z1 = VertexCoordArray[iPt * 3 + 2];
                    double x2 = x1 * rot[0, 0] + y1 * rot[0, 1] + z1 * rot[0, 2];
                    double y2 = x1 * rot[1, 0] + y1 * rot[1, 1] + z1 * rot[1, 2];
                    double z2 = x1 * rot[2, 0] + y1 * rot[2, 1] + z1 * rot[2, 2];
                    maxX = (x2 > maxX) ? x2 : maxX; minX = (x2 < minX) ? x2 : minX;
                    maxY = (y2 > maxY) ? y2 : maxY; minY = (y2 < minY) ? y2 : minY;
                    maxZ = (z2 > maxZ) ? z2 : maxZ; minZ = (z2 < minZ) ? z2 : minZ;
                }

                double c1X = (minX + maxX) * 0.5;
                double c1Y = (minY + maxY) * 0.5;
                double c1Z = (minZ + maxZ) * 0.5;
                double c2X = c1X * rot[0, 0] + c1Y * rot[1, 0] + c1Z * rot[2, 0];
                double c2Y = c1X * rot[0, 1] + c1Y * rot[1, 1] + c1Z * rot[2, 1];
                double c2Z = c1X * rot[0, 2] + c1Y * rot[1, 2] + c1Z * rot[2, 2];
                double hX = (maxX - minX) * 0.5;
                double hY = (maxY - minY) * 0.5;
                double hZ = (maxZ - minZ) * 0.5;
                BoundingBox3D bb = new BoundingBox3D(
                    c2X - hX, c2X + hX,
                    c2Y - hY, c2Y + hY,
                    c2Z - hZ, c2Z + hZ);
                return bb;
            }

            return new BoundingBox3D();
        }

        public void EnableUVMap(bool isUVMap)
        {
            if ((UVCoordArray != null) == isUVMap)
            {
                return;
            }
            if (isUVMap)
            {
                uint ptCnt = PointCount;
                UVCoordArray = new double[ptCnt * 2];
                for (uint i = 0; i < ptCnt * 2; i++)
                {
                    UVCoordArray[i] = 0;
                }
            }
            else
            {
                UVCoordArray = null;
                UVCoordArray = null;
            }
        }

    }
}
