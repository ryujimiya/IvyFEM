using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Edge3D : IObject
    {
        public uint SVId { get; set; } = 0;
        public uint EVId { get; set; } = 0;
        public OpenTK.Vector3d SPt { get; set; } = new OpenTK.Vector3d();
        public OpenTK.Vector3d EPt { get; set; } = new OpenTK.Vector3d();
        public double[] Color { get; set; } = new double[3];

        public Edge3D()
        {
            Color = new double[3] { 0.0, 0.0, 0.0 };
        }

        public Edge3D(uint sVId, uint eVId)
        {
            SVId = sVId;
            EVId = eVId;
            Color = new double[3] { 0.0, 0.0, 0.0 };
        }

        public Edge3D(Edge3D src)
        {
            Copy(src);
        }

        public void Copy(IObject src)
        {
            Edge3D srcEdge = src as Edge3D;

            SVId = srcEdge.SVId;
            EVId = srcEdge.EVId;
            SPt = new OpenTK.Vector3d(srcEdge.SPt.X, srcEdge.SPt.Y, srcEdge.SPt.Z);
            EPt = new OpenTK.Vector3d(srcEdge.EPt.X, srcEdge.EPt.Y, srcEdge.EPt.Z);
            Color = new double[srcEdge.Color.Length];
            srcEdge.Color.CopyTo(Color, 0);
        }

        public void SetVertexCoords(OpenTK.Vector3d sPt, OpenTK.Vector3d ePt)
        {
            SPt = sPt;
            EPt = ePt;
        }

        public OpenTK.Vector3d GetVertexCoord(bool isRoot)
        {
            return isRoot ? SPt : EPt;
        }

        public void SetVertexIds(uint vSId, uint vEId)
        {
            SVId = vSId;
            EVId = vEId;
        }

        public uint GetVertexId(bool isRoot)
        {
            return isRoot ? SVId : EVId;
        }

        public double GetCurveLength()
        {
            OpenTK.Vector3d h0 = SPt - EPt;
            return h0.Length;
        }

        public bool GetCurveAsPolyline(out IList<OpenTK.Vector3d> pts, int div)
        {
            pts = new List<OpenTK.Vector3d>();
            if (div <= 0)
            {
                // do nothing
                return true;
            }
            else
            {
                OpenTK.Vector3d tDiv = (EPt - SPt) * (1.0 / div);
                for (uint iDiv = 1; iDiv < div; iDiv++)
                {
                    OpenTK.Vector3d vec0 = SPt + tDiv * iDiv;
                    pts.Add(vec0);
                }
                System.Diagnostics.Debug.Assert(pts.Count <= div);
            }
            return true;
        }

        public OpenTK.Vector3d GetNearestPoint(OpenTK.Vector3d pt)
        {
            double t = CadUtils3D.FindNearestPointParameterLinePoint(pt, SPt, EPt);
            if (t < 0)
            {
                return SPt;
            }
            if (t > 1)
            {
                return EPt;
            }
            else
            {
                return SPt + (EPt - SPt) * t;
            }

            /*
            System.Diagnostics.Debug.Assert(false);
            OpenTK.Vector3d zeroV = new OpenTK.Vector3d();
            return zeroV;
            */
            throw new InvalidOperationException();
        }

        public bool Split(Edge3D addEdge, OpenTK.Vector3d addPt)
        {
            // TODO:

            return true;
        }
    }
}
