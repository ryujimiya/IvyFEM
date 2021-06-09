using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Loop3D : IObject
    {
        public OpenTK.Vector3d Origin { get; set; }
        public OpenTK.Vector3d Normal { get; set; }
        public OpenTK.Vector3d XDir { get; set; }
        public IList<KeyValuePair<Edge2D, bool>> Edges { get; set; } = new List<KeyValuePair<Edge2D, bool>>();
        public IList<uint> EdgeIndexs { get; set; } = new List<uint>();
        public double[] Color { get; set; } = new double[3] { 0.8, 0.8, 0.8 };
        private BoundingBox3D BB = new BoundingBox3D();

        public Loop3D()
        {
            Origin = new OpenTK.Vector3d(0.0, 0.0, 0.0);
            Normal = new OpenTK.Vector3d(0.0, 0.0, 1.0);
            XDir = new OpenTK.Vector3d(1.0, 0.0, 0.0);
        }

        public Loop3D(OpenTK.Vector3d origin, OpenTK.Vector3d normal, OpenTK.Vector3d xdir)
        {
            Origin = new OpenTK.Vector3d(origin.X, origin.Y, origin.Z);
            Normal = OpenTK.Vector3d.Normalize(normal);
            XDir = OpenTK.Vector3d.Normalize(xdir);
            //--------------------------------------
            //DEBUG
            System.Diagnostics.Debug.Assert(!double.IsNaN(Normal.X));
            System.Diagnostics.Debug.Assert(!double.IsNaN(Normal.Y));
            //--------------------------------------
        }

        public Loop3D(Loop3D src)
        {
            Copy(src);
        }

        public void Copy(IObject src)
        {
            Loop3D srcLoop = src as Loop3D;

            Origin = new OpenTK.Vector3d(srcLoop.Origin.X, srcLoop.Origin.Y, srcLoop.Origin.Z);
            Normal = new OpenTK.Vector3d(srcLoop.Normal.X, srcLoop.Normal.Y, srcLoop.Normal.Z);
            XDir = new OpenTK.Vector3d(srcLoop.XDir.X, srcLoop.XDir.Y, srcLoop.XDir.Z);
            Edges = new List<KeyValuePair<Edge2D, bool>>();
            foreach (var edge in srcLoop.Edges)
            {
                Edge2D e = edge.Key;
                bool value = edge.Value;
                Edges.Add(new KeyValuePair<Edge2D, bool>(e, value));
            }
            EdgeIndexs = new List<uint>();
            foreach (uint index in  srcLoop.EdgeIndexs)
            {
                EdgeIndexs.Add(index);
            }
            Color = new double[srcLoop.Color.Length];
            srcLoop.Color.CopyTo(Color, 0);
            BB = new BoundingBox3D(srcLoop.BB);
        }

        public OpenTK.Vector2d Project(OpenTK.Vector3d p)
        {
            double x = OpenTK.Vector3d.Dot((p - Origin), XDir);
            double y = OpenTK.Vector3d.Dot((p - Origin), OpenTK.Vector3d.Cross(Normal, XDir));
            return new OpenTK.Vector2d(x, y);
        }
        
        public OpenTK.Vector3d UnProject(OpenTK.Vector2d p)
        {
            return Origin + XDir * p.X + OpenTK.Vector3d.Cross(Normal, XDir) * p.Y;
        }

        public BoundingBox3D GetBoundingBox()
        {
            if (BB.IsntEmpty)
            {
                return BB;
            }
            BoundingBox2D bb2 = new BoundingBox2D();
            for(int i = 0; i < Edges.Count; i++)
            {
                var edge = Edges[i];
                Edge2D e = edge.Key;
                bb2 += e.GetBoundingBox();
            }
            double ep0 = (bb2.MaxX - bb2.MinX + bb2.MaxY - bb2.MinY) * 1.0e-10;
            BB.AddPoint(UnProject(new OpenTK.Vector2d(bb2.MinX, bb2.MinY) ), ep0 );
            BB.AddPoint(UnProject(new OpenTK.Vector2d(bb2.MaxX,bb2.MinY) ), ep0 );
            BB.AddPoint(UnProject(new OpenTK.Vector2d(bb2.MinX,bb2.MaxY) ), ep0 );
            BB.AddPoint(UnProject(new OpenTK.Vector2d(bb2.MaxX,bb2.MaxY) ), ep0 );
            return BB;
        }

        public void ClearBoundingBox()
        {
            BB = new BoundingBox3D();
        }
    }
}
