using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Solid : IObject
    {
        public IList<uint> LoopIds { get; set; } = new List<uint>();
        public IList<OpenTK.Vector3d> Holes { get; set; } = new List<OpenTK.Vector3d>();
        public IList<uint> InsideVIds { get; set; } = new List<uint>();

        public Solid()
        {

        }

        public Solid(Solid src)
        {
            Copy(src);
        }

        public Solid(IList<uint> lIds, IList<OpenTK.Vector3d> holes, IList<uint> insideVIds)
        {
            LoopIds = new List<uint>(lIds);
            Holes = new List<OpenTK.Vector3d>(holes);
            InsideVIds = new List<uint>(insideVIds);
        }

        public void Copy(IObject src)
        {
            Solid srcSolid = src as Solid;

            LoopIds = new List<uint>(srcSolid.LoopIds);
            Holes = new List<OpenTK.Vector3d>(srcSolid.Holes);
            InsideVIds = new List<uint>(srcSolid.InsideVIds);
        }
    }
}
