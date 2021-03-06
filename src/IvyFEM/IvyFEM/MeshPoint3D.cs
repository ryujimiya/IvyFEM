﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MeshPoint3D
    {
        public OpenTK.Vector3d Point { get; set; } = new OpenTK.Vector3d();
        public int Element { get; set; } = -1;
        public int Node { get; set; } = -1;
        public int OldNode { get; set; } = -1;

        public MeshPoint3D()
        {

        }

        public MeshPoint3D(MeshPoint3D src)
        {
            Point = src.Point;
            Element = src.Element;
            Node = src.Node;
            OldNode = src.OldNode;
        }
    }

}
