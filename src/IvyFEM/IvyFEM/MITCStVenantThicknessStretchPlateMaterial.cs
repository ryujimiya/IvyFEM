﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MITCStVenantThicknessStretchPlateMaterial : PlateBaseMaterial
    {
        public MITCStVenantThicknessStretchPlateMaterial() : base()
        {

        }

        public MITCStVenantThicknessStretchPlateMaterial(MITCStVenantThicknessStretchPlateMaterial src) : base(src)
        {

        }

        public override void Copy(IObject src)
        {
            base.Copy(src);
        }
    }
}