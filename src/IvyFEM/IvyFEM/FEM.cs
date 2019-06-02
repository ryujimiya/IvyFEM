using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class FEM
    {
        public FEWorld World { get; set; } = null;
        public IvyFEM.Linear.IEquationSolver Solver { get; set; } = null;

        public abstract void Solve();
    }
}
