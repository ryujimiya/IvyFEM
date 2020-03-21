using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DFEMUtils
    {
        // Truss
        public static IvyFEM.Lapack.DoubleMatrix CalcTrussLocalKe(double le, double E, double Ae)
        {
            int n = 2;
            var ke = new IvyFEM.Lapack.DoubleMatrix(n, n);
            ke[0, 0] = 1.0;
            ke[0, 1] = -1.0;
            ke[1, 1] = 1.0;
            for (int i = 0; i < n; i++)
            {
                for (int j = (i + 1); j < n; j++)
                {
                    ke[j, i] = ke[i, j];
                }
            }
            IvyFEM.Lapack.Functions.dscal(ke.Buffer, E * Ae / le);
            return ke;
        }

        public static IvyFEM.Lapack.DoubleMatrix CalcTrussLocalMe(double le, double rho, double Ae)
        {
            int n = 2;
            var me = new IvyFEM.Lapack.DoubleMatrix(n, n);
            me[0, 0] = 2.0;
            me[0, 1] = 1.0;
            me[1, 1] = 2.0;
            for (int i = 0; i < n; i++)
            {
                for (int j = (i + 1); j < n; j++)
                {
                    me[j, i] = me[i, j];
                }
            }
            IvyFEM.Lapack.Functions.dscal(me.Buffer, rho * Ae * le / 6.0);
            return me;
        }
    }
}
