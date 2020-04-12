using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DFEMUtils
    {
        // Beam
        public static IvyFEM.Lapack.DoubleMatrix CalcBeamLocalKe(double le, double E, double I)
        {
            int n = 4;
            var ke = new IvyFEM.Lapack.DoubleMatrix(n, n);
            ke[0, 0] = 12.0;
            ke[0, 1] = 6.0 * le;
            ke[0, 2] = -12.0;
            ke[0, 3] = 6.0 * le;
            ke[1, 1] = 4.0 * le * le;
            ke[1, 2] = -6.0 * le;
            ke[1, 3] = 2.0 * le * le;
            ke[2, 2] = 12.0;
            ke[2, 3] = -6.0 * le;
            ke[3, 3] = 4.0 * le * le;
            for (int i = 0; i < n; i++)
            {
                for (int j = (i + 1); j < n; j++)
                {
                    ke[j, i] = ke[i, j];
                }
            }
            ke = IvyFEM.Lapack.DoubleMatrix.Scal(ke, E * I / (le * le * le));
            return ke;
        }

        public static IvyFEM.Lapack.DoubleMatrix CalcBeamLocalMe(double le, double rho, double Ae)
        {
            int n = 4;
            var me = new IvyFEM.Lapack.DoubleMatrix(n, n);
            me[0, 0] = 156.0;
            me[0, 1] = 22.0 * le;
            me[0, 2] = 54.0;
            me[0, 3] = -13.0 * le;
            me[1, 1] = 4.0 * le * le;
            me[1, 2] = 13.0 * le;
            me[1, 3] = -3.0 * le * le;
            me[2, 2] = 156.0;
            me[2, 3] = -22.0 * le;
            me[3, 3] = 4.0 * le * le;
            for (int i = 0; i < n; i++)
            {
                for (int j = (i + 1); j < n; j++)
                {
                    me[j, i] = me[i, j];
                }
            }
            me = IvyFEM.Lapack.DoubleMatrix.Scal(me, rho * Ae * le / 420.0);
            return me;
        }
    }
}
