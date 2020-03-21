using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DFEMUtils
    {
        // Frame
        public static IvyFEM.Lapack.DoubleMatrix CalcFrameLocalKe(double le, double E, double Ae, double I)
        {
            int n = 6;
            var ke = new IvyFEM.Lapack.DoubleMatrix(n, n);
            double c1 = E * Ae / le;
            double c2 = E * I / (le * le * le);

            ke[0, 0] = c1;
            ke[0, 3] = c1 * (-1.0);
            ke[1, 1] = c2 * 12.0;
            ke[1, 2] = c2 * 6.0 * le;
            ke[1, 4] = c2 * (-12.0);
            ke[1, 5] = c2 * 6.0 * le;
            ke[2, 2] = c2 * 4.0 * le * le;
            ke[2, 4] = c2 * (-6.0 * le);
            ke[2, 5] = c2 * 2.0 * le * le;
            ke[3, 3] = c1;
            ke[4, 4] = c2 * 12.0;
            ke[4, 5] = c2 * (-6.0 * le);
            ke[5, 5] = c2 * 4.0 * le * le;
            for (int i = 0; i < n; i++)
            {
                for (int j = (i + 1); j < n; j++)
                {
                    ke[j, i] = ke[i, j];
                }
            }
            return ke;
        }

        public static IvyFEM.Lapack.DoubleMatrix CalcFrameLocalMe(double le, double rho, double Ae)
        {
            int n = 6;
            var me = new IvyFEM.Lapack.DoubleMatrix(n, n);
            double c1 = rho * Ae * le / 6.0;
            double c2 = rho * Ae * le / 420.0;
            me[0, 0] = c1 * 2.0;
            me[0, 3] = c1;
            me[1, 1] = c2 * 156.0;
            me[1, 2] = c2 * 22.0 * le;
            me[1, 4] = c2 * 54.0;
            me[1, 5] = c2 * (-13.0 * le);
            me[2, 2] = c2 * 4.0 * le * le;
            me[2, 4] = c2 * 13.0 * le;
            me[2, 5] = c2 * (-3.0 * le * le);
            me[3, 3] = c1 * 2.0;
            me[4, 4] = c2 * 156.0;
            me[4, 5] = c2 * (-22.0 * le);
            me[5, 5] = c2 * 4.0 * le * le;
            for (int i = 0; i < n; i++)
            {
                for (int j = (i + 1); j < n; j++)
                {
                    me[j, i] = me[i, j];
                }
            }
            return me;
        }
    }
}
