using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class Elastic2DBaseFEM
    {
        protected void SetMultipointConstraintSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            for (uint quantityId = 0; quantityId < World.GetQuantityCount(); quantityId++)
            {
                int cnt = World.GetMultipointConstraintCount(quantityId);
                if (cnt > 0)
                {
                    SetMultipointConstraintQuantitySpecialBC(quantityId, A, B);
                }
            }
        }

        private void SetMultipointConstraintQuantitySpecialBC(
            uint cQuantityId, IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            uint uQuantityId = 0;
            System.Diagnostics.Debug.Assert(World.GetCoordCount(uQuantityId) ==
                World.GetCoordCount(cQuantityId));
            int coCnt = (int)World.GetCoordCount(cQuantityId);
            System.Diagnostics.Debug.Assert(World.GetDof(uQuantityId) == 2);
            System.Diagnostics.Debug.Assert(World.GetDof(cQuantityId) == 1);
            int uDof = 2;
            int cDof = 1;
            int uNodeCnt = (int)World.GetNodeCount(uQuantityId);
            int cNodeCnt = (int)World.GetNodeCount(cQuantityId);
            int offset = GetOffset(cQuantityId);

            for (int coId = 0; coId < coCnt; coId++)
            {
                int cNodeId = World.Coord2Node(cQuantityId, coId);
                if (cNodeId == -1)
                {
                    continue;
                }
                int uNodeId = World.Coord2Node(uQuantityId, coId);
                System.Diagnostics.Debug.Assert(uNodeId != -1);

                double[] coord = World.GetCoord(cQuantityId, coId);
                double[] curCoord = new double[uDof];
                for (int i = 0; i < uDof; i++)
                {
                    curCoord[i] = coord[i] + U[uNodeId * uDof + i];
                }
                // Lagrangeの未定乗数
                double c = U[offset + cNodeId];

                // 多点拘束
                var mpcs = World.GetMultipointConstraintFromCoord(cQuantityId, coId);
                foreach (MultipointConstraint mpc in mpcs)
                {
                    Constraint constraint = mpc.Constraint;
                    double G = constraint.GetValue(curCoord);
                    double[] dG = new double[uDof];
                    for (int iDof = 0; iDof < uDof; iDof++)
                    {
                        dG[iDof] = constraint.GetDerivative(iDof, curCoord);
                    }
                    double[,] d2G = new double[uDof, uDof];
                    for (int iDof = 0; iDof < uDof; iDof++)
                    {
                        for (int jDof = 0; jDof < uDof; jDof++)
                        {
                            d2G[iDof, jDof] = constraint.Get2ndDerivative(iDof, jDof, curCoord);
                        }
                    }

                    if (constraint.Equality == EqualityType.Eq)
                    {
                        CalcEqConstraitsAB(c, G, dG, d2G,
                            uDof, uNodeId, offset, cNodeId,
                            A, B);
                    }
                    else
                    {
                        CalcLessEqualConstraitsAB(c, G, dG, d2G,
                            uDof, uNodeId, offset, cNodeId,
                            A, B);
                    }
                }
            }
        }

        private void CalcEqConstraitsAB(double c, double G, double[] dG, double[,] d2G,
            int uDof, int uNodeId, int offset, int cNodeId,
            IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            {
                double[,] kuu = new double[uDof, uDof];
                for (int iDof = 0; iDof < uDof; iDof++)
                {
                    for (int jDof = 0; jDof < uDof; jDof++)
                    {
                        kuu[iDof, jDof] = c * d2G[iDof, jDof];
                    }
                }
                for (int rowDof = 0; rowDof < uDof; rowDof++)
                {
                    for (int colDof = 0; colDof < uDof; colDof++)
                    {
                        A[uNodeId * uDof + rowDof, uNodeId * uDof + colDof] +=
                            kuu[rowDof, colDof];
                        B[uNodeId * uDof + rowDof] +=
                            kuu[rowDof, colDof] * U[uNodeId * uDof + colDof];
                    }
                }
            }

            {
                double[,] kuc = new double[uDof, 1];
                double[,] kcu = new double[1, uDof];
                for (int iDof = 0; iDof < uDof; iDof++)
                {
                    kuc[iDof, 0] = dG[iDof];
                    kcu[0, iDof] = dG[iDof];
                }
                for (int rowDof = 0; rowDof < uDof; rowDof++)
                {
                    A[uNodeId * uDof + rowDof, offset + cNodeId] += kuc[rowDof, 0];
                    A[offset + cNodeId, uNodeId * uDof + rowDof] += kcu[0, rowDof];
                    B[uNodeId * uDof + rowDof] +=
                        kuc[rowDof, 0] * U[offset + cNodeId];
                    B[offset + cNodeId] +=
                        kcu[0, rowDof] * U[uNodeId * uDof + rowDof];
                }
            }

            {
                double[] qu = new double[uDof];
                for (int iDof = 0; iDof < uDof; iDof++)
                {
                    qu[iDof] = c * dG[iDof];
                }
                double qc = G;

                for (int iDof = 0; iDof < uDof; iDof++)
                {
                    B[uNodeId * uDof + iDof] += -qu[iDof];
                }
                B[offset + cNodeId] += -qc;
            }
        }

        private void CalcLessEqualConstraitsAB(double c, double G, double[] dG, double[,] d2G,
            int uDof, int uNodeId, int offset, int cNodeId,
            IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            double tolerance = IvyFEM.Linear.Constants.ConvRatioTolerance;
            if (c <= tolerance &&
                G <= tolerance)
            {
                // unconstrained
                {
                    // kcc
                    double kcc = 1.0;
                    A[offset + cNodeId, offset + cNodeId] += kcc;
                    B[offset + cNodeId] += kcc * U[offset + cNodeId];
                }
                {
                    // qc
                    double qc = c;
                    B[offset + cNodeId] += -qc;
                }
            }
            else
            {
                // constrained
                CalcEqConstraitsAB(c, G, dG, d2G,
                    uDof, uNodeId, offset, cNodeId, A, B);
            }
        }
    }
}
