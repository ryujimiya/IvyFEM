using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class FluidDisplacement2DTDFEM
    {
        private void SetMultipointConstraintSpecialBC(IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
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
            double dt = TimeStep;
            uint vQuantityId = 0;
            uint pQuantityId = 1;
            int vDof = 2;
            int pDof = 1;
            int cDof = 1;
            int vNodeCnt = (int)World.GetNodeCount(vQuantityId);
            int pNodeCnt = (int)World.GetNodeCount(pQuantityId);
            int pOffset = vNodeCnt * vDof;
            int cOffset = pOffset + pNodeCnt;
            System.Diagnostics.Debug.Assert(World.GetCoordCount(vQuantityId) ==
                World.GetCoordCount(cQuantityId));
            int coCnt = (int)World.GetCoordCount(cQuantityId);
            System.Diagnostics.Debug.Assert(World.GetDof(vQuantityId) == 2);
            System.Diagnostics.Debug.Assert(World.GetDof(cQuantityId) == 1);

            for (int coId = 0; coId < coCnt; coId++)
            {
                int cNodeId = World.Coord2Node(cQuantityId, coId);
                if (cNodeId == -1)
                {
                    continue;
                }
                int vNodeId = World.Coord2Node(vQuantityId, coId);
                System.Diagnostics.Debug.Assert(vNodeId != -1);

                //double[] coord = World.GetCoord(cQuantityId, coId);
                double[] curCoord = new double[vDof];
                for (int i = 0; i < vDof; i++)
                {
                    curCoord[i] = UpdatedCoord[vNodeId * vDof + i] + U[vNodeId * vDof + i] * dt;
                }
                // Lagrangeの未定乗数
                double c = U[cOffset + cNodeId];

                // 多点拘束
                var mpcs = World.GetMultipointConstraintFromCoord(cQuantityId, coId);
                foreach (MultipointConstraint mpc in mpcs)
                {
                    Constraint constraint = mpc.Constraint;
                    double G = constraint.GetValue(curCoord);
                    double[] dG = new double[vDof];
                    for (int iDof = 0; iDof < vDof; iDof++)
                    {
                        dG[iDof] = constraint.GetDerivative(iDof, curCoord);
                    }
                    double[,] d2G = new double[vDof, vDof];
                    for (int iDof = 0; iDof < vDof; iDof++)
                    {
                        for (int jDof = 0; jDof < vDof; jDof++)
                        {
                            d2G[iDof, jDof] = constraint.Get2ndDerivative(iDof, jDof, curCoord);
                        }
                    }

                    System.Diagnostics.Debug.Assert(constraint.Equality != EqualityType.Eq);
                    CalcLessEqualConstraitsAB(c, G, dG, d2G,
                        vDof, vNodeId, cOffset, cNodeId,
                        A, B);
                }
            }
        }

        // Newton-Raphson
        private void CalcEqConstraitsAB(double c, double G, double[] dG, double[,] d2G,
            int vDof, int vNodeId, int cOffset, int cNodeId,
            IvyFEM.Linear.DoubleSparseMatrix A, double[] B)
        {
            double dt = TimeStep;
            {
                double[,] kuu = new double[vDof, vDof];
                for (int iDof = 0; iDof < vDof; iDof++)
                {
                    for (int jDof = 0; jDof < vDof; jDof++)
                    {
                        kuu[iDof, jDof] = c * d2G[iDof, jDof] * dt * dt;
                    }
                }
                for (int rowDof = 0; rowDof < vDof; rowDof++)
                {
                    for (int colDof = 0; colDof < vDof; colDof++)
                    {
                        A[vNodeId * vDof + rowDof, vNodeId * vDof + colDof] +=
                            kuu[rowDof, colDof];
                        B[vNodeId * vDof + rowDof] +=
                            kuu[rowDof, colDof] * U[vNodeId * vDof + colDof];
                    }
                }
            }

            {
                double[,] kuc = new double[vDof, 1];
                double[,] kcu = new double[1, vDof];
                for (int iDof = 0; iDof < vDof; iDof++)
                {
                    kuc[iDof, 0] = dG[iDof] * dt;
                    kcu[0, iDof] = dG[iDof] * dt;
                }
                for (int rowDof = 0; rowDof < vDof; rowDof++)
                {
                    A[vNodeId * vDof + rowDof, cOffset + cNodeId] += kuc[rowDof, 0];
                    A[cOffset + cNodeId, vNodeId * vDof + rowDof] += kcu[0, rowDof];
                    B[vNodeId * vDof + rowDof] +=
                        kuc[rowDof, 0] * U[cOffset + cNodeId];
                    B[cOffset + cNodeId] +=
                        kcu[0, rowDof] * U[vNodeId * vDof + rowDof];
                }
            }

            {
                double[] qu = new double[vDof];
                for (int iDof = 0; iDof < vDof; iDof++)
                {
                    qu[iDof] = c * dG[iDof] * dt;
                }
                double qc = G;

                for (int iDof = 0; iDof < vDof; iDof++)
                {
                    B[vNodeId * vDof + iDof] += -qu[iDof];
                }
                B[cOffset + cNodeId] += -qc;
            }
        }

        private void CalcLessEqualConstraitsAB(double c, double G, double[] dG, double[,] d2G,
            int vDof, int vNodeId, int cOffset, int cNodeId,
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
                    A[cOffset + cNodeId, cOffset + cNodeId] += kcc;
                    B[cOffset + cNodeId] += kcc * U[cOffset + cNodeId];
                }
                {
                    // qc
                    double qc = c;
                    B[cOffset + cNodeId] += -qc;
                }
            }
            else
            {
                // constrained
                CalcEqConstraitsAB(c, G, dG, d2G,
                    vDof, vNodeId, cOffset, cNodeId, A, B);
            }
        }
    }
}
