using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TriDiaMat3
    {
        private uint n;
        private double[] v;

        public TriDiaMat3(uint n)
        {
            this.n = n;
            v = new double[(n * 3 - 2) * 9];
        }

        public void Clear()
        {
            for (uint i = 0; i < (n * 3 - 2) * 9; i++)
            {
                v[i] = 0;
            }
        }

        // eM[][2][3][3]
        public void Merge(uint idiv, double[][][][] eM)
        {
            //for (uint i = 0; i < 36; i++)
            //{
            //    v[idiv * 27 + i] += (&eM[0][0][0][0])[i];
            //}

            uint index = 0;
            for (uint p = 0; p < eM.Length; p++)
            {
                for (uint q = 0; q < 2; q++)
                {
                    for (uint r = 0; r < 3; r++)
                    {
                        for (uint s = 0; s < 3; s++)
                        {
                            v[idiv * 27 + index] += eM[p][q][r][s];
                            index++;
                        }
                    }
                }
            }
            System.Diagnostics.Debug.Assert(index == 36);
        }

        public void FixBoundaryCondition(uint ino, uint idof)
        {
            System.Diagnostics.Debug.Assert(idof < 3 && ino < n);
            if (ino != 0)
            {
                //double* pvu = v + ino * 27 - 18;
                //double* pvl = v + ino * 27 - 9;

                int ivu = (int)(ino * 27 - 18);
                v[ivu + 0 * 3 + idof] = 0;
                v[ivu + 1 * 3 + idof] = 0;
                v[ivu + 2 * 3 + idof] = 0;

                int ivl = (int)(ino * 27 - 9);
                v[ivl + idof * 3 + 0] = 0;
                v[ivl + idof * 3 + 1] = 0;
                v[ivl + idof * 3 + 2] = 0;
            }
            if (ino != n - 1)
            {
                //double* pvu = v + ino * 27 + 18;
                //double* pvl = v + ino * 27 + 9;

                int ivu = (int)(ino * 27 + 18);
                v[ivu + 0 * 3 + idof] = 0;
                v[ivu + 1 * 3 + idof] = 0;
                v[ivu + 2 * 3 + idof] = 0;

                int ivl = (int)(ino * 27 + 9);
                v[ivl + idof * 3 + 0] = 0;
                v[ivl + idof * 3 + 1] = 0;
                v[ivl + idof * 3 + 2] = 0;
            }

            //double* pvc = v + ino * 27;
            int ivc = (int)(ino * 27);
            v[ivc + 0 * 3 + idof] = 0;
            v[ivc + 1 * 3 + idof] = 0;
            v[ivc + 2 * 3 + idof] = 0;
            v[ivc + idof * 3 + 0] = 0;
            v[ivc + idof * 3 + 1] = 0;
            v[ivc + idof * 3 + 2] = 0;
            v[ivc + idof * 3 + idof] = 1;
        }

        // execute ILU factorization
        public void ILUFrac()
        {
            double[] tmpBlk = new double[9];
            for (uint iblk = 0; iblk < n; iblk++)
            {
                if (iblk != 0)
                {
                    //double* pVal_ik = v + 27 * iblk - 9;
                    //double* pVal_kj = v + 27 * iblk - 18;
                    //double* pVal_ij = v + 27 * iblk;
                    double[] valik = new double[9];
                    for (int p = 0; p < 9; p++)
                    {
                        valik[p] = v[p + 27 * iblk - 9];
                    }
                    double[] valkj = new double[9];
                    for (int p = 0; p < 9; p++)
                    {
                        valkj[p] = v[p + 27 * iblk - 18];
                    }
                    double[] valij = new double[9];
                    for (int p = 0; p < 9; p++)
                    {
                        valij[p] = v[p + 27 * iblk];
                    }
                    for (uint i = 0; i < 3; i++)
                    {
                        valij[i * 3 + 0] -= valik[i * 3 + 0] * valkj[0] + valik[i * 3 + 1] * valkj[3] + valik[i * 3 + 2] * valkj[6];
                        valij[i * 3 + 1] -= valik[i * 3 + 0] * valkj[1] + valik[i * 3 + 1] * valkj[4] + valik[i * 3 + 2] * valkj[7];
                        valij[i * 3 + 2] -= valik[i * 3 + 0] * valkj[2] + valik[i * 3 + 1] * valkj[5] + valik[i * 3 + 2] * valkj[8];
                    }
                    for (int p = 0; p < 9; p++)
                    {
                        v[p + 27 * iblk] = valij[p];
                    }
                }
                {
                    // calc inverse of diagonal
                    //double* pVal_ii = v + 27 * iblk;
                    double[] valii = new double[9];
                    for (int p = 0; p < 9; p++)
                    {
                        valii[p] = v[p + 27 * iblk];
                    }
                    CalcInvMat3(valii, tmpBlk);
                    for (int p = 0; p < 9; p++)
                    {
                        v[p + 27 * iblk] = valii[p];
                    }
                }
                if (iblk != n - 1)
                {
                    //double* pVal_ij = v + 27 * iblk + 9;
                    //double* pVal_ii = v + 27 * iblk;
                    double[] valij = new double[9];
                    for (int p = 0; p < 9; p++)
                    {
                        valij[p] = v[p + 27 * iblk + 9];
                    }
                    double[] valii = new double[9];
                    for (int p = 0; p < 9; p++)
                    {
                        valij[p] = v[p + 27 * iblk];
                    }
                    for (uint i = 0; i < 9; i++)
                    {
                        tmpBlk[i] = valij[i];
                    }
                    for (uint i = 0; i < 3; i++)
                    {
                        valij[i * 3 + 0] = valii[i * 3 + 0] * tmpBlk[0] + valii[i * 3 + 1] * tmpBlk[3] + valii[i * 3 + 2] * tmpBlk[6];
                        valij[i * 3 + 1] = valii[i * 3 + 0] * tmpBlk[1] + valii[i * 3 + 1] * tmpBlk[4] + valii[i * 3 + 2] * tmpBlk[7];
                        valij[i * 3 + 2] = valii[i * 3 + 0] * tmpBlk[2] + valii[i * 3 + 1] * tmpBlk[5] + valii[i * 3 + 2] * tmpBlk[8];
                    }
                    for (int p = 0; p < 9; p++)
                    {
                        v[p + 27 * iblk + 9] = valij[p];
                    }
                }
            }   // end iblk
        }

        // solve matrix
        public void Solve(double[] res)
        {
            double[] tmpVec = new double[3];
            for (uint iblk = 0; iblk < n; iblk++)
            {
                tmpVec[0] = res[iblk * 3 + 0];
                tmpVec[1] = res[iblk * 3 + 1];
                tmpVec[2] = res[iblk * 3 + 2];
                if (iblk != 0)
                {
                    //double* pVal_ij = v + iblk * 27 - 9;
                    double[] valij = new double[9];
                    for (int p = 0; p < 9; p++)
                    {
                        valij[p] = v[p + iblk * 27 - 9];
                    }
                    double valj0 = res[(iblk - 1) * 3 + 0];
                    double valj1 = res[(iblk - 1) * 3 + 1];
                    double valj2 = res[(iblk - 1) * 3 + 2];
                    tmpVec[0] -= valij[0] * valj0 + valij[1] * valj1 + valij[2] * valj2;
                    tmpVec[1] -= valij[3] * valj0 + valij[4] * valj1 + valij[5] * valj2;
                    tmpVec[2] -= valij[6] * valj0 + valij[7] * valj1 + valij[8] * valj2;
                }
                //double* pVal_ii = v + 27 * iblk;
                double[] valii = new double[9];
                for (int p = 0; p < 9; p++)
                {
                    valii[p] = v[p + 27 * iblk];
                }
                res[iblk * 3 + 0] = valii[0] * tmpVec[0] + valii[1] * tmpVec[1] + valii[2] * tmpVec[2];
                res[iblk * 3 + 1] = valii[3] * tmpVec[0] + valii[4] * tmpVec[1] + valii[5] * tmpVec[2];
                res[iblk * 3 + 2] = valii[6] * tmpVec[0] + valii[7] * tmpVec[1] + valii[8] * tmpVec[2];
            }
            for (int iblk = (int)n - 1; iblk >= 0; iblk--)
            {
                tmpVec[0] = res[iblk * 3 + 0];
                tmpVec[1] = res[iblk * 3 + 1];
                tmpVec[2] = res[iblk * 3 + 2];
                if (iblk != (int)n - 1)
                {
                    //double* pVal_ij = v + 27 * iblk + 9;
                    double[] valij = new double[9];
                    for (int p = 0; p < 9; p++)
                    {
                        valij[p] = v[p + 27 * iblk + 9];
                    }
                    double valj0 = res[(iblk + 1) * 3 + 0];
                    double valj1 = res[(iblk + 1) * 3 + 1];
                    double valj2 = res[(iblk + 1) * 3 + 2];
                    tmpVec[0] -= valij[0] * valj0 + valij[1] * valj1 + valij[2] * valj2;
                    tmpVec[1] -= valij[3] * valj0 + valij[4] * valj1 + valij[5] * valj2;
                    tmpVec[2] -= valij[6] * valj0 + valij[7] * valj1 + valij[8] * valj2;
                }
                res[iblk * 3 + 0] = tmpVec[0];
                res[iblk * 3 + 1] = tmpVec[1];
                res[iblk * 3 + 2] = tmpVec[2];
            }
        }

        private static void CalcInvMat3(double[] a, double[] t)
        {
            double det = a[0] * a[4] * a[8] + a[3] * a[7] * a[2] + a[6] * a[1] * a[5]
                - a[0] * a[7] * a[5] - a[6] * a[4] * a[2] - a[3] * a[1] * a[8];
            double inv_det = 1.0 / det;
            for (int i = 0; i < 9; i++)
            {
                t[i] = a[i];
            }
            a[0] = inv_det * (t[4] * t[8] - t[5] * t[7]);
            a[1] = inv_det * (t[2] * t[7] - t[1] * t[8]);
            a[2] = inv_det * (t[1] * t[5] - t[2] * t[4]);
            a[3] = inv_det * (t[5] * t[6] - t[3] * t[8]);
            a[4] = inv_det * (t[0] * t[8] - t[2] * t[6]);
            a[5] = inv_det * (t[2] * t[3] - t[0] * t[5]);
            a[6] = inv_det * (t[3] * t[7] - t[4] * t[6]);
            a[7] = inv_det * (t[1] * t[6] - t[0] * t[7]);
            a[8] = inv_det * (t[0] * t[4] - t[1] * t[3]);
        }
    }

}
