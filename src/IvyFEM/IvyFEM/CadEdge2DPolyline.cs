using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class CadEdge2DPolyline
    {
        private OpenTK.Vector2d PickPos;
        private int PickedDivIndex;
        private uint ECadId;
        ////
        private uint No;
        private double[] Ut;
        private double[] IniX;
        private int[] BCFlag;
        ////////////////
        private double EI;
        private double ARho;
        private double AE;

        private double[] Dut;
        private double[] Res;
        private TriDiaMat3 Mat;

        public CadEdge2DPolyline()
        {
            EI = 30.0;
            ARho = 1;
            AE = 100000.0;
            ////
            No = 0;
            ////
            Ut = null;
            IniX = null;
            BCFlag = null;
            Dut = null;
            Res = null;
            Mat = null;
        }

        private void ClearMemory()
        {
            IniX = null;
            Ut = null;
            Dut = null;
            Res = null;
            Dut = null;
            BCFlag = null;
            Mat = null;
        }

        private void SetFixedBoundaryFlag(uint ino, uint idim)
        {
            System.Diagnostics.Debug.Assert(ino < No);
            System.Diagnostics.Debug.Assert(idim < 3);
            System.Diagnostics.Debug.Assert(BCFlag != null);
            BCFlag[ino * 3 + idim] = 1;
        }

        private void SetDisp(uint ino, uint idim, double disp)
        {
            if (ino >= No)
            {
                return;
            }
            System.Diagnostics.Debug.Assert(idim < 3);
            System.Diagnostics.Debug.Assert(Ut != null);
            Ut[ino * 3 + idim] = disp;
        }

        private void GetValueNode(uint ino, out double x, out double y, out double t)
        {
            System.Diagnostics.Debug.Assert(ino < No);
            x = Ut[ino * 3 + 0] + IniX[ino * 2 + 0];
            y = Ut[ino * 3 + 1] + IniX[ino * 2 + 1];
            t = Ut[ino * 3 + 2];
        }

        public void SetCadEdge(CadObject2D cad2D, uint eId, OpenTK.Vector2d pickPos)
        {
            this.PickPos = pickPos;
            this.ECadId = eId;
            ClearMemory();
            IList<double> relCoPolys = new List<double>();
            {
                System.Diagnostics.Debug.Assert(cad2D.IsElemId(CadElementType.Edge, eId));
                Edge2D e = cad2D.GetEdge(eId);
                relCoPolys = e.GetCurveRelPoint();
            }
            IList<double> xys = new List<double>();
            uint ndiv = (uint)(relCoPolys.Count / 2 + 1);
            No = ndiv + 1;
            for (int i = 0; i < No * 2; i++)
            {
                xys.Add(0);
            }
            uint sVId = cad2D.GetEdgeVertexId(eId, true);
            uint eVId = cad2D.GetEdgeVertexId(eId, false);
            OpenTK.Vector2d sV = cad2D.GetVertex(sVId);
            OpenTK.Vector2d eV = cad2D.GetVertex(eVId);
            OpenTK.Vector2d hse = eV - sV;
            OpenTK.Vector2d vse = new OpenTK.Vector2d(-hse.Y, hse.X);
            xys.Add(sV.X);
            xys.Add(sV.Y);
            for (int i = 0; i < ndiv - 1; i++)
            {
                OpenTK.Vector2d v = sV + relCoPolys[i * 2 + 0] * hse + relCoPolys[i * 2 + 1] * vse;
                xys.Add(v.X);
                xys.Add(v.Y);
            }
            xys.Add(eV.X);
            xys.Add(eV.Y);
            //////////////////////////////////////
            Mat = new TriDiaMat3(No);
            IniX = new double[No * 2];
            Ut = new double[No * 3];
            Res = new double[No * 3];
            Dut = new double[No * 3];
            BCFlag = new int[No * 3];
            ////
            for (int i = 0; i < No * 2; i++)
            {
                IniX[i] = xys[i];
            }
            for (int i = 0; i < No * 3; i++)
            {
                Ut[i] = 0;
            }
            for (int i = 0; i < No * 3; i++)
            {
                BCFlag[i] = 0;
            }
            //////////////////////////////////////
            SetFixedBoundaryFlag(0, 0);
            SetFixedBoundaryFlag(0, 1);
            SetFixedBoundaryFlag(No - 1, 0);
            SetFixedBoundaryFlag(No - 1, 1);
            double alpha;
            double ndist;
            double normX;
            double normY;
            ProjectPoint(pickPos.X, pickPos.Y, ref PickedDivIndex,
                           out alpha, out ndist, out normX, out normY);
            if (PickedDivIndex == -1)
            {
                return;
            }
            SetFixedBoundaryFlag((uint)PickedDivIndex, 0);
            SetFixedBoundaryFlag((uint)PickedDivIndex, 1);
            SetFixedBoundaryFlag((uint)PickedDivIndex + 1, 0);
            SetFixedBoundaryFlag((uint)PickedDivIndex + 1, 1);
        }

        public void Drag(CadObject2D cad2D, OpenTK.Vector2d distPos)
        {
            OpenTK.Vector2d del = distPos - PickPos;
            SetDisp((uint)PickedDivIndex, 0, del.X);
            SetDisp((uint)PickedDivIndex, 1, del.Y);
            SetDisp((uint)PickedDivIndex + 1, 0, del.X);
            SetDisp((uint)PickedDivIndex + 1, 1, del.Y);

            SolveLinearStatic();
            IList<OpenTK.Vector2d> xys = new List<OpenTK.Vector2d>();
            for (uint i = 1; i < No - 1; i++)
            {
                double x;
                double y;
                double t;
                GetValueNode(i, out x, out y, out t);
                xys.Add(new OpenTK.Vector2d(x, y));
            }
            cad2D.SetCurvePolyline(ECadId, xys);
        }

        private void ProjectPoint(double inX, double inY,
            ref int minDivIndex,
            out double alpha, out double ndist, out double normX, out double normY)
        {
            // dummy
            alpha = 0;
            ndist = 0;
            normX = 0;
            normY = 0;

            double[] x0 = { inX, inY };
            double dist = Math.Sqrt((IniX[0] + Ut[0] - x0[0]) * (IniX[0] + Ut[0] - x0[0])
                               + (IniX[1] + Ut[1] - x0[1]) * (IniX[1] + Ut[1] - x0[1]));
            minDivIndex = -1;
            uint ndiv = No - 1;
            for (int idiv = 0; idiv < ndiv; idiv++)
            {
                double[] x1 = { IniX[idiv * 2 + 0] + Ut[idiv * 3 + 0], IniX[idiv * 2 + 1] + Ut[idiv * 3 + 1] };
                double[] x2 = { IniX[idiv * 2 + 2] + Ut[idiv * 3 + 3], IniX[idiv * 2 + 3] + Ut[idiv * 3 + 4] };
                double t = CadUtils.FindNearestPointParameterLinePoint(
                    new OpenTK.Vector2d(x0[0], x0[1]),
                    new OpenTK.Vector2d(x1[0], x1[1]),
                    new OpenTK.Vector2d(x2[0], x2[0]));
                if (t < -0.001 || t > 1.001)
                {
                    continue;
                }
                double[] x3 = { x1[0] * (1 - t) + x2[0] * t, x1[1] * (1 - t) + x2[1] * t };
                double[] d = { x0[0] - x3[0], x0[1] - x3[1] };
                double elen = Math.Sqrt((x1[0] - x2[0]) * (x1[0] - x2[0]) + (x1[1] - x2[1]) * (x1[1] - x2[1]));
                double[] n = { (x2[1] - x1[1]) / elen, (x1[0] - x2[0]) / elen };
                double dlen = Math.Sqrt(d[0] * d[0] + d[1] * d[1]);
                if (dlen < dist)
                {
                    minDivIndex = idiv;
                    dist = dlen;
                    alpha = t;
                    ndist = (n[0] * d[0] + n[1] * d[1]) / elen;
                    normX = n[0];
                    normY = n[1];
                }
            }
        }

        private void SolveLinearStatic()
        {
            for (uint i = 0; i < No * 3; i++)
            {
                Res[i] = 0;
            }
            Mat.Clear();
            ////////////////
            uint ndiv = No - 1;
            for (uint idiv = 0; idiv < ndiv; idiv++)
            {
                double t1 = Ut[idiv * 3 + 2];
                double t2 = Ut[idiv * 3 + 5];
                double[][][][] eC = new double[2][][][]; //double[2][2][3][3];
                for (int p = 0; p < 2; p++)
                {
                    eC[p] = new double[2][][];
                    for (int q = 0; q < 2; q++)
                    {
                        eC[p][q] = new double[3][];
                        for (int r = 0; r < 3; r++)
                        {
                            eC[p][q][r] = new double[3];
                        }
                    }
                }
                double[][] eRes = new double[2][]; //double[2][3];
                for (int p = 0; p < 2; p++)
                {
                    eRes[p] = new double[3];
                }
                double[] g = { 0, 0 };
                {
                    //double* x1 = ini_x + idiv * 2;
                    //double* x2 = ini_x + idiv * 2 + 2;
                    //double* u1 = ut + idiv * 3;
                    //double* u2 = ut + idiv * 3 + 3;
                    double[] x1 = new double[2];
                    for (int i = 0; i < 2; i++)
                    {
                        x1[i] = IniX[i + idiv * 2];
                    }
                    double[] x2 = new double[2];
                    for (int i = 0; i < 2; i++)
                    {
                        x2[i] = IniX[i + idiv * 2 + 2];
                    }
                    double[] u1 = new double[2];
                    for (int i = 0; i < 2; i++)
                    {
                        u1[i] = Ut[i + idiv * 3];
                    }
                    double[] u2 = new double[2];
                    for (int i = 0; i < 2; i++)
                    {
                        u2[i] = Ut[i + idiv * 3 + 3];
                    }
                    GetCoeffMatLinear(EI, AE, ARho, g,
                                       x1, u1, t1,
                                       x2, u2, t2,
                                       ref eRes, ref eC);
                }
                Mat.Merge(idiv, eC);
                for (uint i = 0; i < 3; i++)
                {
                    Res[idiv * 3 + i] += eRes[0][i];
                    Res[idiv * 3 + 3 + i] += eRes[1][i];
                }
            }
            for (uint ino = 0; ino<No;ino++)
            {
                for (uint idof = 0; idof < 3; idof++)
                {
                    if (BCFlag[ino * 3 + idof] == 0)
                    {
                        continue;
                    }
                    Mat.FixBoundaryCondition(ino, idof);
                    Res[ino * 3 + idof] = 0;
                }
            }
            {
                double nres = 0;
                for (uint i = 0; i < No * 3; i++)
                {
                    nres += Res[i] * Res[i];
                }
            }
            Mat.ILUFrac();
            {
                for (uint i = 0; i < No * 3; i++)
                {
                    Dut[i] = Res[i];
                }
                Mat.Solve(Dut);
            }
            for (uint i = 0; i < No * 3; i++)
            {
                Ut[i] += Dut[i];
            }
        }

        // x1[2] u1[2] x2[2] u2[2] eQ[6] eK[][6]
        private void GetEMat_Linear(double EI, double AE,
                    double[] x1, double[] u1, double t1,
                    double[] x2, double[] u2, double t2,
                    ref double[] eQ, ref double[][] eK)
        {
            double eLen = Math.Sqrt((x2[0] - x1[0]) * (x2[0] - x1[0]) + (x2[1] - x1[1]) * (x2[1] - x1[1]));
            double inv_eLen = 1.0 / eLen;

            double[][] eC = new double[6][];//double[6][6];
            for (int p = 0; p < 6; p++)
            {
                eC[p] = new double[6];
            }
            {
                double tmp1 = EI / (eLen * eLen * eLen);
                double tmp2 = AE / eLen;
                eC[0][0] = tmp2; eC[0][1] = 0; eC[0][2] = 0; eC[0][3] = -tmp2; eC[0][4] = 0; eC[0][5] = 0;
                eC[1][0] = 0; eC[1][1] = tmp1 * 12; eC[1][2] = tmp1 * eLen * 6; eC[1][3] = 0; eC[1][4] = -tmp1 * 12; eC[1][5] = tmp1 * eLen * 6;
                eC[2][0] = 0; eC[2][1] = tmp1 * eLen * 6; eC[2][2] = tmp1 * eLen * eLen * 4; eC[2][3] = 0; eC[2][4] = -tmp1 * eLen * 6; eC[2][5] = tmp1 * eLen * eLen * 2;
                eC[3][0] = -tmp2; eC[3][1] = 0; eC[3][2] = 0; eC[3][3] = tmp2; eC[3][4] = 0; eC[3][5] = 0;
                eC[4][0] = 0; eC[4][1] = -tmp1 * 12; eC[4][2] = -tmp1 * eLen * 6; eC[4][3] = 0; eC[4][4] = tmp1 * 12; eC[4][5] = -tmp1 * eLen * 6;
                eC[5][0] = 0; eC[5][1] = tmp1 * eLen * 6; eC[5][2] = tmp1 * eLen * eLen * 2; eC[5][3] = 0; eC[5][4] = -tmp1 * eLen * 6; eC[5][5] = tmp1 * eLen * eLen * 4;
            }
            double[] T = { (x2[0] - x1[0]) * inv_eLen, (x2[1] - x1[1]) * inv_eLen };
            double[][] eR = new double[6][];//double[6][6];
            for (int p = 0; p < 6; p++)
            {
                eR[p] = new double[6];
            }
            {
                eR[0][0] = T[0]; eR[0][1] = -T[1]; eR[0][2] = 0; eR[0][3] = 0; eR[0][4] = 0; eR[0][5] = 0;
                eR[1][0] = T[1]; eR[1][1] = T[0]; eR[1][2] = 0; eR[1][3] = 0; eR[1][4] = 0; eR[1][5] = 0;
                eR[2][0] = 0; eR[2][1] = 0; eR[2][2] = 1; eR[2][3] = 0; eR[2][4] = 0; eR[2][5] = 0;
                eR[3][0] = 0; eR[3][1] = 0; eR[3][2] = 0; eR[3][3] = T[0]; eR[3][4] = -T[1]; eR[3][5] = 0;
                eR[4][0] = 0; eR[4][1] = 0; eR[4][2] = 0; eR[4][3] = T[1]; eR[4][4] = T[0]; eR[4][5] = 0;
                eR[5][0] = 0; eR[5][1] = 0; eR[5][2] = 0; eR[5][3] = 0; eR[5][4] = 0; eR[5][5] = 1;
            }
            for (uint i = 0; i < 6; i++)
            {
                for (uint j = 0; j < 6; j++)
                {
                    eK[i][j] = 0;
                    for (uint k = 0; k < 6; k++)
                    {
                        for (uint l = 0; l < 6; l++)
                        {
                            eK[i][j] += eR[i][k]* eC[k][l]* eR[j][l];
                        }
                    }
                }
            }
            for (uint i = 0; i < 6; i++)
            {
                eQ[i] = eK[i][0] * u1[0] + eK[i][1] * u1[1] + eK[i][2] * t1
                       + eK[i][3] * u2[0] + eK[i][4] * u2[1] + eK[i][5] * t2;
            }
        }

        // x1[2] u1[2] x2[2] u2[2] eRes[2][3] eC[2][2][3][3]
        private void GetCoeffMatLinear(double EI, double AE,
                        double ARho, double[] g,
                        double[] x1, double[] u1, double t1,
                        double[] x2, double[] u2, double t2,
                        ref double[][] eRes, ref double[][][][] eC)
        {
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    eRes[i][j] = 0;
                }
            }
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            eC[i][j][k][l] = 0;
                        }
                    }
                }
            }
            ////////////////
            double eLen = Math.Sqrt((x2[0] - x1[0]) * (x2[0] - x1[0]) + (x2[1] - x1[1]) * (x2[1] - x1[1]));
            double erho = ARho * eLen * 0.5;
            double dtmp_x = g[0] * erho;
            double dtmp_y = g[1] * erho;
            {
                eRes[0][0] += dtmp_x;
                eRes[0][1] += dtmp_y;
                eRes[1][0] += dtmp_x;
                eRes[1][1] += dtmp_y;
            }

            double[] eQ = new double[6];
            double[][] eK = new double[6][]; //double[6][6];
            for (int p = 0; p < 6; p++)
            {
                eK[p] = new double[6];
            }
            GetEMat_Linear(EI, AE, x1, u1, t1, x2, u2, t2, ref eQ, ref eK);
            for (uint i = 0; i < 3; i++)
            {
                for (uint j = 0; j < 3; j++)
                {
                    eC[0][0][i][j] += eK[0 + i][0 + j];
                    eC[0][1][i][j] += eK[0 + i][3 + j];
                    eC[1][0][i][j] += eK[3 + i][0 + j];
                    eC[1][1][i][j] += eK[3 + i][3 + j];
                }
            }
            for (uint i = 0; i < 3; i++)
            {
                eRes[0][i] -= eQ[0 + i];
                eRes[1][i] -= eQ[3 + i];
            }
        }
    }
}
