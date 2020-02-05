using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TriangleFEBellInterpolate : IInterpolate
    {
        public TriangleFE Owner { get; set; }

        public TriangleFEBellInterpolate()
        {

        }

        public TriangleFEBellInterpolate(TriangleFE owner)
        {
            Owner = owner;
        }


        public uint GetNodeCount()
        {
            // 3節点 x 6変数
            return 3;
        }

        protected static int[][] NodeIdsForQuantity = new int[6][]
        {
            new int[] {0, 6, 12},  // φ
            new int[] {1, 7, 13},  // φx
            new int[] {2, 8, 14},  // φy
            new int[] {3, 9, 15},  // φxx
            new int[] {4, 10, 16}, // φxy
            new int[] {5, 11, 17}  // φyy
        };

        protected int[] GetNodeIdsForQuantity()
        {
            System.Diagnostics.Debug.Assert(Owner.QuantityId >= 1 && Owner.QuantityId <= 6);
            System.Diagnostics.Debug.Assert(Owner.QuantityIdBaseOffset == 1);
            int quantityId = Owner.QuantityId - Owner.QuantityIdBaseOffset;
            System.Diagnostics.Debug.Assert(quantityId >= 0 && quantityId < 6);
            return NodeIdsForQuantity[quantityId];
        }

        protected double[] GetValuesForQuantity(double[] values)
        {
            System.Diagnostics.Debug.Assert(values.Length == 18);
            int[] nodeIds = GetNodeIdsForQuantity();
            double[] valuesForQuantity = new double[3];
            for (int i = 0; i < 3; i++)
            {
                int nodeId = nodeIds[i];
                valuesForQuantity[i] = values[nodeId];
            }
            return valuesForQuantity;
        }

        public double[] GetNodeL(int nodeId)
        {
            int[] nodeIndexs = GetNodeIdsForQuantity();
            int _nodeId = nodeIndexs[nodeId];

            double[] _nodeL = _GetNodeL(_nodeId);
            double[] nodeL = GetValuesForQuantity(_nodeL);
            return nodeL;
        }

        protected double[] _GetNodeL(int nodeId)
        {
            double[][] nodeL = new double[18][];

            for (int i = 0; i < 6; i++)
            {
                nodeL[i] = new double[3] { 1.0, 0.0, 0.0 };
                nodeL[i + 6] = new double[3] { 0.0, 1.0, 0.0 };
                nodeL[i + 12] = new double[3] { 0.0, 0.0, 1.0 };
            }
            return nodeL[nodeId];
        }

        protected double[] _GetGlobalCoordShapeFunction(double[] barN)
        {
            double[] a;
            double[] b;
            double[] c;
            Owner.CalcTransMatrix(out a, out b, out c);
            // 1/(2A)の係数をなくしたa, b, cを求める
            double A = Owner.GetArea();
            for (int i = 0; i < 3; i++)
            {
                a[i] *= (2.0 * A);
                b[i] *= (2.0 * A);
                c[i] *= (2.0 * A);
            }

            IvyFEM.Lapack.DoubleMatrix C0 = new IvyFEM.Lapack.DoubleMatrix(6, 6);
            C0[0, 0] = 1.0;
            C0[1, 1] = c[1];
            C0[1, 2] = -b[1];
            C0[2, 1] = -c[0];
            C0[2, 2] = b[0];
            C0[3, 3] = c[1] * c[1];
            C0[3, 4] = -2.0 * b[1] * c[1];
            C0[3, 5] = b[1] * b[1];
            C0[4, 3] = -c[0] * c[1];
            C0[4, 4] = b[0] * c[1] + b[1] * c[0];
            C0[4, 5] = -b[0] * b[1];
            C0[5, 3] = c[0] * c[0];
            C0[5, 4] = -2.0 * b[0] * c[0];
            C0[5, 5] = b[0] * b[0];

            C0.Transpose();

            double[] barN1 = new double[6];
            double[] barN2 = new double[6];
            double[] barN3 = new double[6];
            Array.Copy(barN, 0, barN1, 0, barN1.Length);
            Array.Copy(barN, 6, barN2, 0, barN2.Length);
            Array.Copy(barN, 12, barN3, 0, barN3.Length);

            double[] N1 = C0 * barN1;
            double[] N2 = C0 * barN2;
            double[] N3 = C0 * barN3;

            double[] N = new double[18];
            N1.CopyTo(N, 0);
            N2.CopyTo(N, 6);
            N3.CopyTo(N, 12);
            return N;
        }

        public double[] CalcN(double[] L)
        {
            double[] _barN = _CalcN(L);
            double[] _N = _GetGlobalCoordShapeFunction(_barN);
            double[] N = GetValuesForQuantity(_N);
            return N;
        }

        protected double[] _CalcN(double[] L)
        {
            double[] N = new double[18];

            double r = L[0];
            double s = L[1];
            double[] p = new double[21];
            p[0] = 1.0;
            p[1] = r;
            p[2] = s;
            p[3] = r * r;
            p[4] = r * s;
            p[5] = s * s;
            p[6] = r * r * r;
            p[7] = r * r * s;
            p[8] = r * s * s;
            p[9] = s * s * s;
            p[10] = r * r * r * r;
            p[11] = r * r * r * s;
            p[12] = r * r * s * s;
            p[13] = r * s * s * s;
            p[14] = s * s * s * s;
            p[15] = r * r * r * r * r;
            p[16] = r * r * r * r * s;
            p[17] = r * r * r * s * s;
            p[18] = r * r * s * s * s;
            p[19] = r * s * s * s * s;
            p[20] = s * s * s * s * s;

            double[][] alphas = _CalcAlphas();

            for (int iN = 0; iN < 18; iN++)
            {
                for (int iAlpha = 0; iAlpha < 21; iAlpha++)
                {
                    double[] alpha = alphas[iAlpha];
                    N[iN] += alpha[iN] * p[iAlpha];
                }
            }
            return N;
        }

        protected double[][] _CalcAlphas()
        {
            double[][] alphas = new double[21][];
            double[] a;
            double[] b;
            double[] c;
            Owner.CalcTransMatrix(out a, out b, out c);
            // 1/(2A)の係数をなくしたa, b, cを求める
            double A = Owner.GetArea();
            for (int i = 0; i < 3; i++)
            {
                a[i] *= (2.0 * A);
                b[i] *= (2.0 * A);
                c[i] *= (2.0 * A);
            }

            // 計算順番に依存する
            int[] indexs = new int[21] {
                0, 1, 2, 3, 4, 5,
                6, 10, 15,
                9, 14, 20,
                8, 13, 19,
                7, 11, 16,
                12, 17, 18
            };
            for (int i = 0; i < indexs.Length; i++)
            {
                int index = indexs[i];
                double[] coef = _CalcAlpha(index, alphas, a, b, c);
                alphas[index] = coef;
            }
            return alphas;
        }

        protected double[] _CalcAlpha(int index, double[][] alphas, double[] a, double[] b, double[] c)
        {
            double[] coef = new double[18];

            // 式(19)
            if (index == 0)
            {
                coef[12] = 1.0;
            }
            else if (index == 1)
            {
                coef[13] = 1.0;
            }
            else if (index == 2)
            {
                coef[14] = 1.0;
            }
            else if (index == 3)
            {
                coef[15] = 1.0 / 2.0;
            }
            else if (index == 4)
            {
                coef[16] = 1.0;
            }
            else if (index == 5)
            {
                coef[17] = 1.0 / 2.0;
            }
            // 式(20)
            else if (index == 6)
            {
                coef[0] = 10.0;
                coef[1] = -4.0;
                coef[3] = 1.0 / 2.0;
                coef[12] = -10.0;
                coef[13] = -6.0;
                coef[15] = -3.0 / 2.0;

            }
            // 式(21)
            else if (index == 10)
            {
                coef[0] = -15.0;
                coef[1] = 7.0;
                coef[3] = -1.0;
                coef[12] = 15.0;
                coef[13] = 8.0;
                coef[15] = 3.0 / 2.0;
            }
            // 式(22)
            else if (index == 15)
            {
                coef[0] = 6.0;
                coef[1] = -3.0;
                coef[3] = 1.0 / 2.0;
                coef[12] = -6.0;
                coef[13] = -3.0;
                coef[15] = -1.0 / 2.0;
            }
            // 式(23)
            else if (index == 9)
            {
                coef[6] = 10.0;
                coef[8] = -4.0;
                coef[11] = 1.0 / 2.0;
                coef[12] = -10.0;
                coef[14] = -6.0;
                coef[17] = -3.0 / 2.0;
            }
            // 式(24)
            else if (index == 14)
            {
                coef[6] = -15.0;
                coef[8] = 7.0;
                coef[11] = -1.0;
                coef[12] = 15.0;
                coef[14] = 8.0;
                coef[17] = 3.0 / 2.0;
            }
            // 式(25)
            else if (index == 20)
            {
                coef[6] = 6.0;
                coef[8] = -3.0;
                coef[11] = 1.0 / 2.0;
                coef[12] = -6.0;
                coef[14] = -3.0;
                coef[17] = -1.0 / 2.0;
            }
            // 式(34)
            else if (index == 8)
            {
                double[] k1 = _CalcK(0, alphas, a, b, c);
                coef[7] = -5.0;
                coef[10] = 1.0;
                double[] alpha1 = alphas[1];
                double[] alpha4 = alphas[4];
                for (int i = 0; i < coef.Length; i++)
                {
                    coef[i] += 5.0 * alpha1[i];
                    coef[i] += 4.0 * alpha4[i];
                    coef[i] += k1[i];
                }
            }
            // 式(35)
            else if (index == 13)
            {
                double[] k1 = _CalcK(0, alphas, a, b, c);
                coef[7] = 14.0;
                coef[10] = -3.0;
                double[] alpha1 = alphas[1];
                double[] alpha4 = alphas[4];
                for (int i = 0; i < coef.Length; i++)
                {
                    coef[i] += -14.0 * alpha1[i];
                    coef[i] += -11.0 * alpha4[i];
                    coef[i] += -2.0 * k1[i];
                }
            }
            // 式(36)
            else if (index == 19)
            {
                double[] k1 = _CalcK(0, alphas, a, b, c);
                coef[7] = -8.0;
                coef[10] = 2.0;
                double[] alpha1 = alphas[1];
                double[] alpha4 = alphas[4];
                for (int i = 0; i < coef.Length; i++)
                {
                    coef[i] += 8.0 * alpha1[i];
                    coef[i] += 6.0 * alpha4[i];
                    coef[i] += k1[i];
                }
            }
            // 式(41)
            else if (index == 7)
            {
                double[] k2 = _CalcK(1, alphas, a, b, c);
                coef[2] = -5.0;
                coef[4] = 1.0;
                double[] alpha2 = alphas[2];
                double[] alpha4 = alphas[4];
                for (int i = 0; i < coef.Length; i++)
                {
                    coef[i] += 5.0 * alpha2[i];
                    coef[i] += 4.0 * alpha4[i];
                    coef[i] += k2[i];
                }
            }
            // 式(42)
            else if (index == 11)
            {
                double[] k2 = _CalcK(1, alphas, a, b, c);
                coef[2] = 14.0;
                coef[4] = -3.0;
                double[] alpha2 = alphas[2];
                double[] alpha4 = alphas[4];
                for (int i = 0; i < coef.Length; i++)
                {
                    coef[i] += -14.0 * alpha2[i];
                    coef[i] += -11.0 * alpha4[i];
                    coef[i] += -2.0 * k2[i];
                }
            }
            // 式(43)
            else if (index == 16)
            {
                double[] k2 = _CalcK(1, alphas, a, b, c);
                coef[2] = -8.0;
                coef[4] = 2.0;
                double[] alpha2 = alphas[2];
                double[] alpha4 = alphas[4];
                for (int i = 0; i < coef.Length; i++)
                {
                    coef[i] += 8.0 * alpha2[i];
                    coef[i] += 6.0 * alpha4[i];
                    coef[i] += k2[i];
                }
            }
            // 式(47)
            else if (index == 12)
            {
                double[] k3 = _CalcK(2, alphas, a, b, c);
                double det = -4.0;
                coef[5] = (1.0 / det) * (-6.0);
                coef[9] = (1.0 / det) * (-4.0);
                double[] alpha3 = alphas[3];
                double[] alpha5 = alphas[5];
                double[] alpha7 = alphas[7];
                double[] alpha8 = alphas[8];
                for (int i = 0; i < coef.Length; i++)
                {
                    coef[i] += (1.0 / det) * 8.0 * alpha3[i];
                    coef[i] += (1.0 / det) * 12.0 * alpha5[i];
                    coef[i] += (1.0 / det) * 8.0 * alpha7[i];
                    coef[i] += (1.0 / det) * 12.0 * alpha8[i];
                    coef[i] += (1.0 / det) * 4.0 * k3[i];
                }
            }
            // 式(48)
            else if (index == 17)
            {
                double[] k3 = _CalcK(2, alphas, a, b, c);
                double det = -4.0;
                coef[5] = (1.0 / det) * 4.0;
                coef[9] = (1.0 / det) * 4.0;
                double[] alpha3 = alphas[3];
                double[] alpha5 = alphas[5];
                double[] alpha7 = alphas[7];
                double[] alpha8 = alphas[8];
                for (int i = 0; i < coef.Length; i++)
                {
                    coef[i] += (1.0 / det) * (-8.0) * alpha3[i];
                    coef[i] += (1.0 / det) * (-8.0) * alpha5[i];
                    coef[i] += (1.0 / det) * (-8.0) * alpha7[i];
                    coef[i] += (1.0 / det) * (-8.0) * alpha8[i];
                    coef[i] += (1.0 / det) * (-4.0) * k3[i];
                }
            }
            // 式(49)
            else if (index == 18)
            {
                double[] k3 = _CalcK(2, alphas, a, b, c);
                double det = -4.0;
                coef[5] = (1.0 / det) * 6.0;
                coef[9] = (1.0 / det) * 2.0;
                double[] alpha3 = alphas[3];
                double[] alpha5 = alphas[5];
                double[] alpha7 = alphas[7];
                double[] alpha8 = alphas[8];
                for (int i = 0; i < coef.Length; i++)
                {
                    coef[i] += (1.0 / det) * (-4.0) * alpha3[i];
                    coef[i] += (1.0 / det) * (-12.0) * alpha5[i];
                    coef[i] += (1.0 / det) * (-4.0) * alpha7[i];
                    coef[i] += (1.0 / det) * (-12.0) * alpha8[i];
                    coef[i] += (1.0 / det) * (-4.0) * k3[i];
                }
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return coef;
        }

        protected double[] _CalcK(int KIndex, double[][] alphas, double[] a, double[] b, double[] c)
        {
            double[] coef = new double[18];

            // 式(33)
            if (KIndex == 0)
            {
                coef[7] = 8.0;
                coef[13] = 8.0;
                coef[10] = -2.0;
                coef[16] = 2.0;
                double d01 = b[0] * b[1] + c[0] * c[1];
                double l0 = Math.Sqrt(b[0] * b[0] + c[0] * c[0]);
                double workC = d01 / (l0 * l0);
                coef[8] = workC * 8.0;
                coef[14] = workC * 8.0;
                coef[11] = workC * (-2.0);
                coef[17] = workC * 2.0;

                double[] alpha2 = alphas[2];
                double[] alpha5 = alphas[5];
                double[] alpha9 = alphas[9];
                double[] alpha14 = alphas[14];
                double[] alpha20 = alphas[20];
                for (int i = 0; i < coef.Length; i++)
                {
                    coef[i] += workC * (-16.0) * alpha2[i];
                    coef[i] += workC * (-16.0) * alpha5[i];
                    coef[i] += workC * (-12.0) * alpha9[i];
                    coef[i] += workC * (-8.0) * alpha14[i];
                    coef[i] += workC * (-5.0) * alpha20[i];
                }
            }
            // 式(40)
            else if (KIndex == 1)
            {
                coef[14] = 8.0;
                coef[2] = 8.0;
                coef[16] = 2.0;
                coef[4] = -2.0;
                double d01 = b[0] * b[1] + c[0] * c[1];
                double l1 = Math.Sqrt(b[1] * b[1] + c[1] * c[1]);
                double workC = d01 / (l1 * l1);
                coef[13] = workC * 8.0;
                coef[1] = workC * 8.0;
                coef[15] = workC * 2.0;
                coef[3] = workC * (-2.0);
                double[] alpha1 = alphas[1];
                double[] alpha3 = alphas[3];
                double[] alpha6 = alphas[6];
                double[] alpha10 = alphas[10];
                double[] alpha15 = alphas[15];
                for (int i = 0; i < coef.Length; i++)
                {
                    coef[i] += workC * (-16.0) * alpha1[i];
                    coef[i] += workC * (-16.0) * alpha3[i];
                    coef[i] += workC * (-12.0) * alpha6[i];
                    coef[i] += workC * (-8.0) * alpha10[i];
                    coef[i] += workC * (-5.0) * alpha15[i];
                }
            }
            // 式(46)
            else if (KIndex == 2)
            {
                double d02 = b[0] * b[2] + c[0] * c[2];
                double d12 = b[1] * b[2] + c[1] * c[2];
                double d02d12 = d02 + d12;
                double workC02 = d02 / d02d12;
                double workC12 = d12 / d02d12;
                coef[1] = workC02 * 8.0;
                coef[7] = workC02 * 8.0;
                coef[3] = workC02 * (-2.0);
                coef[4] = workC02 * 2.0;
                coef[9] = workC02 * 2.0;
                coef[10] = workC02 * (-1.0);
                coef[2] = workC12 * 8.0;
                coef[8] = workC12 * 8.0;
                coef[4] += workC12 * (-2.0);
                coef[5] = workC12 * 2.0 * 5.0 / 4.0;
                coef[9] += workC12 * 2.0 * (-1.0 / 4.0);
                coef[10] += workC12 * 2.0;
                coef[11] = workC12 * 2.0 * (-1.0);
                double[] alpha1 = alphas[1];
                double[] alpha3 = alphas[3];
                double[] alpha4 = alphas[4];
                double[] alpha6 = alphas[6];
                double[] alpha7 = alphas[7];
                double[] alpha8 = alphas[8];
                double[] alpha10 = alphas[10];
                double[] alpha11 = alphas[11];
                double[] alpha13 = alphas[13];
                double[] alpha15 = alphas[15];
                double[] alpha16 = alphas[16];
                double[] alpha19 = alphas[19];
                double[] alpha2 = alphas[2];
                double[] alpha5 = alphas[5];
                double[] alpha9 = alphas[9];
                double[] alpha14 = alphas[14];
                double[] alpha20 = alphas[20];
                for (int i = 0; i < coef.Length; i++)
                {
                    coef[i] += workC02 * (-16.0) * alpha1[i];
                    coef[i] += workC02 * (-16.0) * alpha3[i];
                    coef[i] += workC02 * (-8.0) * alpha4[i];
                    coef[i] += workC02 * (-12.0) * alpha6[i];
                    coef[i] += workC02 * (-8.0) * alpha7[i];
                    coef[i] += workC02 * (-4.0) * alpha8[i];
                    coef[i] += workC02 * (-8.0) * alpha10[i];
                    coef[i] += workC02 * (-6.0) * alpha11[i];
                    coef[i] += workC02 * (-2.0) * alpha13[i];
                    coef[i] += workC02 * (-5.0) * alpha15[i];
                    coef[i] += workC02 * (-4.0) * alpha16[i];
                    coef[i] += workC02 * (-1.0) * alpha19[i];

                    coef[i] += workC12 * (-16.0) * alpha2[i];
                    coef[i] += workC12 * 1.0 * alpha3[i];
                    coef[i] += workC12 * (-8.0) * alpha4[i];
                    coef[i] += workC12 * (-17.0) * alpha5[i];
                    coef[i] += workC12 * (-3.0) * alpha7[i];
                    coef[i] += workC12 * (-9.0) * alpha8[i];
                    coef[i] += workC12 * (-12.0) * alpha9[i];
                    coef[i] += workC12 * (-2.0) * alpha11[i];
                    coef[i] += workC12 * (-6.0) * alpha13[i];
                    coef[i] += workC12 * (-8.0) * alpha14[i];
                    coef[i] += workC12 * (-1.0) * alpha16[i];
                    coef[i] += workC12 * (-4.0) * alpha19[i];
                    coef[i] += workC12 * (-5.0) * alpha20[i];
                }
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }

            return coef;
        }

        public double[][] CalcNu(double[] L)
        {
            double[][] _barNu = _CalcNu(L);
            double[][] _Nu = new double[2][];
            for (int i = 0; i < 2; i++)
            {
                _Nu[i] = _GetGlobalCoordShapeFunction(_barNu[i]);
            }
            double[][] Nu = new double[2][];
            for (int i = 0; i < 2; i++)
            {
                Nu[i] = GetValuesForQuantity(_Nu[i]);
            }
            return Nu;
        }

        protected double[][] _CalcNu(double[] L)
        {
            double[][] Nu = new double[2][];

            double[] a;
            double[] b;
            double[] c;
            Owner.CalcTransMatrix(out a, out b, out c);

            double[] Nx = _CalcNx(L, a, b, c);
            Nu[0] = Nx;

            double[] Ny = _CalcNy(L, a, b, c);
            Nu[1] = Ny;

            return Nu;
        }

        protected double[] _CalcNx(double[] L, double[] a, double[] b, double[] c)
        {
            double[] Nx = new double[18];

            double r = L[0];
            double s = L[1];
            double[] p = new double[21];
            p[0] = 0;
            p[1] = b[0];
            p[2] = b[1];
            p[3] = 2.0 * r * b[0];
            p[4] = s * b[0] + r * b[1];
            p[5] = 2.0 * s * b[1];
            p[6] = 3.0 * r * r * b[0];
            p[7] = 2.0 * r * s * b[0] + r * r * b[1];
            p[8] = s * s * b[0] + 2.0 * r * s * b[1];
            p[9] = 3.0 * s * s * b[1];
            p[10] = 4.0 * r * r * r * b[0];
            p[11] = 3.0 * r * r * s * b[0] + r * r * r * b[1];
            p[12] = 2.0 * r * s * s * b[0] + 2.0 * r * r * s * b[1];
            p[13] = s * s * s * b[0] + 3.0 * r * s * s * b[1];
            p[14] = 4.0 * s * s * s * b[1];
            p[15] = 5.0 * r * r * r * r * b[0];
            p[16] = 4.0 * r * r * r * s * b[0] + r * r * r * r * b[1];
            p[17] = 3.0 * r * r * s * s * b[0] + 2.0 * r * r * r * s * b[1];
            p[18] = 2.0 * r * s * s * s * b[0] + 3.0 * r * r * s * s * b[1];
            p[19] = s * s * s * s * b[0] + 4.0 * r * s * s * s * b[1];
            p[20] = 5.0 * s * s * s * s * b[1];

            double[][] alphas = _CalcAlphas();

            for (int iN = 0; iN < 18; iN++)
            {
                for (int iAlpha = 0; iAlpha < 21; iAlpha++)
                {
                    double[] alpha = alphas[iAlpha];
                    Nx[iN] += alpha[iN] * p[iAlpha];
                }
            }
            return Nx;
        }

        protected double[] _CalcNy(double[] L, double[] a, double[] b, double[] c)
        {
            double[] Ny = new double[18];

            double r = L[0];
            double s = L[1];
            double[] p = new double[21];
            p[0] = 0;
            p[1] = c[0];
            p[2] = c[1];
            p[3] = 2.0 * r * c[0];
            p[4] = s * c[0] + r * c[1];
            p[5] = 2.0 * s * c[1];
            p[6] = 3.0 * r * r * c[0];
            p[7] = 2.0 * r * s * c[0] + r * r * c[1];
            p[8] = s * s * c[0] + 2.0 * r * s * c[1];
            p[9] = 3.0 * s * s * c[1];
            p[10] = 4.0 * r * r * r * c[0];
            p[11] = 3.0 * r * r * s * c[0] + r * r * r * c[1];
            p[12] = 2.0 * r * s * s * c[0] + 2.0 * r * r * s * c[1];
            p[13] = s * s * s * c[0] + 3.0 * r * s * s * c[1];
            p[14] = 4.0 * s * s * s * c[1];
            p[15] = 5.0 * r * r * r * r * c[0];
            p[16] = 4.0 * r * r * r * s * c[0] + r * r * r * r * c[1];
            p[17] = 3.0 * r * r * s * s * c[0] + 2.0 * r * r * r * s * c[1];
            p[18] = 2.0 * r * s * s * s * c[0] + 3.0 * r * r * s * s * c[1];
            p[19] = s * s * s * s * c[0] + 4.0 * r * s * s * s * c[1];
            p[20] = 5.0 * s * s * s * s * c[1];

            double[][] alphas = _CalcAlphas();

            for (int iN = 0; iN < 18; iN++)
            {
                for (int iAlpha = 0; iAlpha < 21; iAlpha++)
                {
                    double[] alpha = alphas[iAlpha];
                    Ny[iN] += alpha[iN] * p[iAlpha];
                }
            }
            return Ny;
        }

        public double[,][] CalcNuv(double[] L)
        {
            double[,][] _barNuv = _CalcNuv(L);
            double[,][] _Nuv = new double[2, 2][];
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    _Nuv[i, j] = _GetGlobalCoordShapeFunction(_barNuv[i, j]);
                }
            }

            double[,][] Nuv = new double[2, 2][];
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    Nuv[i, j] = GetValuesForQuantity(_Nuv[i, j]);
                }
            }
            return Nuv;
        }

        protected double[,][] _CalcNuv(double[] L)
        {
            double[,][] Nuv = new double[2, 2][];

            double[] a;
            double[] b;
            double[] c;
            Owner.CalcTransMatrix(out a, out b, out c);

            double[] Nxx = _CalcNxx(L, a, b, c);
            Nuv[0, 0] = Nxx;

            double[] Nyy = _CalcNyy(L, a, b, c);
            Nuv[1, 1] = Nyy;

            double[] Nxy = _CalcNxy(L, a, b, c);
            Nuv[0, 1] = Nxy;
            Nuv[1, 0] = Nxy;

            return Nuv;
        }

        protected double[] _CalcNxx(double[] L, double[] a, double[] b, double[] c)
        {
            double[] Nxx = new double[18];

            double r = L[0];
            double s = L[1];
            double[] p = new double[21];
            p[0] = 0;
            p[1] = 0;
            p[2] = 0;
            p[3] = 2.0 * b[0] * b[0];
            p[4] = 2.0 * b[0] * b[1];
            p[5] = 2.0 * b[1] * b[1];
            p[6] = 6.0 * r * b[0] * b[0];
            p[7] = 4.0 * r * b[0] * b[1] + 2.0 * s * b[0] * b[0];
            p[8] = 4.0 * s * b[0] * b[1] + 2.0 * r * b[1] * b[1];
            p[9] = 6.0 * s * b[1] * b[1];
            p[10] = 12.0 * r * r * b[0] * b[0];
            p[11] = 6.0 * r * s * b[0] * b[0] + 6.0 * r * r * b[0] * b[1];
            p[12] = 2.0 * s * s * b[0] * b[0] + 8.0 * r * s * b[0] * b[1] + 2.0  * r * r * b[1] * b[1];
            p[13] = 6.0 * s * s * b[0] * b[1] + 6.0 * r * s * b[1] * b[1];
            p[14] = 12.0 * s * s * b[1] * b[1];
            p[15] = 20.0 * r * r * r * b[0] * b[0];
            p[16] = 12.0 * r * r * s * b[0] * b[0] + 8.0 * r * r * r * b[0] * b[1];
            p[17] = 6.0 * r * s * s * b[0] * b[0] + 12.0 * r * r * s * b[0] * b[1] + 2.0 * r * r * r * b[1] * b[1];
            p[18] = 2.0 * s * s * s * b[0] * b[0] + 12.0 * r * s * s * b[0] * b[1] + 6.0 * r * r * s * b[1] * b[1];
            p[19] = 8.0 * s * s * s * b[0] * b[1] + 12.0 * r * s * s * b[1] * b[1];
            p[20] = 20.0 * s * s * s * b[1] * b[1];

            double[][] alphas = _CalcAlphas();

            for (int iN = 0; iN < 18; iN++)
            {
                for (int iAlpha = 0; iAlpha < 21; iAlpha++)
                {
                    double[] alpha = alphas[iAlpha];
                    Nxx[iN] += alpha[iN] * p[iAlpha];
                }
            }
            return Nxx;
        }

        protected double[] _CalcNyy(double[] L, double[] a, double[] b, double[] c)
        {
            double[] Nyy = new double[18];

            double r = L[0];
            double s = L[1];
            double[] p = new double[21];
            p[0] = 0;
            p[1] = 0;
            p[2] = 0;
            p[3] = 2.0 * c[0] * c[0];
            p[4] = 2.0 * c[0] * c[1];
            p[5] = 2.0 * c[1] * c[1];
            p[6] = 6.0 * r * c[0] * c[0];
            p[7] = 2.0 * s * c[0] * c[0] + 4.0 * r * c[0] * c[1];
            p[8] = 4.0 * s * c[0] * c[1] + 2.0 * r * c[1] * c[1];
            p[9] = 6.0 * s * c[1] * c[1];
            p[10] = 12.0 * r * r * c[0] * c[0];
            p[11] = 6.0 * r * s * c[0] * c[0] + 6.0 * r * r * c[0] * c[1];
            p[12] = 2.0 * s * s * c[0] * c[0] + 8.0 * r * s * c[0] * c[1] + 2.0 * r * r * c[1] * c[1];
            p[13] = 6.0 * s * s * c[0] * c[1] + 6.0 * r * s * c[1] * c[1];
            p[14] = 12.0 * s * s * c[1] * c[1];
            p[15] = 20.0 * r * r * r * c[0] * c[0];
            p[16] = 12.0 * r * r * s * c[0] * c[0] + 8.0 * r * r * r * c[0] * c[1];
            p[17] = 6.0 * r * s * s * c[0] * c[0] + 12.0 * r * r * s * c[0] * c[1] + 2.0 * r * r * r * c[1] * c[1];
            p[18] = 2.0 * s * s * s * c[0] * c[0] + 12.0 * r * s * s * c[0] * c[1] + 6.0 * r * r * s * c[1] * c[1];
            p[19] = 8.0 * s * s * s * c[0] * c[1] + 12.0 * r * s * s * c[1] * c[1];
            p[20] = 20.0 * s * s * s * c[1] * c[1];

            double[][] alphas = _CalcAlphas();

            for (int iN = 0; iN < 18; iN++)
            {
                for (int iAlpha = 0; iAlpha < 21; iAlpha++)
                {
                    double[] alpha = alphas[iAlpha];
                    Nyy[iN] += alpha[iN] * p[iAlpha];
                }
            }
            return Nyy;
        }

        protected double[] _CalcNxy(double[] L, double[] a, double[] b, double[] c)
        {
            double[] Nxy = new double[18];

            double r = L[0];
            double s = L[1];
            double[] p = new double[21];
            p[0] = 0;
            p[1] = 0;
            p[2] = 0;
            p[3] = 2.0 * b[0] * c[0];
            p[4] = b[1] * c[0] + b[0] * c[1];
            p[5] = 2.0 * b[1] * c[1];
            p[6] = 6.0 * r * b[0] * c[0];
            p[7] = 2.0 * s * b[0] * c[0] + 2.0 * r * (b[1] * c[0] + b[0] * c[1]);
            p[8] = 2.0 * s * (b[1] * c[0] + b[0] * c[1]) + 2.0 * r * b[1] * c[1];
            p[9] = 6.0 * s * b[1] * c[1];
            p[10] = 12.0 * r * r * b[0] * c[0];
            p[11] = 6.0 * r * s * b[0] * c[0] + 3.0 * r * r * (b[1] * c[0] + b[0] * c[1]);
            p[12] = 2.0 * s * s * b[0] * c[0] + 4.0 * r * s * (b[1] * c[0] + b[0] * c[1]) + 2.0 * r * r * b[1] * c[1];
            p[13] = 3.0 * s * s * (b[1] * c[0] + b[0] * c[1]) + 6.0 * r * s * b[1] * c[1];
            p[14] = 12.0 * s * s * b[1] * c[1];
            p[15] = 20.0 * r * r * r * b[0] * c[0];
            p[16] = 12.0 * r * r * s * b[0] * c[0] + 4.0 * r * r * r * (b[1] * c[0] + b[0] * c[1]);
            p[17] = 6.0 * r * s * s * b[0] * c[0] + 6.0 * r * r * s * (b[1] * c[0] + b[0] * c[1]) +
                2.0 * r * r * r * b[1] * c[1];
            p[18] = 2.0 * s * s * s * b[0] * c[0] + 6.0 * r * s * s * (b[1] * c[0] + b[0] * c[1]) +
                6.0 * r * r * s * b[1] * c[1];
            p[19] = 4.0 * s * s * s * (b[1] * c[0] + b[0] * c[1]) + 12.0 * r * s * s * b[1] * c[1];
            p[20] = 20.0 * s * s * s * b[1] * c[1];

            double[][] alphas = _CalcAlphas();

            for (int iN = 0; iN < 18; iN++)
            {
                for (int iAlpha = 0; iAlpha < 21; iAlpha++)
                {
                    double[] alpha = alphas[iAlpha];
                    Nxy[iN] += alpha[iN] * p[iAlpha];
                }
            }
            return Nxy;
        }

        public double[] CalcSN()
        {
            // 対応しない
            throw new NotImplementedException();
        }

        public double[,] CalcSNN()
        {
            // 対応しない
            throw new NotImplementedException();
        }

        public double[,][,] CalcSNuNv()
        {
            // 対応しない
            throw new NotImplementedException();
        }

        public double[][,] CalcSNuN()
        {
            // 対応しない
            throw new NotImplementedException();
        }
    }
}
