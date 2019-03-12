using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class IntegrationPoints
    {
        public int PointCount { get; set; } = 0;
        public double[][] Ls { get; set; } = null;
        public double[] Weights { get; set; } = null;

        public IntegrationPoints()
        {

        }

        public IntegrationPoints(IntegrationPoints src)
        {
            PointCount = src.PointCount;
            Ls = null;
            if (src.Ls != null)
            {
                Ls = new double[src.Ls.Length][];
                for (int i = 0; i < src.Ls.Length; i++)
                {
                    double[] srcPoint = src.Ls[i];
                    Ls[i] = new double[srcPoint.Length];
                    srcPoint.CopyTo(Ls[i], 0);
                }
            }
            Weights = null;
            if (src.Weights != null)
            {
                Weights = new double[src.Weights.Length];
                src.Weights.CopyTo(Weights, 0);
            }
        }

        public static IntegrationPoints[] LineIntegrationPoints =
        {
            new IntegrationPoints{
                PointCount = (int)LineIntegrationPointCount.Point1,
                Ls = new double[(int)LineIntegrationPointCount.Point1][]
                {
                    LineFE.GetLFromXi(0.0)
                },
                Weights = new double[(int)LineIntegrationPointCount.Point1]
                {
                    2.0
                }
            },
            new IntegrationPoints{
                PointCount = (int)LineIntegrationPointCount.Point2,
                Ls = new double[(int)LineIntegrationPointCount.Point2][]
                {
                    LineFE.GetLFromXi(-0.57735027),
                    LineFE.GetLFromXi(0.57735027)
                },
                Weights = new double[(int)LineIntegrationPointCount.Point2]
                {
                    1.0,
                    1.0
                }
            },
            new IntegrationPoints{
                PointCount = (int)LineIntegrationPointCount.Point3,
                Ls = new double[(int)LineIntegrationPointCount.Point3][]
                {
                    LineFE.GetLFromXi(-0.77459667),
                    LineFE.GetLFromXi(0.0),
                    LineFE.GetLFromXi(0.77459667)
                },
                Weights = new double[(int)LineIntegrationPointCount.Point3]
                {
                    0.55555556,
                    0.88888889,
                    0.55555556
                }
            },
            new IntegrationPoints{
                PointCount = (int)LineIntegrationPointCount.Point4,
                Ls = new double[(int)LineIntegrationPointCount.Point4][]
                {
                    LineFE.GetLFromXi(-0.86113631),
                    LineFE.GetLFromXi(-0.33998104),
                    LineFE.GetLFromXi(0.33998104),
                    LineFE.GetLFromXi(0.86113631)
                },
                Weights = new double[(int)LineIntegrationPointCount.Point4]
                {
                    0.34785485,
                    0.65214515,
                    0.65214515,
                    0.34785485
                }
            },
            new IntegrationPoints{
                PointCount = (int)LineIntegrationPointCount.Point5,
                Ls = new double[(int)LineIntegrationPointCount.Point5][]
                {
                    LineFE.GetLFromXi(-0.90617985),
                    LineFE.GetLFromXi(-0.53846931),
                    LineFE.GetLFromXi(0.0),
                    LineFE.GetLFromXi(0.53846931),
                    LineFE.GetLFromXi(0.90617985),
                },
                Weights = new double[(int)LineIntegrationPointCount.Point5]
                {
                    0.23692689,
                    0.47862867,
                    0.56888889,
                    0.47862867,
                    0.23692689
                }
            }
        };

        // α, β, γ
        private static double[] TriangleIP3Alpha = { 1.0 / 3.0, 0.6, 0.2 };
        // α, β, γ, δ, ε
        private static double[] TriangleIP7Alpha = 
            { 1.0 / 3.0, 0.05971587, 0.47014206, 0.79742669, 0.10128651};

        public static IntegrationPoints[] TriangleIntegrationPoints =
        {
            new IntegrationPoints{
                PointCount = (int)TriangleIntegrationPointCount.Point1,
                Ls = new double[(int)TriangleIntegrationPointCount.Point1][]
                {
                    new double[3] {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}
                },
                Weights = new double[(int)TriangleIntegrationPointCount.Point1]
                {
                    1.0
                }
            },
            new IntegrationPoints{
                PointCount = (int)TriangleIntegrationPointCount.Point3,
                Ls = new double[(int)TriangleIntegrationPointCount.Point3][]
                {
                    new double[3] {1.0 / 2.0, 1.0 / 2.0, 0.0},
                    new double[3] {0.0, 1.0 / 2.0, 1.0 / 2.0},
                    new double[3] {1.0 / 2.0, 0.0, 1.0 / 2.0}
                },
                Weights = new double[(int)TriangleIntegrationPointCount.Point3]
                {
                    1.0 / 3.0,
                    1.0 / 3.0,
                    1.0 / 3.0
                }
            },
            new IntegrationPoints{
                PointCount = (int)TriangleIntegrationPointCount.Point4,
                Ls = new double[(int)TriangleIntegrationPointCount.Point4][]
                {
                    new double[3] {TriangleIP3Alpha[0], TriangleIP3Alpha[0], TriangleIP3Alpha[0]},
                    new double[3] {TriangleIP3Alpha[1], TriangleIP3Alpha[2], TriangleIP3Alpha[2]},
                    new double[3] {TriangleIP3Alpha[2], TriangleIP3Alpha[1], TriangleIP3Alpha[2]},
                    new double[3] {TriangleIP3Alpha[2], TriangleIP3Alpha[2], TriangleIP3Alpha[1]}
                },
                Weights = new double[(int)TriangleIntegrationPointCount.Point4]
                {
                    -27.0 / 48.0,
                    25.0 / 48.0,
                    25.0 / 48.0,
                    25.0 / 48.0
                }
            },
            new IntegrationPoints{
                PointCount = (int)TriangleIntegrationPointCount.Point7,
                Ls = new double[(int)TriangleIntegrationPointCount.Point7][]
                {
                    new double[3] {TriangleIP7Alpha[0], TriangleIP7Alpha[0], TriangleIP7Alpha[0]},
                    new double[3] {TriangleIP7Alpha[1], TriangleIP7Alpha[2], TriangleIP7Alpha[2]},
                    new double[3] {TriangleIP7Alpha[2], TriangleIP7Alpha[1], TriangleIP7Alpha[2]},
                    new double[3] {TriangleIP7Alpha[2], TriangleIP7Alpha[2], TriangleIP7Alpha[1]},
                    new double[3] {TriangleIP7Alpha[3], TriangleIP7Alpha[4], TriangleIP7Alpha[4]},
                    new double[3] {TriangleIP7Alpha[4], TriangleIP7Alpha[3], TriangleIP7Alpha[4]},
                    new double[3] {TriangleIP7Alpha[4], TriangleIP7Alpha[4], TriangleIP7Alpha[3]}
                },
                Weights = new double[(int)TriangleIntegrationPointCount.Point7]
                {
                    0.225,
                    0.13239415,
                    0.13239415,
                    0.13239415,
                    0.12593918,
                    0.12593918,
                    0.12593918
                }
            }
        };

    }
}
