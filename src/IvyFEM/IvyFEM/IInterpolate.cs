using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public interface IInterpolate
    {
        uint GetNodeCount();
        double[] GetNodeL(int nodeId);
        double[] CalcN(double[] L);
        double[][] CalcNu(double[] L);
        double[,][] CalcNuv(double[] L);
        double[] CalcSN();
        double[,] CalcSNN();
        double[,][,] CalcSNuNv();
        double[][,] CalcSNuN();
    }
}
