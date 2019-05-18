using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class DistributedPortCondition : PortCondition
    {
        public IList<int> CoIds { get; protected set; } = new List<int>();
        public IList<double[]> DoubleValuess { get; protected set; } = new List<double[]>();
        public IList<System.Numerics.Complex[]> ComplexValuess { get; protected set; } =
            new List<System.Numerics.Complex[]>();
        public IList<double[]> DoubleAdditionalParameterss { get; protected set; } = new List<double[]>();
        public IList<System.Numerics.Complex[]> ComplexAdditionalParameterss { get; protected set; } =
            new List<System.Numerics.Complex[]>();

        public DistributedPortCondition(IList<uint> eIds,
            FieldValueType valueType, IList<uint> fixedDofIndexs,
            uint additionalParametersDof = 0) :
            base(eIds, valueType, fixedDofIndexs, additionalParametersDof)
        {

        }

        public DistributedPortCondition(DistributedPortCondition src)
        {
            Copy(src);
        }

        public override void Copy(PortCondition src)
        {
            base.Copy(src);

            DistributedPortCondition srcDist = src as DistributedPortCondition;
            CoIds = new List<int>(srcDist.CoIds);
            DoubleValuess = new List<double[]>();
            foreach (double[] srcValues in srcDist.DoubleValuess)
            {
                double[] values = new double[srcValues.Length];
                srcValues.CopyTo(values, 0);
                DoubleValuess.Add(values);
            }
            ComplexValuess = new List<System.Numerics.Complex[]>();
            foreach (System.Numerics.Complex[] srcValues in srcDist.ComplexValuess)
            {
                System.Numerics.Complex[] values = new System.Numerics.Complex[srcValues.Length];
                srcValues.CopyTo(values, 0);
                ComplexValuess.Add(values);
            }
        }

        public void InitCoordIds(IList<int> coIds)
        {
            CoIds = new List<int>(coIds);
        }

        public void InitDoubleValues()
        {
            DoubleValuess = new List<double[]>();
            for (int i = 0; i < CoIds.Count; i++)
            {
                double[] values = new double[Dof];
                DoubleValuess.Add(values);
            }
        }

        public void InitComplexValues()
        {
            for (int i = 0; i < CoIds.Count; i++)
            {
                System.Numerics.Complex[] values = new System.Numerics.Complex[Dof];
                ComplexValuess.Add(values);
            }
        }

        public void InitDoubleAdditionalParameters()
        {
            if (AdditionalParametersDof != 0)
            {
                for (int i = 0; i < CoIds.Count; i++)
                {
                    double[] param = new double[AdditionalParametersDof];
                    DoubleAdditionalParameterss.Add(param);
                }
            }
        }

        public void InitComplexAdditionalParameters()
        {
            if (AdditionalParametersDof != 0)
            {
                for (int i = 0; i < CoIds.Count; i++)
                {
                    System.Numerics.Complex[] param = new System.Numerics.Complex[AdditionalParametersDof];
                    ComplexAdditionalParameterss.Add(param);
                }
            }
        }

        public override double[] GetDoubleValues(int coId = -1)
        {
            double[] values = null;
            int index = CoIds.IndexOf(coId);
            if (index != -1)
            {
                values = DoubleValuess[index];
            }
            return values;
        }

        public override System.Numerics.Complex[] GetComplexValues(int coId = -1)
        {
            System.Numerics.Complex[] values = null;
            int index = CoIds.IndexOf(coId);
            if (index != -1)
            {
                values = ComplexValuess[index];
            }
            return values;
        }

        public override double[] GetDoubleAdditionalParameters(int coId = -1)
        {
            double[] param = null;
            int index = CoIds.IndexOf(coId);
            if (index != -1)
            {
                param = DoubleAdditionalParameterss[index];
            }
            return param;
        }

        public override System.Numerics.Complex[] GetComplexAdditionalParameters(int coId = -1)
        {
            System.Numerics.Complex[] param = null;
            int index = CoIds.IndexOf(coId);
            if (index != -1)
            {
                param = ComplexAdditionalParameterss[index];
            }
            return param;
        }
    }
}
