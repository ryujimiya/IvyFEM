using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class ConstFieldFixedCad : FieldFixedCad
    {
        public double[] DoubleValues { get; protected set; } = null;
        public System.Numerics.Complex[] ComplexValues { get; protected set; } = null;

        public ConstFieldFixedCad(uint cadId, CadElementType cadElemType,
            FieldValueType valueType, IList<uint> fixedDofIndexs, double[] values) :
            base(cadId, cadElemType, valueType, fixedDofIndexs)
        {
            DoubleValues = new double[Dof];
            values.CopyTo(DoubleValues, 0);
            ComplexValues = null;
        }

        public ConstFieldFixedCad(uint cadId, CadElementType cadElemType,
            FieldValueType valueType, IList<uint> fixedDofIndexs, System.Numerics.Complex[] values) :
            base(cadId, cadElemType, valueType, fixedDofIndexs)
        {
            DoubleValues = null;
            ComplexValues = new System.Numerics.Complex[Dof];
            values.CopyTo(ComplexValues, 0);
        }

        public ConstFieldFixedCad(ConstFieldFixedCad src)
        {
            Copy(src);
        }

        public override void Copy(FieldFixedCad src)
        {
            base.Copy(src);

            ConstFieldFixedCad srcConst = src as ConstFieldFixedCad;
            DoubleValues = null;
            if (srcConst.DoubleValues != null)
            {
                DoubleValues = new double[srcConst.DoubleValues.Length];
                srcConst.DoubleValues.CopyTo(DoubleValues, 0);
            }
            ComplexValues = null;
            if (srcConst.ComplexValues != null)
            {
                ComplexValues = new System.Numerics.Complex[srcConst.ComplexValues.Length];
                srcConst.ComplexValues.CopyTo(ComplexValues, 0);
            }
        }

        public override double[] GetDoubleValues(int coId = -1)
        {
            return DoubleValues;
        }

        public override System.Numerics.Complex[] GetComplexValues(int coId = -1)
        {
            return ComplexValues;
        }
    }
}
