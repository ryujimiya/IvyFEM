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
        public double[] DoubleAdditionalParameters { get; protected set; } = null;
        public System.Numerics.Complex[] ComplexAdditionalParameters { get; protected set; } = null;

        public ConstFieldFixedCad(uint cadId, CadElementType cadElemType,
            FieldValueType valueType, IList<uint> fixedDofIndexs, IList<double> fixedValues,
            uint additionalParametersDof = 0) :
            base(cadId, cadElemType, valueType, fixedDofIndexs, additionalParametersDof)
        {
            System.Diagnostics.Debug.Assert(FixedDofIndexs.Count == fixedValues.Count);
            DoubleValues = new double[Dof];
            for (int i = 0; i < FixedDofIndexs.Count; i++)
            {
                uint iDof = FixedDofIndexs[i];
                double value = fixedValues[i];
                DoubleValues[iDof] = value;
            }
            ComplexValues = null;
            DoubleAdditionalParameters = new double[AdditionalParametersDof];
            ComplexAdditionalParameters = null;
        }

        public ConstFieldFixedCad(uint cadId, CadElementType cadElemType,
            FieldValueType valueType, IList<uint> fixedDofIndexs, IList<System.Numerics.Complex> fixedValues,
            uint additionalParametersDof = 0) :
            base(cadId, cadElemType, valueType, fixedDofIndexs, additionalParametersDof)
        {
            System.Diagnostics.Debug.Assert(fixedDofIndexs.Count == fixedValues.Count);
            DoubleValues = null;
            ComplexValues = new System.Numerics.Complex[Dof];
            for (int i = 0; i < FixedDofIndexs.Count; i++)
            {
                uint iDof = FixedDofIndexs[i];
                System.Numerics.Complex value = fixedValues[i];
                ComplexValues[iDof] = value;
            }
            DoubleAdditionalParameters = null;
            ComplexAdditionalParameters = new System.Numerics.Complex[additionalParametersDof];
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
            DoubleAdditionalParameters = null;
            if (srcConst.DoubleAdditionalParameters != null)
            {
                DoubleAdditionalParameters = new double[srcConst.DoubleAdditionalParameters.Length];
                srcConst.DoubleAdditionalParameters.CopyTo(DoubleAdditionalParameters, 0);
            }
            ComplexAdditionalParameters = null;
            if (srcConst.ComplexAdditionalParameters != null)
            {
                ComplexAdditionalParameters = new System.Numerics.Complex[srcConst.ComplexAdditionalParameters.Length];
                srcConst.ComplexAdditionalParameters.CopyTo(ComplexAdditionalParameters, 0);
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

        public override double[] GetDoubleAdditionalParameters(int coId = -1)
        {
            return DoubleAdditionalParameters;
        }

        public override System.Numerics.Complex[] GetComplexAdditionalParameters(int coId = -1)
        {
            return ComplexAdditionalParameters;
        }
    }
}
