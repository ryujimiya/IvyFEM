using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class PortCondition
    {
        public IList<uint> EIds { get; protected set; } = null;
        public FieldValueType ValueType { get; protected set; } = FieldValueType.NoValue;
        public uint Dof { get; set; } = 0;
        public IList<uint> FixedDofIndexs { get; protected set; } = new List<uint>();
        public double[] DoubleValues { get; protected set; } = null;
        public System.Numerics.Complex[] ComplexValues { get; protected set; } = null;
        public IList<int> IntAdditionalParameters { get; set; } = new List<int>();
        public IList<double> DoubleAdditionalParameters { get; set; } = new List<double>();
        public IList<System.Numerics.Complex> ComplexAdditionalParameters { get; set; } =
            new List<System.Numerics.Complex>();

        public PortCondition()
        {

        }

        public PortCondition(IList<uint> eIds, FieldValueType valueType)
        {
            EIds = eIds;
            ValueType = valueType;
            Dof = FieldValue.GetDof(ValueType);
            FixedDofIndexs = new List<uint>();
            for (uint iDof = 0; iDof < Dof; iDof++)
            {
                FixedDofIndexs.Add(iDof);
            }
            DoubleValues = null;
            ComplexValues = null;
            IntAdditionalParameters = new List<int>();
            DoubleAdditionalParameters = new List<double>();
            ComplexAdditionalParameters = new List<System.Numerics.Complex>();
        }

        public PortCondition(IList<uint> eIds,
            FieldValueType valueType, IList<uint> fixedDofIndexs, IList<double> fixedValues)
        {
            EIds = eIds;
            ValueType = valueType;
            Dof = FieldValue.GetDof(ValueType);
            FixedDofIndexs = new List<uint>(fixedDofIndexs);

            System.Diagnostics.Debug.Assert(FixedDofIndexs.Count == fixedValues.Count);
            DoubleValues = new double[Dof];
            for (int i = 0; i < FixedDofIndexs.Count; i++)
            {
                uint iDof = FixedDofIndexs[i];
                double value = fixedValues[i];
                DoubleValues[iDof] = value;
            }
            ComplexValues = null;
            IntAdditionalParameters = new List<int>();
            DoubleAdditionalParameters = new List<double>();
            ComplexAdditionalParameters = new List<System.Numerics.Complex>();
        }

        public PortCondition(IList<uint> eIds,
            FieldValueType valueType, IList<uint> fixedDofIndexs, IList<System.Numerics.Complex> fixedValues)
        {
            EIds = eIds;
            ValueType = valueType;
            Dof = FieldValue.GetDof(ValueType);
            FixedDofIndexs = new List<uint>(fixedDofIndexs);

            System.Diagnostics.Debug.Assert(FixedDofIndexs.Count == fixedValues.Count);
            DoubleValues = null;
            ComplexValues = new System.Numerics.Complex[Dof];
            for (int i = 0; i < FixedDofIndexs.Count; i++)
            {
                uint iDof = FixedDofIndexs[i];
                System.Numerics.Complex value = fixedValues[i];
                ComplexValues[iDof] = value;
            }
            IntAdditionalParameters = new List<int>();
            DoubleAdditionalParameters = new List<double>();
            ComplexAdditionalParameters = new List<System.Numerics.Complex>();
        }

        public PortCondition(PortCondition src)
        {
            Copy(src);
        }

        public void Copy(PortCondition src)
        {
            EIds = new List<uint>(src.EIds);
            ValueType = src.ValueType;
            Dof = src.Dof;
            FixedDofIndexs = new List<uint>(src.FixedDofIndexs);
            DoubleValues = null;
            if (src.DoubleValues != null)
            {
                DoubleValues = new double[src.DoubleValues.Length];
                src.DoubleValues.CopyTo(DoubleValues, 0);
            }
            ComplexValues = null;
            if (src.ComplexValues != null)
            {
                ComplexValues = new System.Numerics.Complex[src.ComplexValues.Length];
                src.ComplexValues.CopyTo(ComplexValues, 0);
            }
            IntAdditionalParameters = new List<int>(src.IntAdditionalParameters);
            DoubleAdditionalParameters = new List<double>(src.DoubleAdditionalParameters);
            ComplexAdditionalParameters = new List<System.Numerics.Complex>(src.ComplexAdditionalParameters);
        }

        public double[] GetDoubleValues()
        {
            return DoubleValues;
        }

        public System.Numerics.Complex[] GetComplexValues()
        {
            return ComplexValues;
        }
    }
}
