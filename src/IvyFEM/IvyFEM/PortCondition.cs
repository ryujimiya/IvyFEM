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
        public uint Dof { get; protected set; } = 0;
        public IList<uint> FixedDofIndexs { get; protected set; } = new List<uint>();
        public IList<int> IntAdditionalParameters { get; set; } = new List<int>();
        public uint AdditionalParametersDof { get; protected set; } = 0;

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
            IntAdditionalParameters = new List<int>();
            AdditionalParametersDof = 0;
        }

        public PortCondition(IList<uint> eIds,
            FieldValueType valueType, IList<uint> fixedDofIndexs, uint additionalParametersDof)
        {
            EIds = eIds;
            ValueType = valueType;
            Dof = FieldValue.GetDof(ValueType);
            FixedDofIndexs = new List<uint>(fixedDofIndexs);
            IntAdditionalParameters = new List<int>();
            AdditionalParametersDof = additionalParametersDof;
        }

        public PortCondition(PortCondition src)
        {
            Copy(src);
        }

        public virtual void Copy(PortCondition src)
        {
            EIds = new List<uint>(src.EIds);
            ValueType = src.ValueType;
            Dof = src.Dof;
            FixedDofIndexs = new List<uint>(src.FixedDofIndexs);
            IntAdditionalParameters = new List<int>(src.IntAdditionalParameters);
            AdditionalParametersDof = src.AdditionalParametersDof;
        }

        public virtual double[] GetDoubleValues(int coId = -1)
        {
            return null;
        }

        public virtual System.Numerics.Complex[] GetComplexValues(int coId = -1)
        {
            return null;
        }

        public virtual double[] GetDoubleAdditionalParameters(int coId = -1)
        {
            return null;
        }

        public virtual System.Numerics.Complex[] GetComplexAdditionalParameters(int coId = -1)
        {
            return null;
        }

    }
}
