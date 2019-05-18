using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class FieldFixedCad
    {
        public uint CadId { get; protected set; } = 0;
        public CadElementType CadElemType { get; protected set; } = CadElementType.NotSet;
        public FieldValueType ValueType { get; protected set; } = FieldValueType.NoValue;
        public uint Dof { get; protected set; } = 0;
        public IList<uint> FixedDofIndexs { get; protected set; } = new List<uint>();
        public IList<int> IntAdditionalParameters { get; set; } = new List<int>();
        public uint AdditionalParametersDof { get; protected set; } = 0;

        public FieldFixedCad()
        {

        }

        // Boundary Condition
        public FieldFixedCad(uint cadId, CadElementType cadElemType,
            FieldValueType valueType, IList<uint> fixedDofIndexs, uint additionalParametersDof = 0)
        {
            CadId = cadId;
            CadElemType = cadElemType;
            ValueType = valueType;
            Dof = FieldValue.GetDof(ValueType);
            FixedDofIndexs = new List<uint>(fixedDofIndexs);
            IntAdditionalParameters = new List<int>();
            AdditionalParametersDof = additionalParametersDof;
        }

        // Zero
        public FieldFixedCad(uint cadId, CadElementType cadElemType,
            FieldValueType valueType)
        {
            CadId = cadId;
            CadElemType = cadElemType;
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

        public FieldFixedCad(FieldFixedCad src)
        {
            Copy(src);
        }

        public virtual void Copy(FieldFixedCad src)
        {
            CadId = src.CadId;
            CadElemType = src.CadElemType;
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
