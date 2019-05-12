using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class FieldFixedCad
    {
        public uint CadId { get; set; } = 0;
        public CadElementType CadElemType { get; protected set; } = CadElementType.NotSet;
        public FieldValueType ValueType { get; protected set; } = FieldValueType.NoValue;
        public uint Dof { get; set; } = 0;
        public IList<uint> FixedDofIndexs { get; protected set; } = new List<uint>();
        public IList<int> IntAdditionalParameters { get; set; } = new List<int>();
        public IList<double> DoubleAdditionalParameters { get; set; } = new List<double>();
        public IList<System.Numerics.Complex> ComplexAdditionalParameters { get; set; } =
            new List<System.Numerics.Complex>();

        public FieldFixedCad()
        {

        }

        // Boundary Condition
        public FieldFixedCad(uint cadId, CadElementType cadElemType,
            FieldValueType valueType, IList<uint> fixedDofIndexs)
        {
            CadId = cadId;
            CadElemType = cadElemType;
            ValueType = valueType;
            Dof = FieldValue.GetDof(ValueType);
            FixedDofIndexs = new List<uint>(fixedDofIndexs);
            IntAdditionalParameters = new List<int>();
            DoubleAdditionalParameters = new List<double>();
            ComplexAdditionalParameters = new List<System.Numerics.Complex>();
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
            DoubleAdditionalParameters = new List<double>();
            ComplexAdditionalParameters = new List<System.Numerics.Complex>();
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
            DoubleAdditionalParameters = new List<double>(src.DoubleAdditionalParameters);
            ComplexAdditionalParameters = new List<System.Numerics.Complex>(src.ComplexAdditionalParameters);
        }

        public virtual double[] GetDoubleValues(int coId = -1)
        {
            return null;
        }

        public virtual System.Numerics.Complex[] GetComplexValues(int coId = -1)
        {
            return null;
        }
    }
}
