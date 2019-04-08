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

        public FieldFixedCad()
        {

        }

        // Boundary Condition
        public FieldFixedCad(uint cadId, CadElementType cadElemType,
            FieldValueType valueType, uint dof, IList<uint> fixedDofIndexs)
        {
            CadId = cadId;
            CadElemType = cadElemType;
            ValueType = valueType;
            Dof = dof;
            FixedDofIndexs = new List<uint>(fixedDofIndexs);
        }

        // Zero
        public FieldFixedCad(uint cadId, CadElementType cadElemType,
            FieldValueType valueType, uint dof)
        {
            CadId = cadId;
            CadElemType = cadElemType;
            ValueType = valueType;
            Dof = dof;
            FixedDofIndexs = new List<uint>();
            for (uint iDof = 0; iDof < dof; iDof++)
            {
                FixedDofIndexs.Add(iDof);
            }
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
