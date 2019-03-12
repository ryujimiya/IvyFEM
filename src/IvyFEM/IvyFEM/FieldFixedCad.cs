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
        public CadElementType CadElemType { get; private set; } = CadElementType.NotSet;
        public FieldValueType ValueType { get; private set; } = FieldValueType.NoValue;
        public uint Dof { get; set; } = 0;
        public IList<uint> FixedDofIndexs { get; private set; } = new List<uint>();
        public double[] DoubleValues { get; set; } = null;
        public System.Numerics.Complex[] ComplexValues { get; set; } = null;

        public FieldFixedCad()
        {

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
            DoubleValues = null;
            ComplexValues = null;
        }

        public FieldFixedCad(uint cadId, CadElementType cadElemType,
            FieldValueType valueType, uint dof,
            IList<uint> fixedDofIndexs, double[] values)
        {
            CadId = cadId;
            CadElemType = cadElemType;
            ValueType = valueType;
            Dof = dof;
            FixedDofIndexs = new List<uint>(fixedDofIndexs);
            DoubleValues = new double[Dof];
            values.CopyTo(DoubleValues, 0);
            ComplexValues = null;
        }

        public FieldFixedCad(uint cadId, CadElementType cadElemType,
            FieldValueType valueType, uint dof,
            IList<uint> fixedDofIndexs, System.Numerics.Complex[] values)
        {
            CadId = cadId;
            CadElemType = cadElemType;
            ValueType = valueType;
            Dof = dof;
            FixedDofIndexs = new List<uint>(fixedDofIndexs);
            DoubleValues = null;
            ComplexValues = new System.Numerics.Complex[dof];
            values.CopyTo(ComplexValues, 0);
        }

        public FieldFixedCad(FieldFixedCad src)
        {
            CadId = src.CadId;
            CadElemType = src.CadElemType;
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
        }
    }
}
