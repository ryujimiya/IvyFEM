using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class PortCondition
    {
        public bool IsPeriodic { get; protected set; } = false;
        public IList<uint> CadIds { get; protected set; } = null;
        public CadElementType CadElemType { get; protected set; } = CadElementType.Edge;
        public IList<uint> EIds => CadElemType == CadElementType.Edge ? CadIds : null;
        public IList<uint> LIds => CadElemType == CadElementType.Loop ? CadIds : null;
        // 2D
        public IList<uint> LoopIdsForPeriodic { get; protected set; } = null;
        public IList<uint> BcEIdsForPeriodic1 { get; protected set; } = null;
        public IList<uint> BcEIdsForPeriodic2 { get; protected set; } = null;
        public IList<uint> BcEIdsForPeriodic3 { get; protected set; } = null; // for photonic band
        public IList<uint> BcEIdsForPeriodic4 { get; protected set; } = null; // for photonic band

        public FieldValueType ValueType { get; protected set; } = FieldValueType.NoValue;
        public uint Dof { get; protected set; } = 0;
        public IList<uint> FixedDofIndexs { get; protected set; } = new List<uint>();
        public IList<int> IntAdditionalParameters { get; set; } = new List<int>();
        public uint AdditionalParametersDof { get; protected set; } = 0;

        public PortCondition()
        {

        }

        public PortCondition(
            IList<uint> cadIds, CadElementType cadElemType, FieldValueType valueType)
        {
            System.Diagnostics.Debug.Assert(
                cadElemType == CadElementType.Edge || cadElemType == CadElementType.Loop);

            IsPeriodic = false;
            CadIds = cadIds;
            CadElemType = cadElemType;
            LoopIdsForPeriodic = null;
            BcEIdsForPeriodic1 = null;
            BcEIdsForPeriodic2 = null;
            BcEIdsForPeriodic3 = null;
            BcEIdsForPeriodic4 = null;
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

        public PortCondition(IList<uint> cadIds, CadElementType cadElemType,
            FieldValueType valueType, IList<uint> fixedDofIndexs, uint additionalParametersDof)
        {
            System.Diagnostics.Debug.Assert(
                cadElemType == CadElementType.Edge || cadElemType == CadElementType.Loop);

            IsPeriodic = false;
            CadIds = cadIds;
            CadElemType = cadElemType;
            LoopIdsForPeriodic = null;
            BcEIdsForPeriodic1 = null;
            BcEIdsForPeriodic2 = null;
            BcEIdsForPeriodic3 = null;
            BcEIdsForPeriodic4 = null;
            ValueType = valueType;
            Dof = FieldValue.GetDof(ValueType);
            FixedDofIndexs = new List<uint>(fixedDofIndexs);
            IntAdditionalParameters = new List<int>();
            AdditionalParametersDof = additionalParametersDof;
        }

        // 周期構造導波路のポート条件
        public PortCondition(
            CadElementType cadElemType,
            IList<uint> lIdsForPeriodic, IList<uint> eIdsForPeriodic1, IList<uint> eIdsForPeriodic2,
            FieldValueType valueType, IList<uint> fixedDofIndexs, uint additionalParametersDof)
        {
            System.Diagnostics.Debug.Assert(cadElemType == CadElementType.Edge);

            IsPeriodic = true; // 周期構造
            CadIds = null;
            CadElemType = cadElemType;
            LoopIdsForPeriodic = lIdsForPeriodic;
            BcEIdsForPeriodic1 = eIdsForPeriodic1;
            BcEIdsForPeriodic2 = eIdsForPeriodic2;
            BcEIdsForPeriodic3 = null;
            BcEIdsForPeriodic4 = null;
            ValueType = valueType;
            Dof = FieldValue.GetDof(ValueType);
            FixedDofIndexs = new List<uint>(fixedDofIndexs);
            IntAdditionalParameters = new List<int>();
            AdditionalParametersDof = additionalParametersDof;
        }

        // 周期構造導波路のポート条件(photonic band)
        public PortCondition(
            CadElementType cadElemType,
            IList<uint> lIdsForPeriodic, IList<uint> eIdsForPeriodic1, IList<uint> eIdsForPeriodic2,
            IList<uint> eIdsForPeriodic3, IList<uint> eIdsForPeriodic4,
            FieldValueType valueType, IList<uint> fixedDofIndexs, uint additionalParametersDof)
        {
            System.Diagnostics.Debug.Assert(cadElemType == CadElementType.Edge);

            IsPeriodic = true; // 周期構造
            CadIds = null;
            CadElemType = cadElemType;
            LoopIdsForPeriodic = lIdsForPeriodic;
            BcEIdsForPeriodic1 = eIdsForPeriodic1;
            BcEIdsForPeriodic2 = eIdsForPeriodic2;
            BcEIdsForPeriodic3 = eIdsForPeriodic3;
            BcEIdsForPeriodic4 = eIdsForPeriodic4;
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
            IsPeriodic = src.IsPeriodic;
            CadIds = null;
            if (src.CadIds != null)
            {
                CadIds = new List<uint>(src.CadIds);
            }
            CadElemType = src.CadElemType;
            LoopIdsForPeriodic = null;
            if (src.LoopIdsForPeriodic != null)
            {
                LoopIdsForPeriodic = new List<uint>(src.LoopIdsForPeriodic);
            }
            BcEIdsForPeriodic1 = null;
            if (src.BcEIdsForPeriodic1 != null)
            {
                BcEIdsForPeriodic1 = new List<uint>(src.BcEIdsForPeriodic1);
            }
            BcEIdsForPeriodic2 = null;
            if (src.BcEIdsForPeriodic2 != null)
            {
                BcEIdsForPeriodic2 = new List<uint>(src.BcEIdsForPeriodic2);
            }
            BcEIdsForPeriodic3 = null;
            if (src.BcEIdsForPeriodic3 != null)
            {
                BcEIdsForPeriodic3 = new List<uint>(src.BcEIdsForPeriodic3);
            }
            BcEIdsForPeriodic4 = null;
            if (src.BcEIdsForPeriodic4 != null)
            {
                BcEIdsForPeriodic4 = new List<uint>(src.BcEIdsForPeriodic4);
            }
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
