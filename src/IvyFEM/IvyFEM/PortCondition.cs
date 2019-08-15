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
        public IList<uint> EIds { get; protected set; } = null;
        public IList<uint> LoopIdsForPeriodic { get; protected set; } = null;
        public IList<uint> BcEIdsForPeriodic1 { get; protected set; } = null;
        public IList<uint> BcEIdsForPeriodic2 { get; protected set; } = null;
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
            IsPeriodic = false;
            EIds = eIds;
            LoopIdsForPeriodic = null;
            BcEIdsForPeriodic1 = null;
            BcEIdsForPeriodic2 = null;
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
            IsPeriodic = false;
            EIds = eIds;
            LoopIdsForPeriodic = null;
            BcEIdsForPeriodic1 = null;
            BcEIdsForPeriodic2 = null;
            ValueType = valueType;
            Dof = FieldValue.GetDof(ValueType);
            FixedDofIndexs = new List<uint>(fixedDofIndexs);
            IntAdditionalParameters = new List<int>();
            AdditionalParametersDof = additionalParametersDof;
        }

        // 周期構造導波路のポート条件
        public PortCondition(
            IList<uint> lIdsForPeriodic, IList<uint> eIdsForPeriodic1, IList<uint> eIdsForPeriodic2,
            FieldValueType valueType, IList<uint> fixedDofIndexs, uint additionalParametersDof)
        {
            IsPeriodic = true; // 周期構造
            EIds = null;
            LoopIdsForPeriodic = lIdsForPeriodic;
            BcEIdsForPeriodic1 = eIdsForPeriodic1;
            BcEIdsForPeriodic2 = eIdsForPeriodic2;
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
            EIds = null;
            if (src.EIds != null)
            {
                EIds = new List<uint>(src.EIds);
            }
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
