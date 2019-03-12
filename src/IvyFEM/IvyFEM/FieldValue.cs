using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class FieldValue : IObject
    {
        public FieldValueType Type { get; set; } = FieldValueType.NoValue;
        public FieldDerivativeType DerivativeType { get; set; } = 0;
        public uint QuantityId { get; set; } = 0;
        public uint Dof { get; set; } = 1;
        public bool IsBubble { get; set; } = false;
        public FieldShowType ShowType { get; set; } = FieldShowType.Real;
        public double[] DoubleValues { get; set; } = null;
        public double[] DoubleVelocityValues { get; set; } = null;
        public double[] DoubleAccelerationValues { get; set; } = null;
        public System.Numerics.Complex[] ComplexValues { get; set; } = null;
        public System.Numerics.Complex[] ComplexVelocityValues { get; set; } = null;
        public System.Numerics.Complex[] ComplexAccelerationValues { get; set; } = null;

        public FieldValue()
        {

        }

        public FieldValue(FieldValue src)
        {
            Copy(src);
        }

        public void Copy(IObject src)
        {
            FieldValue srcFV = src as FieldValue;
            Type = srcFV.Type;
            DerivativeType = srcFV.DerivativeType;
            QuantityId = srcFV.QuantityId;
            Dof = srcFV.Dof;
            IsBubble = srcFV.IsBubble;
            ShowType = srcFV.ShowType;
            CopyValues(srcFV);
        }

        public void CopyValues(FieldValue src)
        {
            DoubleValues = null;
            if (src.DoubleValues != null)
            {
                DoubleValues = new double[src.DoubleValues.Length];
                src.DoubleValues.CopyTo(DoubleValues, 0);
            }
            DoubleVelocityValues = null;
            if (src.DoubleVelocityValues != null)
            {
                DoubleVelocityValues = new double[src.DoubleVelocityValues.Length];
                src.DoubleVelocityValues.CopyTo(DoubleVelocityValues, 0);
            }
            DoubleAccelerationValues = null;
            if (src.DoubleAccelerationValues != null)
            {
                DoubleAccelerationValues = new double[src.DoubleAccelerationValues.Length];
                src.DoubleAccelerationValues.CopyTo(DoubleAccelerationValues, 0);
            }

            ComplexValues = null;
            if (src.ComplexValues != null)
            {
                ComplexValues = new System.Numerics.Complex[src.ComplexValues.Length];
                src.ComplexValues.CopyTo(ComplexValues, 0);
            }
            ComplexVelocityValues = null;
            if (src.ComplexVelocityValues != null)
            {
                ComplexVelocityValues = new System.Numerics.Complex[src.ComplexVelocityValues.Length];
                src.ComplexVelocityValues.CopyTo(ComplexVelocityValues, 0);
            }
            ComplexAccelerationValues = null;
            if (src.ComplexAccelerationValues != null)
            {
                ComplexAccelerationValues = new System.Numerics.Complex[src.ComplexAccelerationValues.Length];
                src.ComplexAccelerationValues.CopyTo(ComplexAccelerationValues, 0);
            }
        }

        public void AllocValues(uint dof, uint pointCnt)
        {
            Dof = dof;
            if (Type == FieldValueType.ZScalar)
            {
                // complex
                if (DerivativeType.HasFlag(FieldDerivativeType.Value))
                {
                    ComplexValues = new System.Numerics.Complex[pointCnt * Dof];
                }
                if (DerivativeType.HasFlag(FieldDerivativeType.Velocity))
                {
                    ComplexVelocityValues = new System.Numerics.Complex[pointCnt * Dof];
                }
                if (DerivativeType.HasFlag(FieldDerivativeType.Acceleration))
                {
                    ComplexAccelerationValues = new System.Numerics.Complex[pointCnt * Dof];
                }
            }
            else
            {
                // double
                if (DerivativeType.HasFlag(FieldDerivativeType.Value))
                {
                    DoubleValues = new double[pointCnt * Dof];
                }
                if (DerivativeType.HasFlag(FieldDerivativeType.Velocity))
                {
                    DoubleVelocityValues = new double[pointCnt * Dof];
                }
                if (DerivativeType.HasFlag(FieldDerivativeType.Acceleration))
                {
                    DoubleAccelerationValues = new double[pointCnt * Dof];
                }
            }
        }

        public uint GetPointCount()
        {
            if (Type == FieldValueType.ZScalar)
            {
                if (ComplexValues == null)
                {
                    return 0;
                }
                return (uint)(ComplexValues.Length / Dof);
            }
            else
            {
                if (DoubleValues == null)
                {
                    return 0;
                }
                return (uint)(DoubleValues.Length / Dof);
            }
            throw new InvalidOperationException();
        }

        public double[] GetDoubleValues(FieldDerivativeType dt)
        {
            double[] values = null;
            if (dt.HasFlag(FieldDerivativeType.Value) && DoubleValues != null)
            {
                values = DoubleValues;
            }
            else if (dt.HasFlag(FieldDerivativeType.Velocity) && DoubleVelocityValues != null)
            {
                values = DoubleVelocityValues;
            }
            else if (dt.HasFlag(FieldDerivativeType.Acceleration) && DoubleAccelerationValues != null)
            {
                values = DoubleAccelerationValues;
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return values;
        }

        public System.Numerics.Complex[] GetComplexValues(FieldDerivativeType dt)
        {
            System.Numerics.Complex[] values = null;
            if (dt.HasFlag(FieldDerivativeType.Value) && ComplexValues != null)
            {
                values = ComplexValues;
            }
            else if (dt.HasFlag(FieldDerivativeType.Velocity) && ComplexVelocityValues != null)
            {
                values = ComplexVelocityValues;
            }
            else if (dt.HasFlag(FieldDerivativeType.Acceleration) && ComplexAccelerationValues != null)
            {
                values = ComplexAccelerationValues;
            }
            else
            {
                System.Diagnostics.Debug.Assert(false);
            }
            return values;
        }

        public double[] GetDoubleValue(int coId, FieldDerivativeType dt)
        {
            double[] values = GetDoubleValues(dt);
            double[] value = new double[Dof];
            for (int iDof = 0; iDof < Dof; iDof++)
            {
                value[iDof] = values[coId * Dof + iDof];
            }
            return value;
        }

        public System.Numerics.Complex[] GetComplexValue(int coId, FieldDerivativeType dt)
        {
            System.Numerics.Complex[] values = GetComplexValues(dt);
            System.Numerics.Complex[] value = new System.Numerics.Complex[Dof];
            for (int iDof = 0; iDof < Dof; iDof++)
            {
                value[iDof] = values[coId * Dof + iDof];
            }
            return value;
        }

        public double GetShowValue(int coId, int iDof, FieldDerivativeType dt)
        {
            double value = 0;
            switch(ShowType)
            {
                case FieldShowType.Real:
                    {
                        double[] values = GetDoubleValues(dt);
                        value = values[coId * Dof + iDof];
                    }
                    break;
                case FieldShowType.Abs:
                    {
                        double[] values = GetDoubleValues(dt);
                        value = Math.Abs(values[coId * Dof + iDof]);
                    }
                    break;
                case FieldShowType.ZReal:
                    {
                        System.Numerics.Complex[] values = GetComplexValues(dt);
                        value = values[coId * Dof].Real;
                    }
                    break;
                case FieldShowType.ZImaginary:
                    {
                        System.Numerics.Complex[] values = GetComplexValues(dt);
                        value = values[coId * Dof].Imaginary;
                    }
                    break;
                case FieldShowType.ZAbs:
                    {
                        System.Numerics.Complex[] values = GetComplexValues(dt);
                        value = values[coId * Dof].Magnitude;
                    }
                    break;
            }
            return value;
        }

        public void GetMinMaxShowValue(out double min, out double max, int iDof, FieldDerivativeType dt)
        {
            min = Double.MaxValue;
            max = Double.MinValue;

            uint ptCnt = GetPointCount();
            for (int coId = 0; coId < ptCnt; coId++)
            {
                double value = GetShowValue(coId, iDof, dt);
                if (value < min)
                {
                    min = value;
                }
                if (value > max)
                {
                    max = value;
                }
            }
        }

    }
}
