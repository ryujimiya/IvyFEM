using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Hyperelastic2DTDFEM : Hyperelastic2DDerivedBaseFEM
    {
        public double TimeStep { get; private set; } = 0;
        public double NewmarkBeta { get; private set; } = 1.0 / 4.0;
        public double NewmarkGamma { get; private set; } = 1.0 / 2.0;
        public uint UValueId { get; set; } = 0;
        public uint PrevUValueId { get; set; } = 0;
        public uint LValueId { get; set; } = 0;

        public Hyperelastic2DTDFEM(FEWorld world,
            double timeStep,
            double newmarkBeta, double newmarkGamma,
            uint uValueId, uint prevUValueId,
            uint lValueId)
        {
            World = world;
            TimeStep = timeStep;
            NewmarkBeta = newmarkBeta;
            NewmarkGamma = newmarkGamma;
            UValueId = uValueId;
            PrevUValueId = prevUValueId;
            LValueId = lValueId;
            SetupCalcABs();
        }

        protected void SetupCalcABs()
        {
            CalcElementABs.Clear();
            CalcElementABs.Add(CalcMooneyRivlinHyperelasticElementAB);
            CalcElementABs.Add(CalcOgdenRivlinHyperelasticElementAB);
            //CalcElementABs.Add(CalcOgdenOriginalRivlinIncompressibleHyperelasticElementAB);
        }

        public void UpdateFieldValues()
        {
            double dt = TimeStep;
            double beta = NewmarkBeta;
            double gamma = NewmarkGamma;

            var uFV = World.GetFieldValue(UValueId);
            var prevUFV = World.GetFieldValue(PrevUValueId);
            prevUFV.Copy(uFV);

            World.UpdateFieldValueValuesFromNodeValues(UValueId, FieldDerivativeType.Value, U);
            World.UpdateFieldValueValuesFromNodeValues(LValueId, FieldDerivativeType.Value, U);

            double[] u = uFV.GetDoubleValues(FieldDerivativeType.Value);
            double[] velU = uFV.GetDoubleValues(FieldDerivativeType.Velocity);
            double[] accU = uFV.GetDoubleValues(FieldDerivativeType.Acceleration);
            double[] prevU = prevUFV.GetDoubleValues(FieldDerivativeType.Value);
            double[] prevVelU = prevUFV.GetDoubleValues(FieldDerivativeType.Velocity);
            double[] prevAccU = prevUFV.GetDoubleValues(FieldDerivativeType.Acceleration);

            uint uCoCnt = uFV.GetPointCount();
            uint uQuantityId = uFV.QuantityId;
            int uDof = (int)uFV.Dof;
            System.Diagnostics.Debug.Assert(uCoCnt == World.GetCoordCount(uQuantityId));
            System.Diagnostics.Debug.Assert(uDof == 2);
            System.Diagnostics.Debug.Assert(u.Length == uCoCnt * uDof);
            for (int iPt = 0; iPt < uCoCnt; iPt++)
            {
                for (int iDof = 0; iDof < uDof; iDof++)
                {
                    int index = iPt * uDof + iDof;
                    velU[index] =
                        (gamma / (beta * dt)) * (u[index] - prevU[index]) +
                        (1.0 - gamma / beta) * prevVelU[index] +
                        dt * (1.0 - gamma / (2.0 * beta)) * prevAccU[index];
                    accU[index] =
                        (1.0 / (beta * dt * dt)) * (u[index] - prevU[index]) -
                        (1.0 / (beta * dt)) * prevVelU[index] -
                        (1.0 / (2.0 * beta) - 1.0) * prevAccU[index];
                }
            }
        }
    }
}
