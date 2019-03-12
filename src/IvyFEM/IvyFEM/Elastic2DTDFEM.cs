using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public partial class Elastic2DTDFEM : Elastic2DDerivedBaseFEM
    {
        public double TimeStep { get; private set; } = 0;
        public double NewmarkBeta { get; private set; } = 1.0 / 4.0;
        public double NewmarkGamma { get; private set; } = 1.0 / 2.0;
        uint ValueId { get; set; } = 0;
        uint PrevValueId { get; set; } = 0;

        public Elastic2DTDFEM(FEWorld world, double timeStep,
            double newmarkBeta, double newmarkGamma,
            uint valueId, uint prevValueId)
        {
            World = world;
            TimeStep = timeStep;
            NewmarkBeta = newmarkBeta;
            NewmarkGamma = newmarkGamma;
            ValueId = valueId;
            PrevValueId = prevValueId;
            SetupCalcABs();
        }

        protected void SetupCalcABs()
        {
            CalcElementABs.Clear();
            CalcElementABs.Add(CalcLinearElasticElementAB);
            CalcElementABs.Add(CalcSaintVenantHyperelasticElementAB);
        }

        public void UpdateFieldValues()
        {
            double dt = TimeStep;
            double beta = NewmarkBeta;
            double gamma = NewmarkGamma;

            var FV = World.GetFieldValue(ValueId);
            var prevFV = World.GetFieldValue(PrevValueId);
            prevFV.Copy(FV);

            World.UpdateFieldValueValuesFromNodeValues(ValueId, FieldDerivativeType.Value, U);

            double[] u = FV.GetDoubleValues(FieldDerivativeType.Value);
            double[] vel = FV.GetDoubleValues(FieldDerivativeType.Velocity);
            double[] acc = FV.GetDoubleValues(FieldDerivativeType.Acceleration);
            double[] prevU = prevFV.GetDoubleValues(FieldDerivativeType.Value);
            double[] prevVel = prevFV.GetDoubleValues(FieldDerivativeType.Velocity);
            double[] prevAcc = prevFV.GetDoubleValues(FieldDerivativeType.Acceleration);

            uint coCnt = FV.GetPointCount();
            uint quantityId = FV.QuantityId;
            int dof = (int)FV.Dof;
            System.Diagnostics.Debug.Assert(coCnt == World.GetCoordCount(quantityId));
            System.Diagnostics.Debug.Assert(dof == 2);
            System.Diagnostics.Debug.Assert(u.Length == coCnt * dof);
            for (int iPt = 0; iPt < coCnt; iPt++)
            {
                for(int iDof = 0; iDof < dof; iDof++)
                {
                    int index = iPt * dof + iDof;
                    vel[index] =
                        (gamma / (beta * dt)) * (u[index] - prevU[index]) +
                        (1.0 - gamma / beta) * prevVel[index] +
                        dt * (1.0 - gamma / (2.0 * beta)) * prevAcc[index];
                    acc[index] =
                        (1.0 / (beta * dt * dt)) * (u[index] - prevU[index]) -
                        (1.0 / (beta * dt)) * prevVel[index] -
                        (1.0 / (2.0 * beta) - 1.0) * prevAcc[index];
                }
            }
        }
    }
}
