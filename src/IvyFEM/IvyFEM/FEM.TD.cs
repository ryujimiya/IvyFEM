using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract partial class FEM
    {
        protected void UpdateFieldValuesNewmarkBetaTimeDomain(
            double[] U,
            uint valueId, uint prevValueId,
            double timeStep,
            double newmarkBeta, double newmarkGamma)
        {
            double dt = timeStep;
            double beta = newmarkBeta;
            double gamma = newmarkGamma;

            var FV = World.GetFieldValue(valueId);
            var prevFV = World.GetFieldValue(prevValueId);
            prevFV.Copy(FV);

            World.UpdateFieldValueValuesFromNodeValues(valueId, FieldDerivativeType.Value, U);

            double[] u = FV.GetDoubleValues(FieldDerivativeType.Value);
            double[] vel = FV.GetDoubleValues(FieldDerivativeType.Velocity);
            double[] acc = FV.GetDoubleValues(FieldDerivativeType.Acceleration);
            double[] prevU = prevFV.GetDoubleValues(FieldDerivativeType.Value);
            double[] prevVel = prevFV.GetDoubleValues(FieldDerivativeType.Velocity);
            double[] prevAcc = prevFV.GetDoubleValues(FieldDerivativeType.Acceleration);

            uint pointCnt = FV.GetPointCount();
            uint quantityId = FV.QuantityId;
            int dof = (int)FV.Dof;
            if (FV.NodeType == FieldValueNodeType.ElementEdge)
            {
                // ElementEdge
                System.Diagnostics.Debug.Assert(pointCnt == World.GetEdgeCount(quantityId));
            }
            else
            {
                // ElementNode
                System.Diagnostics.Debug.Assert(pointCnt == World.GetCoordCount(quantityId));
            }
            System.Diagnostics.Debug.Assert(u.Length == pointCnt * dof);
            for (int iPt = 0; iPt < pointCnt; iPt++)
            {
                for (int iDof = 0; iDof < dof; iDof++)
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
