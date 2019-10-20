using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public abstract class EMWaveguide1DEigenBaseFEM : FEM
    {
        public uint QuantityId { get; protected set; } = 0;
        public uint PortId { get; protected set; } = 0;

        /// <summary>
        /// TEモードで実装した式をTMモードに流用するため
        ///   TEモードの場合は μ0
        ///   TMモードの場合は ε0
        /// </summary>
        public double ReplacedMu0 { get; set; } = Constants.Mu0;

        public IvyFEM.Lapack.DoubleMatrix Txx { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Ryy { get; protected set; } = null;
        public IvyFEM.Lapack.DoubleMatrix Uzz { get; protected set; } = null;

        // Solve
        // Input
        public double Frequency { get; set; }
        // Output
        public System.Numerics.Complex[] Betas { get; protected set; }
        public System.Numerics.Complex[][] EzEVecs { get; protected set; }

        public EMWaveguide1DEigenBaseFEM(FEWorld world, uint quantityId, uint portId)
        {
            World = world;
            QuantityId = quantityId;
            PortId = portId;
        }

        public abstract IvyFEM.Lapack.ComplexMatrix CalcBoundaryMatrix(
            double omega, System.Numerics.Complex[] betas, System.Numerics.Complex[][] ezEVecs);

        public abstract System.Numerics.Complex[] CalcIncidentVec(
            System.Numerics.Complex beta0, System.Numerics.Complex[] ezEVec0);

        public abstract System.Numerics.Complex[] CalcSMatrix(double omega, int incidentModeId,
            System.Numerics.Complex[] betas, System.Numerics.Complex[][] ezEVecs,
            System.Numerics.Complex[] Ez);

        public abstract System.Numerics.Complex CalcModeAmp(double omega, int modeIndex,
            System.Numerics.Complex[] betas, System.Numerics.Complex[][] ezEVecs,
            System.Numerics.Complex[] Ez);
    }
}
