using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Constants
    {
        public const double C0 = 2.99792458e+8;
        public const double Mu0 = 4.0e-7 * Math.PI;
        public const double Ep0 = 8.85418782e-12;//1.0 / (Mu0 * C0 * C0);

        public const double PrecisionLowerLimit = 1.0e-16;
    }

    public enum RotMode
    {
        RotModeNotSet,
        RotMode2D,
        RotMode2DH,
        RotMode3D
    }

    public enum CurveType
    {
        CurveEndPoint,
        CurveLine,
        CurveArc,
        CurvePolyline,
        CurveBezier
    }

    public enum CadElementType
    {
        NotSet,
        Vertex,
        Edge,
        Loop,
        Solid
    }

    public enum MeshType
    {
        NotSet,
        Vertex,
        Bar,
        Tri,
        Tet
    }

    public enum ElementType
    {
        NotSet,
        Point,
        Line,
        Tri,
        Tet
    }

    public enum FieldValueType
    {
        NoValue,
        Scalar,
        Vector2,
        Vector3,
        SymmetricTensor2,
        ZScalar,
        ZVector2,
        ZVector3
    }

    public enum FieldValueNodeType
    {
        Node,
        Bubble,
        ElementNode,
        ElementEdge
    }

    [Flags]
    public enum FieldDerivativeType
    {
        Value = 1,
        Velocity = 2,
        Acceleration = 4
    }

    public enum FieldShowType
    {
        Real,
        Abs,
        ZReal,
        ZImaginary,
        ZAbs
    }

    public enum VectorFieldDrawerType
    {
        NotSet,
        Vector,
        SymmetricTensor2
    }

    public enum LineIntegrationPointCount
    {
        Point1 = 1,
        Point2 = 2,
        Point3 = 3,
        Point4 = 4,
        Point5 = 5,
        Point10 = 10
    }

    public enum TriangleIntegrationPointCount
    {
        Point1 = 1,
        Point3 = 3,
        Point4 = 4,
        Point7 = 7,
        Point25 = 25
    }

    public enum TetrahedronIntegrationPointCount
    {
        Point1 = 1,
        Point4 = 4,
        Point5 = 5
    }

    public enum FiniteElementType
    {
        ScalarLagrange,
        ScalarConstant,
        ScalarBell,
        ScalarHermite,
        Edge
    }

    public enum EqualityType
    {
        Eq,
        LessEq,
        GreaterEq
    }

    public enum GaussianType
    {
        Normal, // 素のガウシアンパルス
        SinModulation // 正弦波変調
    }

    public enum ElasticBCType
    {
        ExternalForce
    }

    public enum FluidEquationType
    {
        Stokes,
        StdGNavierStokes,
        SUPGNavierStokes,
        StdGVorticity,
        SUPGVorticity,
        StdGPressurePoisson,
        StdGPressurePoissonWithBell,
        SUPGPressurePoissonWithBell
    }

    public enum FlowVorticityBCType
    {
        TangentialFlow,
        Outflow
    }

    public enum FlowPressureBCType
    {
        NoConstraint,
        NormalInflow,
        Outflow
    }

    public enum EMWaveguideType
    {
        HPlane2D,
        EPlane2D
    }
}
