using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class MultipointConstraint
    {
        public IList<FieldFixedCad >FixedCads { get; private set; } = new List<FieldFixedCad>();
        public Constraint Constraint { get; set; } = null;

        public MultipointConstraint()
        {

        }

        public MultipointConstraint(
            IList<KeyValuePair<uint, CadElementType>> cadIdTypes, Constraint constraint)
        {
            Constraint = constraint;

            foreach (var pair in cadIdTypes)
            {
                uint cadId = pair.Key;
                CadElementType cadElemType = pair.Value;
                // Lagrangeの未定乗数は常に1自由度
                var fixedCad = new FieldFixedCad(cadId, cadElemType, FieldValueType.Scalar);
                FixedCads.Add(fixedCad);
            }
        }
    }
}
