using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using LisScalar = System.Numerics.Complex;

namespace IvyFEM.Lis
{
    unsafe public class LisVector : IDisposable
    {
        internal NativeLisVector* Native = null;

        public LisVector(int comm = IvyFEM.Lis.Constants.LisCommWorld)
        {
            int ret = IvyFEM.Lis.Functions.VectorCreate(comm, this);
            System.Diagnostics.Debug.Assert(ret == 0);
        }

        #region IDisposable Support
        private bool disposedValue = false; // 重複する呼び出しを検出するには

        protected virtual void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                if (disposing)
                {

                }

                int ret = IvyFEM.Lis.Functions.VectorDestroy(this);
                System.Diagnostics.Debug.Assert(ret == 0);

                disposedValue = true;
            }
        }

        ~LisVector()
        {
            Dispose(false);
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }
        #endregion

        public LisVector(double[] values)
        {
            int n = values.Length;
            int ret;
            ret = SetSize(0, n);
            System.Diagnostics.Debug.Assert(ret == 0);
            for (int i = 0; i < n; i++)
            {
                ret = SetValue(SetValueFlag.LisInsValue, i, values[i]);
                System.Diagnostics.Debug.Assert(ret == 0);
            }
        }

        public LisVector(System.Numerics.Complex[] values)
        {
            int n = values.Length;
            int ret;
            ret = SetSize(0, n);
            System.Diagnostics.Debug.Assert(ret == 0);
            for (int i = 0; i < n; i++)
            {
                ret = SetValue(SetValueFlag.LisInsValue, i, values[i]);
                System.Diagnostics.Debug.Assert(ret == 0);
            }
        }

        public int SetSize(int localN, int globalN)
        {
            return IvyFEM.Lis.Functions.VectorSetSize(this, localN, globalN);
        }

        public int GetSize(out int localN, out int globalN)
        {
            return IvyFEM.Lis.Functions.VectorGetSize(this, out localN, out globalN);
        }

        public int GetRange(out int @is, out int ie)
        {
            return IvyFEM.Lis.Functions.VectorGetRange(this, out @is, out ie);
        }

        public int GetValue(int i, out LisScalar value)
        {
            return IvyFEM.Lis.Functions.VectorGetValue(this, i, out value);
        }

        public int GetValues(int start, int count, LisScalar[] value)
        {
            return IvyFEM.Lis.Functions.VectorGetValues(this, start, count, value);
        }

        public int SetValue(SetValueFlag flag, int i, LisScalar value)
        {
            return IvyFEM.Lis.Functions.VectorSetValue(flag, i, value, this);
        }

        public int SetValues(SetValueFlag flag, int count, int[] index, LisScalar[] value)
        {
            return IvyFEM.Lis.Functions.VectorSetValues(flag, count, index, value, this);
        }

        public int SetValues2(SetValueFlag flag, int start, int count, LisScalar[] value)
        {
            return IvyFEM.Lis.Functions.VectorSetValues2(flag, start, count, value, this);
        }

        public int SetAll(LisScalar alpha)
        {
            return IvyFEM.Lis.Functions.VectorSetAll(alpha, this);
        }

        public int Conjugate()
        {
            return IvyFEM.Lis.Functions.VectorConjugate(this);
        }
    }
}
