using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Lis
{
    unsafe public class LisSolver : IDisposable
    {
        internal NativeLisSolver* Native = null;

        public LisSolver()
        {
            int ret = IvyFEM.Lis.Functions.SolverCreate(this);
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

                int ret = IvyFEM.Lis.Functions.SolverDestroy(this);
                System.Diagnostics.Debug.Assert(ret == 0);

                disposedValue = true;
            }
        }

        ~LisSolver()
        {
            Dispose(false);
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }
        #endregion


        public int GetIter(out int iter)
        {
            return IvyFEM.Lis.Functions.SolverGetIter(this, out iter);
        }

        public int SetOption(string text)
        {
            return IvyFEM.Lis.Functions.SolverSetOption(text, this);
        }

        public int SetOptionC()
        {
            return IvyFEM.Lis.Functions.SolverSetOptionC(this);
        }

        public int Solve(LisMatrix A, LisVector b, LisVector x)
        {
            return IvyFEM.Lis.Functions.Solve(A, b, x, this);
        }
    }
}
