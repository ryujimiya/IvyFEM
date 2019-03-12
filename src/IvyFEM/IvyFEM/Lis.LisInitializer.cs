using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM.Lis
{
    public class LisInitializer : IDisposable
    {
        public LisInitializer()
        {
            int ret = IvyFEM.Lis.Functions.Initialize();
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

                int ret = IvyFEM.Lis.Functions.Finalize();
                System.Diagnostics.Debug.Assert(ret == 0);

                disposedValue = true;
            }
        }

        ~LisInitializer()
        {
            Dispose(false);
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }
        #endregion
    }
}
