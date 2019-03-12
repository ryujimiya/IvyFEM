using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class Serializer : IDisposable
    {
        private string FileName;
        private StreamReader Reader;
        private StreamWriter Writer;
        private long Pos;
        public uint Depth { get; private set; }
        public bool IsLoading { get; private set; }

        public Serializer(string fileName, bool isLoading)
        {
            FileName = fileName;
            IsLoading = isLoading;
            Depth = 1;
            ////////////////
            Pos = 0;

            Reader = null;
            Writer = null;

            Open();
        }

        #region IDisposable Support
        private bool disposedValue = false; // 重複する呼び出しを検出するには

        protected virtual void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                if (disposing)
                {
                    // TODO: マネージ状態を破棄します (マネージ オブジェクト)。
                    Close();
                }

                // TODO: アンマネージ リソース (アンマネージ オブジェクト) を解放し、下のファイナライザーをオーバーライドします。
                // TODO: 大きなフィールドを null に設定します。

                disposedValue = true;
            }
        }

        // TODO: 上の Dispose(bool disposing) にアンマネージ リソースを解放するコードが含まれる場合にのみ、ファイナライザーをオーバーライドします。
        // ~Serializer() {
        //   // このコードを変更しないでください。クリーンアップ コードを上の Dispose(bool disposing) に記述します。
        //   Dispose(false);
        // }

        // このコードは、破棄可能なパターンを正しく実装できるように追加されました。
        public void Dispose()
        {
            // このコードを変更しないでください。クリーンアップ コードを上の Dispose(bool disposing) に記述します。
            Dispose(true);
            // TODO: 上のファイナライザーがオーバーライドされる場合は、次の行のコメントを解除してください。
            // GC.SuppressFinalize(this);
        }
        #endregion

        private void Open()
        {
            if (IsLoading)
            {
                Reader = new StreamReader(FileName);
            }
            else
            {
                Writer = new StreamWriter(FileName);
            }
        }

        private void Close()
        {
            if (IsLoading)
            {
                Reader.Dispose();
            }
            else
            {
                Writer.Dispose();
            }
        }

        public string ReadDepthClassName()
        {
            System.Diagnostics.Debug.Assert(IsLoading);
            string line;
            string workStr;

            line = Reader.ReadLine();
            char[] buff = line.ToCharArray();
            System.Diagnostics.Debug.Assert(buff[0] == '#');
            workStr = line.Substring(1);
            uint idepth0 = uint.Parse(workStr);

            line = Reader.ReadLine();
            string className = line;
            System.Diagnostics.Debug.Assert(Depth == idepth0);

            return className;
        }

        public void WriteDepthClassName(string className)
        {
            System.Diagnostics.Debug.Assert(!IsLoading);
            string line;

            line = string.Format("#{0}", Depth);
            Writer.WriteLine(line);

            line = className;
            Writer.WriteLine(line);
        }

        public void ShiftDepth(bool isAdd)
        {
            if (isAdd)
            {
                Depth++;
                return;
            }
            System.Diagnostics.Debug.Assert(Depth > 1);
            Depth -= 1;
        }

        public string[] GetValues()
        {
            System.Diagnostics.Debug.Assert(IsLoading);
            string line;

            line = Reader.ReadLine();

            string[] values = line.Split(
                new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

            return values;
        }

        public string ReadLine()
        {
            System.Diagnostics.Debug.Assert(IsLoading);
            string line;

            line = Reader.ReadLine();

            return line;
        }

        public void WriteLine(string line)
        {
            System.Diagnostics.Debug.Assert(!IsLoading);

            Writer.WriteLine(line);
        }
    }
}
