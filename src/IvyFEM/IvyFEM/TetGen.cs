using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class TetGen
    {
        // 拡張子を除いたファイルパス
        public string FileNameBase { get; set; }
        // Input
        public double ELen { get; set; }
        public IList<OpenTK.Vector3d> Vecs0 { get; set; }
        public IList<MeshTri3D> Tris0 { get; set; }
        public IList<OpenTK.Vector3d> Holes0 { get; set; }
        public IList<OpenTK.Vector3d> InsideVecs0 { get; set; }
        // Output
        public bool Success { get; private set; }
        public IList<OpenTK.Vector3d> Vecs { get; private set; }
        public IList<MeshTet> Tets { get; private set; }

        public TetGen()
        {

        }

        public void Execute()
        {
            Success = true;
            bool ret;

            DeleteFiles(); // 前のファイルを削除

            SaveToSmeshFile();
            SaveToAdditionalNodeFile();

            ExecuteProcess();

            ret = LoadFromNodeFile();
            if (!ret)
            {
                Success = false;
            }
            ret = LoadFromEleFile();
            if (!ret)
            {
                Success = false;
            }
            
            DeleteFiles(); // 後片付け
        }

        private void SaveToSmeshFile()
        {
            string fileName = FileNameBase + ".poly";

            using (StreamWriter sw = new StreamWriter(fileName))
            {
                // node list
                sw.WriteLine("{0} {1} {2} {3}", Vecs0.Count, 3, 0, 0);
                for (int iVec = 0; iVec < Vecs0.Count; iVec++)
                {
                    OpenTK.Vector3d vec = Vecs0[iVec];
                    uint nodeId = (uint)(iVec + 1);
                    sw.WriteLine("{0} {1} {2} {3}", nodeId, vec.X, vec.Y, vec.Z);
                }

                // facet list
                sw.WriteLine("{0} {1}", Tris0.Count, 0);
                for (int iTri = 0; iTri < Tris0.Count; iTri++)
                {
                    MeshTri3D tri = Tris0[iTri];
                    uint nodeId1 = tri.V[0] + 1;
                    uint nodeId2 = tri.V[1] + 1;
                    uint nodeId3 = tri.V[2] + 1;
                    sw.WriteLine("{0}", 1);
                    sw.WriteLine("{0} {1} {2} {3}", 3, nodeId1, nodeId2, nodeId3);
                }

                // hole list
                sw.WriteLine("{0}", Holes0.Count);
                for (int iHole = 0; iHole < Holes0.Count; iHole++)
                {
                    OpenTK.Vector3d hole = Holes0[iHole];
                    uint holeId = (uint)(iHole + 1);
                    sw.WriteLine("{0} {1} {2} {3}", holeId, hole.X, hole.Y, hole.Z);
                }

                // region list
                sw.WriteLine("{0}", 0);
            }
        }

        private void SaveToAdditionalNodeFile()
        {
            string fileName = FileNameBase + ".a.node";

            using (StreamWriter sw = new StreamWriter(fileName))
            {
                // node list
                sw.WriteLine("{0} {1} {2} {3}", InsideVecs0.Count, 3, 0, 0);
                for (int iVec = 0; iVec < InsideVecs0.Count; iVec++)
                {
                    OpenTK.Vector3d vec = InsideVecs0[iVec];
                    uint nodeId = (uint)(iVec + 1);
                    sw.WriteLine("{0} {1} {2} {3}", nodeId, vec.X, vec.Y, vec.Z);
                }
            }
        }

        private void ExecuteProcess()
        {
            string fileName = FileNameBase + ".poly";
            string dir = Path.GetDirectoryName(FileNameBase);
            string processFileName = dir + "\\tetgen.exe";
            string arg = string.Format("-pYqi {0}", fileName);

            if (!File.Exists(processFileName))
            {
                System.Diagnostics.Debug.Assert(false);
                return;
            }
            System.Diagnostics.Process process = System.Diagnostics.Process.Start(processFileName, arg);
            process.WaitForExit();
            return;
        }

        private bool LoadFromNodeFile()
        {
            Vecs = new List<OpenTK.Vector3d>();

            string fileName = FileNameBase + ".1" + ".node";
            if (!File.Exists(fileName))
            {
                System.Diagnostics.Debug.Assert(false);
                return false;
            }
            using (StreamReader sr = new StreamReader(fileName))
            {
                string[] values;

                values = GetValues(sr);
                System.Diagnostics.Debug.Assert(values.Length == 4);
                int nodeCnt = int.Parse(values[0]);
                int dim = int.Parse(values[1]);
                int attrCnt = int.Parse(values[2]);
                int boundaryMaker = int.Parse(values[3]);
                System.Diagnostics.Debug.Assert(dim == 3);
                System.Diagnostics.Debug.Assert(attrCnt == 0);
                System.Diagnostics.Debug.Assert(boundaryMaker == 0);

                for (int iVec = 0; iVec < nodeCnt; iVec++)
                {
                    values = GetValues(sr);
                    System.Diagnostics.Debug.Assert(values.Length == (1 + 3));
                    int nodeId = int.Parse(values[0]) - 1;
                    System.Diagnostics.Debug.Assert(nodeId == iVec);
                    double x = double.Parse(values[1]);
                    double y = double.Parse(values[2]);
                    double z = double.Parse(values[3]);

                    Vecs.Add(new OpenTK.Vector3d(x, y, z));
                }
            }
            return true;
        }

        private bool LoadFromEleFile()
        {
            Tets = new List<MeshTet>();

            string fileName = FileNameBase + ".1" + ".ele";

            if (!File.Exists(fileName))
            {
                System.Diagnostics.Debug.Assert(false);
                return false;
            }
            using (StreamReader sr = new StreamReader(fileName))
            {
                string[] values;

                values = GetValues(sr);
                System.Diagnostics.Debug.Assert(values.Length == 3);
                int tetCnt = int.Parse(values[0]);
                int tetNodeNo = int.Parse(values[1]);
                int regionAttr = int.Parse(values[2]);
                System.Diagnostics.Debug.Assert(tetNodeNo == 4);
                System.Diagnostics.Debug.Assert(regionAttr == 0 || regionAttr == 1);

                for (int iTet = 0; iTet < tetCnt; iTet++)
                {
                    values = GetValues(sr);
                    System.Diagnostics.Debug.Assert(values.Length >= (1 + tetNodeNo));
                    int attrValue = 1; // interior
                    if (values.Length == (1 + tetNodeNo + 1))
                    {
                        attrValue = int.Parse(values[1 + tetNodeNo]);
                    }
                    if (attrValue == -1)
                    {
                        // exterior
                        continue;
                    }
                    MeshTet tet = new MeshTet();
                    Tets.Add(tet);
                    int tetId = int.Parse(values[0]) - 1;
                    System.Diagnostics.Debug.Assert(tetId == iTet);
                    for (int iNode = 0; iNode < tetNodeNo; iNode++)
                    {
                        int nodeId = int.Parse(values[iNode + 1]) - 1;
                        System.Diagnostics.Debug.Assert(nodeId >= 0);
                        tet.V[iNode] = (uint)nodeId;
                    }
                }
            }
            return true;
        }

        private void DeleteFiles()
        {
            string[] fileNames = {
                FileNameBase + ".poly",
                FileNameBase + ".a.node",
                FileNameBase + ".1" + ".edge",
                FileNameBase + ".1" + ".ele",
                FileNameBase + ".1" + ".face",
                FileNameBase + ".1" + ".node"
            };

            foreach (string fileName in fileNames)
            {
                if (File.Exists(fileName))
                {
                    File.Delete(fileName);
                }
            }
        }

        private string[] GetValues(StreamReader sr)
        {
            string line;
            string[] values = new string[0];
            const char delimiter = ' ';

            line = sr.ReadLine();
            if (line == null || line == "")
            {
                return values;
            }

            line = Regex.Replace(line, @"#.*$", "");
            line = line.Trim();
            line = Regex.Replace(line, @"\s+", " ");
            if (line == null || line == "")
            {
                return values;
            }

            values = line.Split(delimiter);
            return values;
        }
    }
}
