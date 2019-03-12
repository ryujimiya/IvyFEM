using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public interface IObject
    {
        void Copy(IObject src);
    }

    public class ObjectArray<T> where T : IObject, new()
    {
        private IList<uint> Index2Ids = new List<uint>();
        private IList<int> Id2Indexs = new List<int>();
        private IList<T> Objects = new List<T>();

        public ObjectArray()
        {

        }

        public ObjectArray(ObjectArray<T> src)
        {
            Copy(src);
        }

        public void Copy(ObjectArray<T> src)
        {
            Clear();
            foreach (uint id in src.Index2Ids)
            {
                Index2Ids.Add(id);
            }
            foreach (int index in src.Id2Indexs)
            {
                Id2Indexs.Add(index);
            }
            foreach (T obj in src.Objects)
            {
                T tmpObj = new T();
                tmpObj.Copy(obj);
                Objects.Add(tmpObj);
            }
        }

        public void Clear()
        {
            Index2Ids.Clear();
            Id2Indexs.Clear();
            Objects.Clear();
        }

        public T GetObject(uint id)
        {
            int index = GetArrayInd(id);
            System.Diagnostics.Debug.Assert(index >= 0 && index < Objects.Count);
            return Objects[index];
        }

        public uint AddObject(KeyValuePair<uint, T> idObj)
        {
            uint id1 = AddId(idObj.Key);
            System.Diagnostics.Debug.Assert(IsObjectId(id1));
            Objects.Add(idObj.Value);
            System.Diagnostics.Debug.Assert(Index2Ids.Count == Objects.Count);
            return id1;
        }

        public bool IsObjectId(uint id)
        {
            if (id == 0)
            {
                return false;
            }
            int index = GetArrayInd(id);
            if (index != -1)
            {
                return true;
            }
            return false;
        }

        public IList<uint> GetObjectIds()
        {
            return Index2Ids;
        }

        public uint GetFreeObjectId()
        {
            if (Id2Indexs.Count == 0)
            {
                return 1;
            }
            uint id;
            for (id = 1; id < Id2Indexs.Count; id++)
            {
                if (Id2Indexs[(int)id] == -1)
                {
                    return id;
                }
            }
            return (uint)Id2Indexs.Count;
        }

        public IList<uint> GetFreeObjectIds(uint size)
        {
            IList<uint> res = new List<uint>();
            for (int i = 0; i < size; i++)
            {
                res.Add(0);
            }
            uint isize = 0;
            System.Diagnostics.Debug.Assert(Id2Indexs.Count != 1);
            if (Id2Indexs.Count == 0)
            {
                for (;;)
                {
                    res[(int)isize] = isize + 1;
                    isize++;
                    if (isize == size)
                    {
                        return res;
                    }
                }
            }
            for (uint id = 1; id < Id2Indexs.Count; id++)
            {
                if (Id2Indexs[(int)id] == -1)
                {
                    res[(int)isize] = id;
                    isize++;
                    if (isize == size)
                    {
                        return res;
                    }
                }
            }
            for (int i = 0; ; i++)
            {
                res[(int)isize] = (uint)(Id2Indexs.Count + i);
                isize++;
                if (isize == size)
                {
                    return res;
                }
            }

            throw new InvalidOperationException();
            //return res;
        }

        public bool DeleteObject(uint id)
        {
            if (!IsObjectId(id))
            {
                return false;
            }
            int index = GetArrayInd(id);
            System.Diagnostics.Debug.Assert(index >= 0 && index < (int)Objects.Count);
            Objects.RemoveAt(index);
            System.Diagnostics.Debug.Assert(index >= 0 && index < (int)Index2Ids.Count);
            Index2Ids.RemoveAt(index);
            for (int i = 0; i < Id2Indexs.Count; i++)
            {
                Id2Indexs[i] = -1;
            }
            for (int i = 0; i < Index2Ids.Count; i++)
            {
                uint workId = Index2Ids[i];
                System.Diagnostics.Debug.Assert(workId < Id2Indexs.Count);
                Id2Indexs[(int)workId] = i;
            }
            return true;
        }

        public uint MaxId()
        {
            return (uint)Id2Indexs.Count;
        }

        private uint AddId(uint id)
        {
		    uint id1 = 0;
		    if ( !IsObjectId(id) && id > 0 && id < 255)
            {
			    id1 = id;
			    if (Id2Indexs.Count <= id1)
                {
                    int cnt = Id2Indexs.Count;
                    for (int i = cnt; i < id1 + 1; i++)
                    {
                        Id2Indexs.Add(-1);
                    }
                    Id2Indexs[(int)id1] = Id2Indexs.Count;
                }
			    Index2Ids.Add(id1);
			    Id2Indexs[(int)id1] = Index2Ids.Count - 1;
			}
            else
            {
                id1 = GetFreeObjectId();
                System.Diagnostics.Debug.Assert(!IsObjectId(id1));
                if (Id2Indexs.Count <= id1)
                {
                    int cnt = Id2Indexs.Count;
                    for (int i = cnt; i < id1 + 1; i++)
                    {
                        Id2Indexs.Add(-1);
                    }
                    Id2Indexs[(int)id1] = Index2Ids.Count;
                 }
                 Index2Ids.Add(id1);
                 Id2Indexs[(int)id1] = Index2Ids.Count - 1;

            }
            System.Diagnostics.Debug.Assert(IsObjectId(id1));
            return id1;
		}

        private int GetArrayInd(uint id)
        {
            int index = -1;
            if (Id2Indexs.Count <= id)
            {
                return index;
            }
            index = Id2Indexs[(int)id];
            return index;
        }

    }
}
