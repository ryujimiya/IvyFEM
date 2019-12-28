using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class PCWaveguidePortInfo
    {
        /// <summary>
        /// 固有値問題を反復で解く？
        ///   true: 取得するモード数が1または2の場合に固有値問題を反復で解く
        ///     ModeCountToSolve == 1の場合：線形固有方程式の反復を基本モードについて行う
        ///     ModeCountToSolve == 2の場合：線形固有方程式の反復を基本モードと高次モードについて行う
        ///   false: 全モードを一括で解く
        /// </summary>
        public bool IsSolveIter { get; set; } = false;
        /// <summary>
        /// 反復で解くモードの数
        /// </summary>
        public int ModeCountToSolveIter { get; set; } = 0;
        /// <summary>
        /// 緩慢変化包絡線近似？
        ///    true:  Φ = φ(x, y) exp(-jβx) の場合
        ///    false: Φ(x, y)を直接解く場合
        /// </summary>
        public bool IsSVEA { get; set; } = true;
        /// <summary>
        /// 格子定数
        /// </summary>
        public double LatticeA { get; set; } = 0.0;
        /// <summary>1
        /// 周期距離(X方向)
        /// </summary>
        public double PeriodicDistanceX { get; set; } = 0.0;
        /// <summary>1
        /// 周期距離(Y方向)
        /// for photonic band
        /// </summary>
        public double PeriodicDistanceY { get; set; } = 0.0;
        /// <summary>
        /// X方向に斜め
        /// </summary>
        public bool IsAslantX { get; set; } = false;
        /// <summary>
        /// 周期構造入出力導波路 ループID
        /// </summary>
        public IList<uint> LoopIds { get; set; } = new List<uint>();
        /// <summary>
        /// 周期構造入出力導波路 境界辺ID
        /// </summary>
        public IList<uint> BcEdgeIds1 { get; set; } = new List<uint>();
        /// <summary>
        /// 周期構造入出力導波路 境界辺ID
        /// </summary>
        public IList<uint> BcEdgeIds2 { get; set; } = new List<uint>();
        /// <summary>
        /// 周期構造上下導波路 境界辺ID
        /// for photonic band
        /// </summary>
        public IList<uint> BcEdgeIds3 { get; set; } = new List<uint>();
        /// <summary>
        /// 周期構造上下導波路 境界辺ID
        /// for photonic band
        /// </summary>
        public IList<uint> BcEdgeIds4 { get; set; } = new List<uint>();
        /// <summary>
        /// 周期構造の方向
        /// false: X方向
        /// true: Y方向
        /// </summary>
        public bool IsYDirectionPeriodic { get; set; } = false;
        /// <summary>
        /// 境界2の節点順番が境界1の順番の逆？
        /// </summary>
        public bool IsPortBc2Reverse { get; set; } = false;
        /// <summary>
        /// 境界4の節点順番が境界1の順番の逆？
        /// for photonic band
        /// </summary>
        public bool IsPortBc4Reverse { get; set; } = false;
        /// <summary>
        /// 境界の節点
        /// </summary>
        public IList<IList<int>> BcNodess { get; set; } = null;
        /// <summary>
        /// 周期構造入出力導波路 フォトニック導波路のチャンネル(欠陥部)節点座標IDリストのリスト
        /// </summary>
        public IList<IList<int>> PCChannelCoIds { get; set; } = new List<IList<int>>();
        /// <summary>
        /// 最小屈折率
        /// </summary>
        public double MinEffN { get; set; } = 0.0;
        /// <summary>
        /// 最大屈折率
        /// </summary>
        public double MaxEffN { get; set; } = 1.0;
        /// <summary>
        /// 考慮する波数ベクトルの最小値
        /// </summary>
        public double MinWaveNum { get; set; } = 0.0;
        /// <summary>
        /// 考慮する波数ベクトルの最大値
        /// </summary>
        public double MaxWaveNum { get; set; } = 0.5;
        /// <summary>
        /// 考慮する最小周波数(βを与えてk0を解く場合)
        /// </summary>
        public double MinFrequency { get; set; } = 0.0;
        /// <summary>
        /// 考慮する最大周波数(βを与えてk0を解く場合)
        /// </summary>
        public double MaxFrequency { get; set; } = double.MaxValue;

        /// <summary>
        /// モード追跡
        /// </summary>
        public System.Numerics.Complex[][] PrevModeEVecs { get; set; } = null;

        public PCWaveguidePortInfo()
        {

        }

        public void SetupAfterMakeElements(FEWorld world, uint quantityId, uint portId)
        {
            System.Diagnostics.Debug.Assert(BcEdgeIds1.Count > 0);
            System.Diagnostics.Debug.Assert(BcEdgeIds1.Count == BcEdgeIds2.Count);
            System.Diagnostics.Debug.Assert(BcEdgeIds3.Count == BcEdgeIds4.Count);
            int bcCnt = BcEdgeIds3.Count == 0 ? 2 : 4; // 4 for photonic band

            // 境界上の節点
            BcNodess = new List<IList<int>>();
            for (uint bcIndex = 0; bcIndex < bcCnt; bcIndex++)
            {
                IList<int> coIds = world.GetPeriodicPortBcCoIds(quantityId, portId, bcIndex);
                IList<int> bcNodes = new List<int>();
                BcNodess.Add(bcNodes);
                foreach (int coId in coIds)
                {
                    int nodeId = world.PortCoord2Node(quantityId, portId, coId);
                    if (nodeId == -1)
                    {
                        continue;
                    }
                    bcNodes.Add(nodeId);
                }
            }
        }
    }
}
