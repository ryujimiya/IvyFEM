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
        ///     PropModeCntToSolve == 1の場合：線形固有方程式の反復を基本モードについて行う
        ///     PropModeCntToSolve == 2の場合：線形固有方程式の反復を基本モードと高次モードについて行う
        ///   false: 全モードを一括で解く
        /// </summary>
        public bool IsSolveEigenItr = true;
        /// <summary>
        /// 緩慢変化包絡線近似？
        ///    true:  Φ = φ(x, y) exp(-jβx) の場合
        ///    false: Φ(x, y)を直接解く場合
        /// </summary>
        public bool IsSVEA = true;
        /// <summary>
        /// モード追跡する？
        /// </summary>
        public bool IsModeTrace = true;

        /// <summary>
        /// 格子定数
        /// </summary>
        public double LatticeA = 0.0;
        /// <summary>1
        /// 周期距離
        /// </summary>
        public double PeriodicDistance = 0.0;
        /// <summary>
        /// 周期構造入出力導波路 ループID
        /// </summary>
        public IList<uint> LoopIds = new List<uint>();
        /// <summary>
        /// 周期構造入出力導波路 境界辺ID
        /// </summary>
        public IList<uint> BcEdgeIds1 = new List<uint>();
        /// <summary>
        /// 周期構造入出力導波路 境界辺ID
        /// </summary>
        public IList<uint> BcEdgeIds2 = new List<uint>();
        /// <summary>
        /// 周期構造の方向
        /// false: X方向
        /// true: Y方向
        /// </summary>
        public bool IsYDirectionPeriodic = false;
        /// <summary>
        /// 境界2の節点順番が境界1の順番の逆？
        /// </summary>
        public bool IsPortBc2Reverse = false;
        /// <summary>
        /// 境界の節点
        /// </summary>
        public IList<IList<int>> BcNodess = null;
        /// <summary>
        /// 周期構造入出力導波路 フォトニック導波路のチャンネル(欠陥部)節点座標IDリストのリスト
        /// </summary>
        public IList<IList<int>> PCChannelCoIds = new List<IList<int>>();
        /// <summary>
        /// 最小屈折率
        /// </summary>
        public double MinEffN = 0.0;
        /// <summary>
        /// 最大屈折率
        /// </summary>
        public double MaxEffN = 1.0;
        /// <summary>
        /// 考慮する波数ベクトルの最小値
        /// </summary>
        public double MinWaveNum = 0.0;
        /// <summary>
        /// 考慮する波数ベクトルの最大値
        /// </summary>
        public double MaxWaveNum = 0.5;
        /// <summary>
        /// TEモードで実装した式をTMモードに流用するため
        ///   TEモードの場合は μ0
        ///   TMモードの場合は ε0
        /// </summary>
        public double ReplacedMu0 = Constants.Mu0;
        
        /// <summary>
        /// モード追跡
        /// </summary>
        public System.Numerics.Complex[][] PrevModeEVecs = null;

        public PCWaveguidePortInfo()
        {

        }

        public void SetupAfterMakeElements(FEWorld world, uint quantityId, uint portId)
        {
            // 境界上の節点
            BcNodess = new List<IList<int>>();
            for (uint bcIndex = 0; bcIndex < 2; bcIndex++)
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
            int bcNodeCnt = BcNodess[0].Count;
            System.Diagnostics.Debug.Assert(BcNodess[0].Count == BcNodess[1].Count);
            // Y方向に周期構造?
            bool isYDirectionPeriodic = false;
            {
                IList<int> bcNodes = BcNodess[0];
                int nodeId1 = bcNodes[0];
                int coId1 = world.PortNode2Coord(quantityId, portId, nodeId1);
                double[] coord1 = world.GetCoord(quantityId, coId1);
                int nodeId2 = bcNodes[bcNodeCnt - 1];
                int coId2 = world.PortNode2Coord(quantityId, portId, nodeId2);
                double[] coord2 = world.GetCoord(quantityId, coId2);
                if (Math.Abs(coord1[1] - coord2[1]) < 1.0e-12)
                {
                    isYDirectionPeriodic = true;
                }
            }
            IsYDirectionPeriodic = isYDirectionPeriodic;
        }
    }
}
