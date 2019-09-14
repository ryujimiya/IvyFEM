using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IvyFEM
{
    public class PCWaveguideUtils
    {
        /// <summary>
        /// 境界を分割する
        /// </summary>
        /// <param name="cad2D"></param>
        /// <param name="eId"></param>
        /// <param name="divCnt"></param>
        /// <param name="x1"></param>
        /// <param name="y1"></param>
        /// <param name="x2"></param>
        /// <param name="y2"></param>
        public static void DivideBoundary(
            CadObject2D cad2D, uint eId, int divCnt, double x1, double y1, double x2, double y2)
        {
            double signedWidthX = x2 - x1;
            double signedWidthY = y2 - y1;
            for (int i = divCnt - 1; i >= 1; i--)
            {
                double x = x1 + i * signedWidthX / divCnt;
                double y = y1 + i * signedWidthY / divCnt;
                AddVertexRes resAddVertex = cad2D.AddVertex(CadElementType.Edge, eId, new OpenTK.Vector2d(x, y));
                uint addVId = resAddVertex.AddVId;
                uint addEId = resAddVertex.AddEId;
                System.Diagnostics.Debug.Assert(addVId != 0);
            }
        }

        /// <summary>
        /// ロッドを追加する
        /// </summary>
        /// <param name="cad2D"></param>
        /// <param name="baseLoopId"></param>
        /// <param name="x0"></param>
        /// <param name="y0"></param>
        /// <param name="rodRadius"></param>
        /// <param name="rodCircleDiv"></param>
        /// <param name="rodRadiusDiv"></param>
        /// <returns></returns>
        public static uint AddRod(
            CadObject2D cad2D,
            uint baseLoopId, double x0, double y0, double rodRadius, int rodCircleDiv, int rodRadiusDiv)
        {
            IList<OpenTK.Vector2d> pts = new List<OpenTK.Vector2d>();
            {
                // メッシュ形状を整えるためにロッドの中心に頂点を追加
                uint centerVId = cad2D.AddVertex(CadElementType.Loop, baseLoopId, new OpenTK.Vector2d(x0, y0)).AddVId;
                System.Diagnostics.Debug.Assert(centerVId != 0);
            }
            // ロッドの分割数調整
            for (int k = 1; k < rodRadiusDiv; k++)
            {
                for (int itheta = 0; itheta < rodCircleDiv; itheta++)
                {
                    double theta = itheta * 2.0 * Math.PI / rodCircleDiv;
                    double x = x0 + (k * rodRadius / rodRadiusDiv) * Math.Cos(theta);
                    double y = y0 + (k * rodRadius / rodRadiusDiv) * Math.Sin(theta);
                    uint addVId = cad2D.AddVertex(CadElementType.Loop, baseLoopId, new OpenTK.Vector2d(x, y)).AddVId;
                    System.Diagnostics.Debug.Assert(addVId != 0);
                }
            }
            // ロッド
            for (int itheta = 0; itheta < rodCircleDiv; itheta++)
            {
                double theta = itheta * 2.0 * Math.PI / rodCircleDiv;
                double x = x0 + rodRadius * Math.Cos(theta);
                double y = y0 + rodRadius * Math.Sin(theta);
                pts.Add(new OpenTK.Vector2d(x, y));
            }
            uint lId = cad2D.AddPolygon(pts, baseLoopId).AddLId;
            System.Diagnostics.Debug.Assert(lId != 0);
            return lId;
        }

        /// <summary>
        /// ロッド(半分)を追加する
        /// </summary>
        /// <param name="cad2D"></param>
        /// <param name="baseLoopId"></param>
        /// <param name="x0"></param>
        /// <param name="y0"></param>
        /// <param name="rodRadius"></param>
        /// <param name="rodCircleDiv"></param>
        /// <param name="rodRadiusDiv"></param>
        /// <returns></returns>
        public static uint AddHalfRod(
            CadObject2D cad2D,
            uint baseLoopId, double x0, double y0, double rodRadius, int rodCircleDiv, int rodRadiusDiv,
            double startAngle, double additionalAngle = 0.0,
            bool isReverseAddVertex = false,
            uint stVId = uint.MaxValue, uint edVId = uint.MaxValue)
        {
            System.Diagnostics.Debug.Assert(additionalAngle < 360.0 / rodCircleDiv);
            IList<OpenTK.Vector2d> pts = new List<OpenTK.Vector2d>();
            System.Diagnostics.Debug.Assert(
                (startAngle == 0.0) || (startAngle == 90.0) || (startAngle == 180.0) || (startAngle == 270.0));
            if (Math.Abs(additionalAngle) >= IvyFEM.Constants.PrecisionLowerLimit)
            {
                // メッシュ形状を整えるためにロッドの中心に頂点を追加
                uint centerVId = cad2D.AddVertex(CadElementType.Loop, baseLoopId, new OpenTK.Vector2d(x0, y0)).AddVId;
                System.Diagnostics.Debug.Assert(centerVId != 0);
            }
            // ロッドの分割数調整
            for (int k = 1; k < rodRadiusDiv; k++)
            {
                for (int itheta = 0; itheta <= (rodCircleDiv / 2); itheta++)
                {
                    if (Math.Abs(additionalAngle) < IvyFEM.Constants.PrecisionLowerLimit &&
                        (itheta == 0 || itheta == (rodCircleDiv / 2)))
                    {
                        continue;
                    }
                    double theta = 0;
                    if (isReverseAddVertex)
                    {
                        theta = startAngle * Math.PI / 180.0 - itheta * 2.0 * Math.PI / rodCircleDiv;
                    }
                    else
                    {
                        theta = startAngle * Math.PI / 180.0 + itheta * 2.0 * Math.PI / rodCircleDiv;
                    }
                    double x = x0 + (k * rodRadius / rodRadiusDiv) * Math.Cos(theta);
                    double y = y0 + (k * rodRadius / rodRadiusDiv) * Math.Sin(theta);
                    uint id_v_add = cad2D.AddVertex(CadElementType.Loop, baseLoopId, new OpenTK.Vector2d(x, y)).AddVId;
                    System.Diagnostics.Debug.Assert(id_v_add != 0);
                }
            }
            // ロッドの分割数調整: ロッド1/4円から超えた部分
            if (Math.Abs(additionalAngle) >= IvyFEM.Constants.PrecisionLowerLimit)
            {
                double theta = 0;
                if (isReverseAddVertex)
                {
                    theta = (startAngle + additionalAngle) * Math.PI / 180.0;
                }
                else
                {
                    theta = (startAngle - additionalAngle) * Math.PI / 180.0;
                }
                for (int k = 1; k < rodRadiusDiv; k++)
                {
                    double x = x0 + (k * rodRadius / rodRadiusDiv) * Math.Cos(theta);
                    double y = y0 + (k * rodRadius / rodRadiusDiv) * Math.Sin(theta);
                    uint addVId = cad2D.AddVertex(CadElementType.Loop, baseLoopId, new OpenTK.Vector2d(x, y)).AddVId;
                    System.Diagnostics.Debug.Assert(addVId != 0);
                }
            }
            // ロッドの分割数調整: ロッド1/4円から超えた部分
            if (Math.Abs(additionalAngle) >= IvyFEM.Constants.PrecisionLowerLimit)
            {
                double theta = 0;
                if (isReverseAddVertex)
                {
                    theta = (startAngle - 180.0 - additionalAngle) * Math.PI / 180.0;
                }
                else
                {
                    theta = (startAngle + 180.0 + additionalAngle) * Math.PI / 180.0;
                }
                for (int k = 1; k < rodRadiusDiv; k++)
                {
                    double x = x0 + (k * rodRadius / rodRadiusDiv) * Math.Cos(theta);
                    double y = y0 + (k * rodRadius / rodRadiusDiv) * Math.Sin(theta);
                    uint addVId = cad2D.AddVertex(CadElementType.Loop, baseLoopId, new OpenTK.Vector2d(x, y)).AddVId;
                    System.Diagnostics.Debug.Assert(addVId != 0);
                }
            }

            uint retLoopId = 0;
            if (stVId != uint.MaxValue && edVId != uint.MaxValue)
            {
                uint prevVId = stVId;

                // ロッド半円から超えた部分
                if (Math.Abs(additionalAngle) >= IvyFEM.Constants.PrecisionLowerLimit)
                {
                    double theta = 0;
                    if (isReverseAddVertex)
                    {
                        theta = (startAngle + additionalAngle) * Math.PI / 180.0;
                    }
                    else
                    {
                        theta = (startAngle - additionalAngle) * Math.PI / 180.0;
                    }
                    double x = x0 + rodRadius * Math.Cos(theta);
                    double y = y0 + rodRadius * Math.Sin(theta);
                    uint addVId = cad2D.AddVertex(CadElementType.Loop, baseLoopId, new OpenTK.Vector2d(x, y)).AddVId;
                    System.Diagnostics.Debug.Assert(addVId != 0);
                    var connectVertexRes = cad2D.ConnectVertexLine(prevVId, addVId);
                    uint addEId = connectVertexRes.AddEId;
                    System.Diagnostics.Debug.Assert(addEId != 0);
                    prevVId = addVId;
                }

                // ロッド半円
                for (int itheta = 0; itheta <= (rodCircleDiv / 2); itheta++)
                {
                    double theta = 0;
                    if (isReverseAddVertex)
                    {
                        theta = startAngle * Math.PI / 180.0 - itheta * 2.0 * Math.PI / rodCircleDiv;
                    }
                    else
                    {
                        theta = startAngle * Math.PI / 180.0 + itheta * 2.0 * Math.PI / rodCircleDiv;
                    }
                    double x = x0 + rodRadius * Math.Cos(theta);
                    double y = y0 + rodRadius * Math.Sin(theta);
                    uint addVId = cad2D.AddVertex(CadElementType.Loop, baseLoopId, new OpenTK.Vector2d(x, y)).AddVId;
                    System.Diagnostics.Debug.Assert(addVId != 0);
                    var connectVertexRes = cad2D.ConnectVertexLine(prevVId, addVId);
                    uint addEId = connectVertexRes.AddEId;
                    System.Diagnostics.Debug.Assert(addEId != 0);
                    prevVId = addVId;
                }
                // ロッド半円から超えた部分
                if (Math.Abs(additionalAngle) >= IvyFEM.Constants.PrecisionLowerLimit)
                {
                    double theta = 0;
                    if (isReverseAddVertex)
                    {
                        theta = (startAngle - 180.0 - additionalAngle) * Math.PI / 180.0;
                    }
                    else
                    {
                        theta = (startAngle + 180.0 + additionalAngle) * Math.PI / 180.0;
                    }
                    double x = x0 + rodRadius * Math.Cos(theta);
                    double y = y0 + rodRadius * Math.Sin(theta);
                    uint addVId = cad2D.AddVertex(CadElementType.Loop, baseLoopId, new OpenTK.Vector2d(x, y)).AddVId;
                    System.Diagnostics.Debug.Assert(addVId != 0);
                    var connectVertexRes = cad2D.ConnectVertexLine(prevVId, addVId);
                    uint addEId = connectVertexRes.AddEId;
                    System.Diagnostics.Debug.Assert(addEId != 0);
                    prevVId = addVId; //!!!!!!!!!!!!!!!
                }
                uint lastVId = edVId;
                {
                    var resConnectVertex = cad2D.ConnectVertexLine(prevVId, lastVId);
                    uint addEId = resConnectVertex.AddEId;
                    uint lId = resConnectVertex.AddLId;
                    System.Diagnostics.Debug.Assert(addEId != 0);
                    System.Diagnostics.Debug.Assert(lId != 0);
                    retLoopId = lId;
                }
            }
            else
            {
                System.Diagnostics.Debug.Assert(isReverseAddVertex == false); // 逆順未対応
                // ロッド半円から超えた部分
                if (Math.Abs(additionalAngle) >= IvyFEM.Constants.PrecisionLowerLimit)
                {
                    double theta = (startAngle - additionalAngle) * Math.PI / 180.0;
                    double x = x0 + rodRadius * Math.Cos(theta);
                    double y = y0 + rodRadius * Math.Sin(theta);
                    pts.Add(new OpenTK.Vector2d(x, y));
                }
                // ロッド半円
                for (int itheta = 0; itheta <= (rodCircleDiv / 2); itheta++)
                {
                    double theta = startAngle * Math.PI / 180.0 + itheta * 2.0 * Math.PI / rodCircleDiv;
                    double x = x0 + rodRadius * Math.Cos(theta);
                    double y = y0 + rodRadius * Math.Sin(theta);
                    pts.Add(new OpenTK.Vector2d(x, y));
                }
                // ロッド半円から超えた部分
                if (Math.Abs(additionalAngle) >= IvyFEM.Constants.PrecisionLowerLimit)
                {
                    double theta = (startAngle + 180.0 + additionalAngle) * Math.PI / 180.0;
                    double x = x0 + rodRadius * Math.Cos(theta);
                    double y = y0 + rodRadius * Math.Sin(theta);
                    pts.Add(new OpenTK.Vector2d(x, y));
                }
                uint lId = cad2D.AddPolygon(pts, baseLoopId).AddLId;
                retLoopId = lId;
            }
            return retLoopId;
        }

        /// <summary>
        /// 左のロッド
        /// </summary>
        /// <param name="cad2D"></param>
        /// <param name="baseLoopId"></param>
        /// <param name="vId0"></param>
        /// <param name="vId1"></param>
        /// <param name="vId2"></param>
        /// <param name="x0"></param>
        /// <param name="y0"></param>
        /// <param name="rodRadius"></param>
        /// <param name="rodCircleDiv"></param>
        /// <param name="rodRadiusDiv"></param>
        /// <returns></returns>
        public static uint AddLeftRod(
            CadObject2D cad2D, uint baseLoopId,
            uint vId0, uint vId1, uint vId2,
            double x0, double y0, double rodRadius, int rodCircleDiv, int rodRadiusDiv)
        {
            return AddExactlyHalfRod(cad2D, baseLoopId,
                vId0, vId1, vId2,
                x0, y0, rodRadius, rodCircleDiv, rodRadiusDiv,
                90.0, true);
        }

        /// <summary>
        /// 半円（余剰角度なし)ロッドの追加
        /// </summary>
        /// <param name="cad2D"></param>
        /// <param name="baseLoopId"></param>
        /// <param name="vId0"></param>
        /// <param name="vId1"></param>
        /// <param name="vId2"></param>
        /// <param name="x0"></param>
        /// <param name="y0"></param>
        /// <param name="rodRadius"></param>
        /// <param name="rodCircleDiv"></param>
        /// <param name="rodRadiusDiv"></param>
        /// <param name="startAngle"></param>
        /// <param name="isReverseAddVertex"></param>
        /// <returns></returns>
        public static uint AddExactlyHalfRod(
            CadObject2D cad2D, uint baseLoopId,
            uint vId0, uint vId1, uint vId2,
            double x0, double y0, double rodRadius, int rodCircleDiv, int rodRadiusDiv,
            double startAngle, bool isReverseAddVertex)
        {
            uint retLoopId = 0;

            OpenTK.Vector2d pt_center = cad2D.GetVertexCoord(vId1);
            double th = 1.0e-12;
            System.Diagnostics.Debug.Assert(Math.Abs(x0 - pt_center.X) < th);
            System.Diagnostics.Debug.Assert(Math.Abs(y0 - pt_center.Y) < th);
            // check
            //CVector2D pt0 = cad2d.GetVertexCoord(id_v0);
            //double x_pt0 = pt0.x;
            //double y_pt0 = pt0.y;
            //CVector2D pt2 = cad2d.GetVertexCoord(id_v2);
            //double x_pt2 = pt2.x;
            //double y_pt2 = pt2.y;

            // ロッドの分割数調整
            for (int k = 1; k < rodRadiusDiv; k++)
            {
                for (int itheta = 1; itheta < (rodCircleDiv / 2); itheta++)
                {
                    double theta = 0.0;
                    if (isReverseAddVertex)
                    {
                        theta = startAngle * Math.PI / 180.0 - itheta * 2.0 * Math.PI / rodCircleDiv;
                    }
                    else
                    {
                        theta = startAngle * Math.PI / 180.0 + itheta * 2.0 * Math.PI / rodCircleDiv;
                    }
                    double x = x0 + (k * rodRadius / rodRadiusDiv) * Math.Cos(theta);
                    double y = y0 + (k * rodRadius / rodRadiusDiv) * Math.Sin(theta);
                    uint addVId = cad2D.AddVertex(CadElementType.Loop, baseLoopId, new OpenTK.Vector2d(x, y)).AddVId;
                    System.Diagnostics.Debug.Assert(addVId != 0);
                }
            }
            // ロッド
            uint prevVId = vId2;
            for (int itheta = 1; itheta < (rodCircleDiv / 2); itheta++)
            {
                double theta = 0.0;
                if (isReverseAddVertex)
                {
                    theta = startAngle * Math.PI / 180.0 - itheta * 2.0 * Math.PI / rodCircleDiv;
                }
                else
                {
                    theta = startAngle * Math.PI / 180.0 + itheta * 2.0 * Math.PI / rodCircleDiv;
                }
                double x = x0 + rodRadius * Math.Cos(theta);
                double y = y0 + rodRadius * Math.Sin(theta);
                uint addVId = cad2D.AddVertex(CadElementType.Loop, baseLoopId, new OpenTK.Vector2d(x, y)).AddVId;
                System.Diagnostics.Debug.Assert(addVId != 0);
                var connectVertexRes = cad2D.ConnectVertexLine(prevVId, addVId);
                uint addEId = connectVertexRes.AddEId;
                System.Diagnostics.Debug.Assert(addEId != 0);
                prevVId = addVId;
            }
            uint lastVId = vId0;
            {
                var connectVertexRes = cad2D.ConnectVertexLine(prevVId, lastVId);
                uint addEId = connectVertexRes.AddEId;
                uint lId = connectVertexRes.AddLId;
                System.Diagnostics.Debug.Assert(addEId != 0);
                System.Diagnostics.Debug.Assert(lId != 0);
                retLoopId = lId;
            }
            return retLoopId;
        }

        /// <summary>
        /// 右のロッド
        /// </summary>
        /// <param name="cad2D"></param>
        /// <param name="baseLoopId"></param>
        /// <param name="vId0"></param>
        /// <param name="vId1"></param>
        /// <param name="vId2"></param>
        /// <param name="x0"></param>
        /// <param name="y0"></param>
        /// <param name="rodRadius"></param>
        /// <param name="rodCircleDiv"></param>
        /// <param name="rodRadiusDiv"></param>
        /// <returns></returns>
        public static uint AddRightRod(
            CadObject2D cad2D, uint baseLoopId,
            uint vId0, uint vId1, uint vId2,
            double x0, double y0, double rodRadius, int rodCircleDiv, int rodRadiusDiv)
        {
            return AddExactlyHalfRod(cad2D, baseLoopId,
                vId0, vId1, vId2,
                x0, y0, rodRadius, rodCircleDiv, rodRadiusDiv,
                270.0, true);
        }

        /// <summary>
        /// 上のロッド
        /// </summary>
        /// <param name="cad2D"></param>
        /// <param name="baseLoopId"></param>
        /// <param name="vId0"></param>
        /// <param name="vId1"></param>
        /// <param name="vId2"></param>
        /// <param name="x0"></param>
        /// <param name="y0"></param>
        /// <param name="rodRadius"></param>
        /// <param name="rodCircleDiv"></param>
        /// <param name="rodRadiusDiv"></param>
        /// <returns></returns>
        public static uint AddTopRod(
            CadObject2D cad2D, uint baseLoopId,
            uint vId0, uint vId1, uint vId2,
            double x0, double y0, double rodRadius, int rodCircleDiv, int rodRadiusDiv)
        {
            return AddExactlyHalfRod(cad2D, baseLoopId,
                vId0, vId1, vId2,
                x0, y0, rodRadius, rodCircleDiv, rodRadiusDiv,
                0.0, true);
        }

        /// <summary>
        /// 下のロッド
        /// </summary>
        /// <param name="cad2D"></param>
        /// <param name="baseLoopId"></param>
        /// <param name="vId0"></param>
        /// <param name="vId1"></param>
        /// <param name="vId2"></param>
        /// <param name="x0"></param>
        /// <param name="y0"></param>
        /// <param name="rodRadius"></param>
        /// <param name="rodCircleDiv"></param>
        /// <param name="rodRadiusDiv"></param>
        /// <returns></returns>
        public static uint AddBottomRod(
            CadObject2D cad2D, uint baseLoopId,
            uint vId0, uint vId1, uint vId2,
            double x0, double y0, double rodRadius, int rodCircleDiv, int rodRadiusDiv)
        {
            // 注意：id_v0とid_v2が逆になる
            return AddExactlyHalfRod(cad2D, baseLoopId,
                vId2, vId1, vId0,
                x0, y0, rodRadius, rodCircleDiv, rodRadiusDiv,
                180.0, true);
        }


        /// <summary>
        /// ロッド1/4円を追加する
        /// </summary>
        /// <param name="cad2D"></param>
        /// <param name="baseLoopId"></param>
        /// <param name="x0"></param>
        /// <param name="y0"></param>
        /// <param name="rodRadius"></param>
        /// <param name="rodCircleDiv"></param>
        /// <param name="rodRadiusDiv"></param>
        /// <returns></returns>
        public static uint AddQuarterRod(
            CadObject2D cad2D, uint baseLoopId,
            uint vId0, uint vId1, uint vId2,
            double x0, double y0, double rodRadius, int rodCircleDiv, int rodRadiusDiv,
            double startAngle, double endAngle, bool isReverseAddVertex = false)
        {
            OpenTK.Vector2d centerPt = cad2D.GetVertexCoord(vId1);
            System.Diagnostics.Debug.Assert(Math.Abs(x0 - centerPt.X) < IvyFEM.Constants.PrecisionLowerLimit);
            System.Diagnostics.Debug.Assert(Math.Abs(y0 - centerPt.Y) < IvyFEM.Constants.PrecisionLowerLimit);
            IList<OpenTK.Vector2d> pts = new List<OpenTK.Vector2d>();

            // ロッドの分割数調整
            for (int k = 1; k < rodRadiusDiv; k++)
            {
                for (int itheta = 0; itheta <= rodCircleDiv; itheta++)
                {
                    double workAngle = 0.0;
                    if (isReverseAddVertex)
                    {
                        workAngle = 360.0 - itheta * 360.0 / rodCircleDiv;
                        if (workAngle <= endAngle || workAngle >= startAngle)
                        {
                            continue;
                        }
                    }
                    else
                    {
                        workAngle = itheta * 360.0 / rodCircleDiv;
                        if (workAngle <= startAngle || workAngle >= endAngle)
                        {
                            continue;
                        }
                    }
                    double theta = workAngle * Math.PI / 180.0;
                    double x = x0 + (k * rodRadius / rodRadiusDiv) * Math.Cos(theta);
                    double y = y0 + (k * rodRadius / rodRadiusDiv) * Math.Sin(theta);
                    uint addVId = cad2D.AddVertex(CadElementType.Loop, baseLoopId, new OpenTK.Vector2d(x, y)).AddVId;
                    System.Diagnostics.Debug.Assert(addVId != 0);
                }
            }
            uint retLoopId = 0;
            uint prevVId = vId0;

            // ロッド1/4円
            for (int itheta = 0; itheta <= rodCircleDiv; itheta++)
            {
                double workAngle = 0.0;
                if (isReverseAddVertex)
                {
                    workAngle = 360.0 - itheta * 360.0 / rodCircleDiv;
                    if (workAngle <= endAngle || workAngle >= startAngle)
                    {
                        continue;
                    }
                }
                else
                {
                    workAngle = itheta * 360.0 / rodCircleDiv;
                    if (workAngle <= startAngle || workAngle >= endAngle)
                    {
                        continue;
                    }
                }

                double theta = workAngle * Math.PI / 180.0;
                double x = x0 + rodRadius * Math.Cos(theta);
                double y = y0 + rodRadius * Math.Sin(theta);
                uint addVId = cad2D.AddVertex(CadElementType.Loop, baseLoopId, new OpenTK.Vector2d(x, y)).AddVId;
                System.Diagnostics.Debug.Assert(addVId != 0);
                var connectVertexRes = cad2D.ConnectVertexLine(prevVId, addVId);
                uint addEId = connectVertexRes.AddEId;
                System.Diagnostics.Debug.Assert(addEId != 0);
                prevVId = addVId;
            }
            uint lastVId = vId2;
            {
                var connectVertexRes = cad2D.ConnectVertexLine(prevVId, lastVId);
                uint addEId = connectVertexRes.AddEId;
                uint lId = connectVertexRes.AddLId;
                System.Diagnostics.Debug.Assert(addEId != 0);
                System.Diagnostics.Debug.Assert(lId != 0);
                retLoopId = lId;
            }

            return retLoopId;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="cad2D"></param>
        /// <param name="baseLoopId"></param>
        /// <param name="x0"></param>
        /// <param name="y0"></param>
        /// <param name="rodRadius"></param>
        /// <param name="rodCircleDiv"></param>
        /// <param name="rodRadiusDiv"></param>
        /// <param name="vId0"></param>
        /// <param name="vId1"></param>
        /// <param name="vId2"></param>
        /// <param name="startAngle"></param>
        /// <param name="isReverseAddVertex"></param>
        /// <returns></returns>
        public static uint AddExactlyQuarterRod(
            CadObject2D cad2D,
            uint baseLoopId, double x0, double y0, double rodRadius, int rodCircleDiv, int rodRadiusDiv,
            uint vId0, uint vId1, uint vId2, double startAngle, bool isReverseAddVertex)
        {
            OpenTK.Vector2d centerPt = cad2D.GetVertexCoord(vId1);
            System.Diagnostics.Debug.Assert(Math.Abs(x0 - centerPt.X) < IvyFEM.Constants.PrecisionLowerLimit);
            System.Diagnostics.Debug.Assert(Math.Abs(y0 - centerPt.Y) < IvyFEM.Constants.PrecisionLowerLimit);
            IList<OpenTK.Vector2d> pts = new List<OpenTK.Vector2d>();
            System.Diagnostics.Debug.Assert(
                (startAngle == 0.0) || (startAngle == 90.0) || (startAngle == 180.0) || (startAngle == 270.0));

            // ロッドの分割数調整
            for (int k = 1; k < rodRadiusDiv; k++)
            {
                for (int itheta = 1; itheta < (rodCircleDiv / 4); itheta++)
                {
                    double theta = 0;
                    if (isReverseAddVertex)
                    {
                        theta = startAngle * Math.PI / 180.0 - itheta * 2.0 * Math.PI / rodCircleDiv;
                    }
                    else
                    {
                        theta = startAngle * Math.PI / 180.0 + itheta * 2.0 * Math.PI / rodCircleDiv;
                    }
                    double x = x0 + (k * rodRadius / rodRadiusDiv) * Math.Cos(theta);
                    double y = y0 + (k * rodRadius / rodRadiusDiv) * Math.Sin(theta);
                    uint addVId = cad2D.AddVertex(CadElementType.Loop, baseLoopId, new OpenTK.Vector2d(x, y)).AddVId;
                    System.Diagnostics.Debug.Assert(addVId != 0);
                }
            }
            uint retLoopId = 0;
            uint prevVId = vId0;

            // ロッド1/4円
            for (int itheta = 1; itheta < (rodCircleDiv / 4); itheta++)
            {
                double theta = 0;
                if (isReverseAddVertex)
                {
                    theta = startAngle * Math.PI / 180.0 - itheta * 2.0 * Math.PI / rodCircleDiv;
                }
                else
                {
                    theta = startAngle * Math.PI / 180.0 + itheta * 2.0 * Math.PI / rodCircleDiv;
                }
                double x = x0 + rodRadius * Math.Cos(theta);
                double y = y0 + rodRadius * Math.Sin(theta);
                uint addVId = cad2D.AddVertex(CadElementType.Loop, baseLoopId, new OpenTK.Vector2d(x, y)).AddVId;
                System.Diagnostics.Debug.Assert(addVId != 0);
                var connectVertexRes = cad2D.ConnectVertexLine(prevVId, addVId);
                uint addEId = connectVertexRes.AddEId;
                System.Diagnostics.Debug.Assert(addEId != 0);
                prevVId = addVId;
            }
            uint lastVId = vId2;
            {
                var connectVertexRes = cad2D.ConnectVertexLine(prevVId, lastVId);
                uint addEId = connectVertexRes.AddEId;
                uint lId = connectVertexRes.AddLId;
                System.Diagnostics.Debug.Assert(addEId != 0);
                System.Diagnostics.Debug.Assert(lId != 0);
                retLoopId = lId;
            }

            return retLoopId;
        }

        /// <summary>
        /// Y軸からの回転移動原点、角度を算出する
        /// </summary>
        /// <param name="world"></param>
        /// <param name="bcEIds"></param>
        /// <param name="rotAngle"></param>
        /// <param name="rotOrigin"></param>
        public static void GetRotOriginRotAngleFromY(
            FEWorld world, uint[] bcEIds, out double rotAngle, out double[] rotOrigin)
        {
            uint eId1 = bcEIds[0];
            uint eId2 = bcEIds[bcEIds.Length - 1];
            Edge2D e1 = world.Mesh.Cad2D.GetEdge(eId1);
            Edge2D e2 = world.Mesh.Cad2D.GetEdge(eId2);
            OpenTK.Vector2d firstPt = e1.GetVertexCoord(true);
            OpenTK.Vector2d lastPt = e2.GetVertexCoord(false);
            double[] firstCoord = { firstPt.X, firstPt.Y };
            double[] lasCoord = { lastPt.X, lastPt.Y };

            // 回転移動
            rotOrigin = null;
            rotAngle = 0.0; // ラジアン
            // 境界の傾きから回転角度を算出する
            if (Math.Abs(firstCoord[0] - lasCoord[0]) >= 1.0e-12)
            {
                // X軸からの回転角
                rotAngle = Math.Atan2((lasCoord[1] - firstCoord[1]), (lasCoord[0] - firstCoord[0]));
                // Y軸からの回転角に変換 (境界はY軸に平行、X方向周期構造)
                rotAngle = rotAngle - 0.5 * Math.PI;
                rotOrigin = firstCoord;
                System.Diagnostics.Debug.WriteLine("rotAngle: {0} rotOrigin:{1} {2}", rotAngle * 180.0 / Math.PI, rotOrigin[0], rotOrigin[1]);
            }
            else
            {
                // 角度0ラジアン
                rotAngle = 0.0;
                rotOrigin = null;
            }
        }
    }
}
