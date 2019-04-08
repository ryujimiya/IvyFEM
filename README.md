IvyFEM  
======  

有限要素法(FEM)を用いたCAEライブラリです。  
.Net WPFアプリケーションで使用することを想定しています。  
　  
　**いまできること**  
　  
　　☑ 単純な2D(ポリゴン)の図面作成  
　　☑ 有限要素（三角形要素）分割  
　　☑ 有限要素行列の作成 （*1）  
　　☑ リニアシステムを解く（LAPACKE, Lis、独自実装）  
　　☑ サーモグラフィーのような分布図  
　  
　　*1 いま用意しているのは  
　　　　弾性体力学：線形弾性体  
　　　　　　　　　　超弾性体  
　　　　　　　　　　Saint Venant Kirchhoff  
　　　　　　　　　　Mooney-Rivlin (非圧縮、微圧縮)  
　　　　　　　　　　Ogden (非圧縮、微圧縮)  
　　　　　　　　　　多点拘束(Multipoint Constraint, MPC)(直線)  
　　　　　　　　　　剛体との接触(直線、円)  
　　　　　　　　　　弾性体二体接触  
　　　　流体力学　：Navier-Stokesの方程式  
　　　　電磁気学　：H面導波管の伝達問題  
　　　　各種方程式：Poisson方程式  
　　　　　　　　　　熱拡散方程式  
　　　　　　　　　　移流拡散方程式  
　　　　　　　　　　Helmholtz方程式  
　  
　**バイナリ（2019-04-08更新）**  
　  
　　IvyFEM.dll version 0.0.0.5  
　　※開発途中なのでDebug版しかありません。  
　  
　　**プラットフォームターゲット:　x64**  
　　[IvyFEM](https://github.com/ryujimiya/IvyFEM/blob/master/publish/)  
　  
　**使い方**  
　  
　　IvyFEMライブラリを使ったサンプルアプリケーションIvyFEMProtoAppをご参照ください。  
　　[IvyFEMProtoApp](https://github.com/ryujimiya/IvyFEMProtoApp/)  
　  
　　また、実装概略を次のページにまとめました。  
　　[.NET向けCAEライブラリIvyFEMを用いて弾性体の曲げの有限要素法シミュレーションをする](https://qiita.com/ryujimiya2361/items/a573ee7d7060a576f304)  
　  
![IvyFEMProtoAppアプリ](https://pbs.twimg.com/media/DjHvvKfUcAEMU_H.jpg)  
　  
　  
　  
