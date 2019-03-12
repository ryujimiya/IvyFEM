IvyFEM  
======  

有限要素法(FEM)を用いたCAEライブラリです。  
.Net WPFアプリケーションで使用することを前提としています。  
　  
IvyFEMは[DelFEM4Net](https://code.google.com/p/delfem4net/) の後継ライブラリとなります。(DelFEM4Netは開発中止)  
　  
　  
　**いまできること（2019-03-12更新）**  
　  
　  ☑ 単純な2D(ポリゴン)の図面作成  
　  ☑ 有限要素（三角形要素）分割  
　  ☑ 有限要素行列の作成 （*1）  
　  ☑ リニアシステムを解く（LAPACKE, Lis、独自実装）（*1）  
　  ☑ サーモグラフィーのような分布図  
　  
　  *1 いま用意しているのは  
　　　　　 電磁気学：H面導波管の伝達問題  
　　　　　 力学： 弾性体の構造解析  
　　　　　　　　　  線形弾性体  
　　　　　　　　　  超弾性体  
　　　　　　　　　　  Saint Venant Kirchhoff  
　　　　　　　　　　  Mooney-Rivlin (非圧縮、微圧縮)  
　　　　　　　　　　  Ogden (非圧縮、微圧縮)  
　　　　　　　　　  多点拘束(Multipoint Constraint, MPC)(直線)  
　　　　　　　　　  剛体との接触(直線、円)  
　　　　　　　　　  弾性体二体接触  
　　　 これからいろんな問題を解けるようにしていきます。  
　  
　**バイナリ**  
　  
　  IvyFEM.dll version 0.0.0.1  
　  ※開発途中なのでDebug版しかありません。  
　  
　  [IvyFEM](https://github.com/ryujimiya/IvyFEM/blob/master/publish/)  
　  
　**使い方**  
　  
　  IvyFEMライブラリを使ったサンプルアプリケーションIvyFEMProtoAppをご参照ください。  
　  
　[IvyFEMProtoApp](https://github.com/ryujimiya/IvyFEMProtoApp/)  
　  
![IvyFEMProtoAppアプリ](https://pbs.twimg.com/media/DjHvvKfUcAEMU_H.jpg)  
　  
　  
　  
