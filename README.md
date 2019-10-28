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
　　　　　　　　　　多点拘束(Multipoint Constraint [MPC])(直線)  
　　　　　　　　　　剛体との接触(直線、円)  
　　　　　　　　　　弾性体二体接触※  
　　　　流体力学　：Navier-Stokesの方程式  
　　　　　　　　　　　標準  
　　　　　　　　　　　SUPG [Streamline Upwind Petrov-Galerkin]安定化  
　　　　　　　　　　Vorticity / Stream Funciton定式化  
　　　　　　　　　　　標準  
　　　　　　　　　　　SUPG[Streamline Upwind Petrov-Galerkin]安定化※  
　　　　　　　　　　　分離型解法(Runge-Kutta)  
　　　　　　　　　　Pressure Poisson定式化※  
　　　　　　　　　　　標準  
　　　　　　　　　　　分離型解法(Runge-Kutta)  
　　　　電磁気学　：電磁導波路の伝達問題  
　　　　　　　　　　　H面TEモード/E面TEモード/2次元TMモードの導波路  
　　　　　　　　　　　周波数領域FEM  
　　　　　　　　　　　　（固有モード展開）  
　　　　　　　　　　　　（Givoli-Neta-Patlashenkoの高次ABC）（Evanescent Wave ABC）  
　　　　　　　　　　　　（PML）  
　　　　　　　　　　　時間領域FEM  
　　　　　　　　　　　　（Givoli-Neta-Patlashenkoの高次ABC）（Evanescent Wave ABC）  
　　　　　　　　　　　　（Jiao-Jin-Michielssen-RileyのPML）  
　　　　　　　　　　フォトニック結晶導波路の伝達問題  
　　　　　　　　　　　正方格子、三角形格子  
　　　　　　　　　　　誘電体ロッド型、エアホール型  
　　　　　　　　　　　周波数領域FEM  
　　　　　　　　　　　　（固有モード展開）  
　　　　各種方程式：Poisson方程式  
　　　　　　　　　　熱拡散方程式  
　　　　　　　　　　移流拡散方程式  
　　　　　　　　　　Helmholtz方程式  
　　　　三角形要素：Lagrange (1次、2次)  
　　　　　　　　　　Bell Triangle(5次)  
　　　　※印：実験的または未完  
　  
　**バイナリ（2019-10-28更新）**  
　  
　　IvyFEM.dll version 0.0.0.25  
　  
　　**プラットフォームターゲット:　x64**  
　　[IvyFEM](https://github.com/ryujimiya/IvyFEM/blob/master/publish/)  
　  
　**依存ライブラリ**  
　  
　　OpenTK.GLControlをNuGetでインストールしてください。  
　  
　**インストールおよび使い方**  
　  
　　インストール、実装概略を次のページにまとめました。  
　　[.NET向けCAEライブラリIvyFEMを用いて弾性体の曲げの有限要素法シミュレーションをする](https://qiita.com/ryujimiya2361/items/a573ee7d7060a576f304)  
　　[.NET向け有限要素法CAEライブラリIvyFEMでカスタマイズ方程式を実装する](https://qiita.com/ryujimiya2361/items/f003c10cfc222378ad5a)  
　  
　**サンプルアプリケーション**  
　  
　　IvyFEM.dllを使ったサンプルアプリケーション  
　　[IvyFEMProtoApp](https://github.com/ryujimiya/IvyFEMProtoApp/)  
　  
　**計算例**  
　　〇メッシュ  
　　![メッシュ](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190421/20190421122831.jpg)  
　　〇弾性体  
　　![弾性体](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190421/20190421123314.jpg)  
　　〇弾性体多点拘束(MPC)  
　　![弾性体多点拘束(MPC)](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190421/20190421123915.jpg)  
　　〇直線、円との接触問題  
　　![直線、円との接触問題](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190421/20190421124711.jpg)  
　　〇弾性体二体接触問題  
　　![弾性体二体接触問題](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190421/20190421125420.jpg)  
　　〇Poissonの方程式  
　　![Poisson](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190408/20190408221503.jpg)  
　　〇Helmholtzの方程式  
　　![Helmholtz](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190408/20190408221935.jpg)  
　　〇拡散方程式  
　　![拡散](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190408/20190408222500.jpg)  
　　〇移流拡散方程式  
　　![移流拡散](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190408/20190408222735.jpg)  
　　〇流体(Navier-Stokesの方程式)  
　　Standard Galerkin Cavity  
　　![Standard Galerkin Cavity](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190428/20190428094803.jpg)  
　　SUPG Cavity  
　　![SUPG Cavity](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190512/20190512174448.jpg) 
　　SUPG Back-step  
　　![SUPG Back-step 1](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190430/20190430193941.jpg)  
　　![SUPG Back-step 2](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190430/20190430193957.jpg)  
　　![SUPG Back-step 3](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190430/20190430194019.jpg)  
　　![SUPG Back-step 4](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190430/20190430194032.jpg)  
　　![SUPG Back-step 5](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190430/20190430194029.jpg)  
　　![SUPG Back-step 6](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190430/20190430194025.jpg)  
　　![SUPG Back-step 7](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190430/20190430194022.jpg)  
　　〇電磁導波路の周波数領域FEM(固有モード展開)  
　　H面導波路直角コーナーベンド  
　　![ベンド](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190421/20190421130127.jpg)  
　　遮蔽形誘電体導波路の終端  
　　![終端](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20191014/20191014221852.jpg)  
　　〇電磁導波路の周波数領域FEM(高次ABC)  
　　開放形誘電体導波路の終端  
　　![終端](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20191020/20191020134009.jpg)  
　　〇電磁導波路の時間領域FEM(高次ABC)  
　　H面導波路直角コーナーベンド  
　　![ベンド](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190920/20190920201958.jpg)  
　　開放形誘電体導波路の終端  
　　![終端1](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20191014/20191014224143.jpg)  
　　![終端2](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20191014/20191014224410.jpg)  
　　![終端3](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20191014/20191014224518.jpg)  
　　![終端4](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20191014/20191014224619.jpg)  
　　〇フォトニック結晶導波路の周波数領域FEM(固有モード展開)  
　　![正方格子ベンド](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190914/20190914205034.jpg)  
　　![三角形格子ベンド](https://cdn-ak.f.st-hatena.com/images/fotolife/r/ryujimiya/20190914/20190914205057.jpg)  
　  
