% IDWT2   単一レベルの逆離散2次元ウェーブレット変換
%
% IDWT2 は、('wname',WFILTERS を参照)で指定された特定のウェーブレット
% または指定された特定のウェーブレット再構成フィルタ(Lo_R と Hi_R)の
% いずれかに関連した単一レベルの2次元ウェーブレット再構成を行います。
%
% X = IDWT2(CA,CH,CV,CD,'wname') は、ウェーブレット 'wname' を使って、
% Approximation 係数ベクトル CA と Detail 行列 CH, CV, CD(水平、垂直、対角)を
% ベースに、単一レベルの再構成 Approximation 係数行列 X を計算します。
%
% X = IDWT2(CA,CH,CV,CD,Lo_R,Hi_R) は、設定したフィルタを使って上述の
% 処理を行います。
%   Lo_R は、再構成ローパスフィルタです。
%   Hi_R は、再構成ハイパスフィルタです。
%   Lo_R と Hi_R は、同じ長さでなければなりません。
%
% SA = size(CA) = size(CH) = size(CV) = size(CD) で、LF がフィルタの長さの場合、
% size(X) = SX  が成立します。ここで、DWT 拡張モードが周期的なモードである場合、
% SX = 2*SA です。その他の拡張モードの場合、SX = 2*SA-LF+2 です。
%
% X = IDWT2(CA,CH,CV,CD,'wname',S)、または、X = IDWT2(CA,CH,CV,CD,Lo_R,Hi_R,S) 
% は、IDWT2(CA,CH,CV,CD,'wname') を使って得られる結果の中心部のサイズ S の部分を
% 出力します。S は、SX よりも小さくなければなりません。
%
% X = IDWT2(CA,CH,CV,CD,'wname','mode',mode) 
% X = IDWT2CA,CH,CV,CD,Lo_R,Hi_R,'mode',mode) 
% X = IDWT2(CA,CH,CV,CD,'wname',L,'mode',mode) 
% X = IDWT2(CA,CH,CV,CD,Lo_R,Hi_R,L,'mode',mode) のいずれも、指定可能な
% 拡張モードで、ウェーブレット再構成の係数を計算します。
%
% X = IDWT2(CA,[],[],[]、... ) は、Approximation 係数行列 CA をベースに、
% 単一レベルに再構成された Approximation 係数行列 X を出力します。
%
% X = IDWT2([],CH,[],[]、... ) は、水平の Detail 係数行列 CH をベースに、
% 単一レベルに再構成された Detail 係数行列 X を出力します。
%
% X = IDWT2([],[],CV,[]、... ) と X = IDWT2([],[],[],CD、... )は、
% 同じ結果を出力します。
% 
%   参考 DWT2, DWTMODE, UPWLEV2.


%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision: 19-May-2003.
%   Copyright 1995-2004 The MathWorks, Inc.
