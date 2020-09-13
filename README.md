# Euler_implicit_LU-SGS

sod,t=0.2.csv：0.2sにおけるsodの衝撃波の結果

scatter_kakunin.py：sodの衝撃波結果との比較コード

sod(FVS,muscl)_im_LU-SGS.jl：juliaの解析コード

sod(FVS,muscl)_im_LU-SGS.py：pythonの解析コード

・使い方

1 juliaまたはpythonの解析コードを実行する

2 scatter_kakunin.pyとsod,t=0.2.csvを結果フォルダ内に置く

3 scatter_kakunin.pyを実行すると解析結果と比較した画像が出力される
