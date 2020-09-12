# Euler_implicit_LU-SGS

sod,t=0.2.csv：0.2sにおけるsodの衝撃波の結果

scatter_kakunin.py：sodの衝撃波結果との比較コード

FVM(FVS,muscl)_in_LDU.jl：juliaの解析コード

FVM(FVS,t(2nd),muscl)_in_LDU.py：pythonの解析コード

・使い方

1 juliaまたはpythonの解析コードを実行する

2 scatter_kakunin.pyとsod,t=0.2.csvを結果フォルダ内に置く

3 scatter_kakunin.pyを実行すると解析結果と比較した画像が出力される
