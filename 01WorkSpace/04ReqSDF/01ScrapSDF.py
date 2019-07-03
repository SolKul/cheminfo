# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # PCCDBからSDFをダウンロード  
# PCCDBからjson形式でSDFがダウンロード可能  

# +
import requests
import pandas as pd
import json
import time
import os
import os.path
import sys

#自作モジュールのインポート
p_mod=os.path.join(os.environ['HOME'],'notebooks/99MyModules')
sys.path.append(p_mod)
import descarray as da

#保存先の設定
DataDir=os.path.join(os.environ['HOME'],'notebooks/50Data/')
# -

df_id=pd.read_pickle('./df_pccdb_id.pkl')
da.descDf(df_id)

UA='Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/74.0.3729.157 Safari/537.36'
headers = {"User-Agent": UA}

RowPC=0
PC_id=int(df_id.iloc[RowPC,1])
PCCDURL='http://pccdb.org/search_pubchemqc/get_sdf/ver0.2/'
TarURL=PCCDURL+str(PC_id)
resp = requests.get(TarURL, timeout=120, headers=headers)

sdf_txt=resp.text
di_sdf=json.loads(sdf_txt)
print([x for x in di_sdf])
TarDir=DataDir+'PCCDB_SDF/'
FN='PCCID_{:0=8}'.format(PC_id)
FP=TarDir+FN+'.sdf'
with open(FP, mode='w') as f:
    f.write(di_sdf['sdf'])

for RowPC in range(47984):
    #RowPCは何行目かを表す。
    #PC_idはPCCDBのid
    PC_id=int(df_id.iloc[RowPC,1])

    PCCDURL='http://pccdb.org/search_pubchemqc/get_sdf/ver0.2/'
    TarURL=PCCDURL+str(PC_id)
    print(TarURL)
    try:
        time.sleep(60)
        resp = requests.get(TarURL, timeout=120, headers=headers)

        #requestsしたtext、これはjson形式になっている
        #属性はpcccdb_idとsdfの2つ
        #json.loadsすると辞書型に変換される。
        sdf_txt=resp.text
        di_sdf=json.loads(sdf_txt)
        print([x for x in di_sdf])

        #ファイル名とファイルパスを作成
        TarDir=DataDir+'PCCDB_SDF/'
        FN='PCCID_{:0=8}'.format(PC_id)
        FP=TarDir+FN+'.sdf'
        print(FP)
        print()

        with open(FP, mode='w') as f:
            f.write(di_sdf['sdf'])
    except:
        print('Error')

# +
# import re

# pattern = 'sdf'
# re.search(sdf_txt,pattern)

# n_txt=eval('{}'.format(sdf_txt))

# tes_txt='000\\\\n aa'
# print(tes_txt)
# te_txt2=tes_txt.replace('\\\\','\\')
# print(te_txt2)

# x = 'a\\b\\c'
# print(x) # a\b\\c
# y = x.replace('\\','')
# print(y) # a\b\c

# import json

# di_sdf=json.loads(sdf_txt)
# print([x for x in di_sdf])

# FP

# df_pccd.iloc[1,0]
