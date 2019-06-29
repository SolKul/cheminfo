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

# +
import pandas as pd
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

#ダウンロード可能な分子のリストのCSVをまとめる
csv_dir='./csv/'
l_csv=[]
for f in os.listdir(csv_dir):
    d_csv={}
    if f.endswith('.csv'):
        d_csv['FN']=f
        d_csv['RP']=os.path.join(csv_dir,f)
        l_csv.append(d_csv)
print(len(l_csv))

l_csv

l_df_id=[]
for d_csv in l_csv:
    l_df_id.append(pd.read_csv(d_csv['RP']))
df_id=pd.concat(l_df_id)
df_id.drop_duplicates(inplace=True)
df_id.reset_index(inplace=True)
da.descDf(df_id)

#csvがおかしいところがあり、PCCDB-IDに数字が入ってないことがあるので、
#それ以外は全部int型に
#変換できないものは残る
df_id_int=df_id.copy()
for i in range(df_id.shape[0]):
    try:
        df_id_int.iloc[i,1]=int(df_id.iloc[i]['PCCDB-ID'])
    except:
        print('Error')

#型がintのものだけを残す
df_id_num=df_id_int[df_id_int['PCCDB-ID'].map(type)==int]

df_id_num.to_pickle('./df_pccdb_id.pkl')

df_id_num.sort_values(by=['PCCDB-ID'])

df_id.iloc[47960]

s_test='328df64'
int(s_test)
