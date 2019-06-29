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

import pandas as pd
import os
import os.path

#保存先の設定
DataDir=os.path.join(os.environ['HOME'],'notebooks/50Data/')

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

l_df_id=[]
for d_csv in l_csv:
    l_df_id.append(pd.read_csv(d_csv['RP']))
df_id=pd.concat(l_df_id)
df_id.drop_duplicates(inplace=True)
df_id.reset_index(inplace=True)
