# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.2
#   kernelspec:
#     display_name: Python scrap
#     language: python
#     name: scrap
# ---

import requests

UA='Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/74.0.3729.157 Safari/537.36'
headers = {"User-Agent": UA}
resp = requests.get('https://news.yahoo.co.jp/', timeout=1, headers=headers)

print(resp.text)

resp = requests.get('http://pccdb.org/search_pubchemqc/get_sdf/ver0.2/1', timeout=120, headers=headers)

import pandas as pd

df_pccd=pd.read_csv('./04ReqSDF/search_result_ids.csv')

TarDir='./04ReqSDF/'
FN='PCCID{}'.format(1)
FP=TarDir+FN+'.sdf'

with open(FP, mode='w') as f:
    f.write(resp.text)

df_pccd.iloc[1,0]
