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

import requests
import pandas as pd
import json
import time

# +
df_pccd=pd.read_csv('./search_result_ids.csv')
RowPC=1

UA='Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/74.0.3729.157 Safari/537.36'
headers = {"User-Agent": UA}
# -

for RowPC in range(100,1000):
    PC_id=int(df_pccd.iloc[RowPC,0])

    PCCDURL='http://pccdb.org/search_pubchemqc/get_sdf/ver0.2/'
    TarURL=PCCDURL+str(PC_id)
    print(TarURL)
    try:
        time.sleep(60)
        resp = requests.get(TarURL, timeout=120, headers=headers)

        sdf_txt=resp.text
        di_sdf=json.loads(sdf_txt)
        print([x for x in di_sdf])

        TarDir='./SDF/'
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
