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
#正規化したSDFをmolで読み込み、再度SDF化しpandasで読み込み

# +
## 1. ライブラリのインポート
from rdkit import rdBase, Chem
from rdkit.Chem import AllChem, Draw, PandasTools, Descriptors
 
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
print(rdBase.rdkitVersion) # 2019.03.2

# +
import os
import os.path
import sys

#自作モジュールディレクトリをパスに追加
p_mod=os.path.join(os.environ['HOME'],'notebooks/99MyModules')
sys.path.append(p_mod)

import descarray as da

#データディレクトリ
DataDir=os.path.join(os.environ['HOME'],'notebooks/50Data/')
# -

#正規化SDFをロード
TarDir=DataDir+'SDF_R/'
L_FP=[]
for f in os.listdir(TarDir):
    D_sdf={}
    if f.endswith('.sdf'):
        D_sdf['FN']=f
        D_sdf['RP']=os.path.join(TarDir,f)
        L_FP.append(D_sdf)
len(L_FP)

#sdfをmolオブジェクトに
L_sdf=[]
L_mol=[]
for D_sdf in L_FP:
    with open(D_sdf['RP'],'rb') as f:
        suppl=Chem.ForwardSDMolSupplier(f,removeHs=False)
        for el in suppl:
            if el is not None:
                L_mol.append(el)
len(L_mol)

#molオブジェクトをまとめてsdfに
WriFP=DataDir+'PCCDB.sdf'
writer = Chem.SDWriter(WriFP)
for el in L_mol:
    writer.write(el)
writer.close()

#まとめたSDFをpandasにロード
df_PCCDB=PandasTools.LoadSDF(WriFP)
da.descDf(df_PCCDB)



Chem.MolToCXSmiles(df_PCCDB['ROMol'][3],isomericSmiles=False)

print(Chem.MolToMolBlock(L_mol[0]))

with open('./SDF/PCCID_00000003.sdf','rb') as f:
    suppl=Chem.ForwardSDMolSupplier(f)
    for el in suppl:
        if el is not None:
            print(el)
