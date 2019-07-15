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

#molオブジェクトをまとめてsdfに
WriFP=DataDir+'PCCDB.sdf'

L_mol=[]
with open(WriFP,'rb') as f:
    suppl=Chem.ForwardSDMolSupplier(f,removeHs=False)
    for el in suppl:
        if el is not None:
            L_mol.append(el)

print(Chem.MolToMolBlock(L_mol[0]))





mol=L_mol[0]
for a in mol.GetAtoms():
    print(a.)

#まとめたSDFをpandasにロード
df_PCCDB=PandasTools.LoadSDF(WriFP,embedProps=True)
da.descDf(df_PCCDB)

mol=df_PCCDB.iloc[6,3]

print(Chem.MolToMolBlock(mol))

for a in mol.GetAtoms():
    print(a.GetSymbol())


