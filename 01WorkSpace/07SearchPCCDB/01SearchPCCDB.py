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

p_mod=os.path.join(os.environ['HOME'],'notebooks/99MyModules')
sys.path.append(p_mod)

import descarray as da

DataDir=os.path.join(os.environ['HOME'],'notebooks/50Data/')
# -

df_AMES=pd.read_csv(DataDir+'ci900161g_si_001/smiles_cas_N6512.smi',sep='\t', header=None)
df_AMES.columns = ['smiles', 'CAS_NO', 'activity']
da.descDf(df_AMES)

# +
TarDir=DataDir+'SDF_R'

L_FP=[]
for f in os.listdir(TarDir):
    D_sdf={}
    if f.endswith('.sdf'):
        D_sdf['FN']=f
        D_sdf['RP']=os.path.join(TarDir,f)
        L_FP.append(D_sdf)
len(L_FP)
# -

D_sdf=L_FP[1]
l_m=[]
for D_sdf in L_FP:
    with open(D_sdf['RP'],'rb') as f:
        suppl=Chem.ForwardSDMolSupplier(f)
    #     print(f.read().decode('utf8'))
        for m in suppl:
            if m is not None:
                l_m.append(m)
    #             print(list(m.GetPropNames()))

s_allsdf=''
for i in range(1,10):
    D_sdf=L_FP[i]
    with open(D_sdf['RP'],'rb') as f:
        s_allsdf=s_allsdf+f.read().decode()

WriFN='01PCCDBall.sdf'
with open(WriFN,'w') as f:
    f.write(s_allsdf)



# D_sdf=L_FP[1]
PandasTools.LoadSDF(WriFN)

D_sdf['RP']

for i in range(10):
    display(Draw.MolToImage(l_m[i]))

from rdkit.Chem.Draw import IPythonConsole

Tes_F='./PCCID_00000004_NS.sdf'
with open(Tes_F,'rb') as f:
    suppl=Chem.ForwardSDMolSupplier(f)
#     print(f.read().decode('utf8'))
    for m in suppl:
        if m is not None:
            l_PN=list(m.GetPropNames())

l_PN

D_sdf['RP']
