# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.2
#   kernelspec:
#     display_name: Python cheminfo
#     language: python
#     name: cheminfo
# ---

import pubchempy as pcp
from rdkit import Chem
import py3Dmol as p3d
import os
from rdkit.Chem import AllChem, Draw

TarDir='./SDF_R/'
L_FP=[]
for f in os.listdir(TarDir):
    D_sdf={}
    if f.endswith('.sdf'):
        D_sdf['FN']=f
        D_sdf['RP']=os.path.join(TarDir,f)
        L_FP.append(D_sdf)
len(L_FP)

L_sdf=[]
for D_sdf in L_FP:
    with open(D_sdf['RP'],'rb') as f:
        suppl=Chem.ForwardSDMolSupplier(f)
        D_sdf_all=D_sdf.copy()
        D_sdf_all['sdf']=f.read().decode()
        L_sdf.append(D_sdf_all)
len(L_sdf)
