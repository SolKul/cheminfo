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

import pubchempy as pcp
from rdkit import Chem
import py3Dmol as p3d
import os
from rdkit.Chem import AllChem, Draw

TarDir='./SDF/'

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
        D_sdf_all=D_sdf.copy()
        D_sdf_all['sdf']=f.read().decode()
        L_sdf.append(D_sdf_all)
len(L_sdf)

N_lim=10
L_mol=list(range(40,50))
view = p3d.view(width=300, height=250*N_lim, linked=False,viewergrid=(N_lim,1))
for i in range(N_lim):    
    j=L_mol[i]
    view.addModel(L_sdf[j]['sdf'],'sdf', viewer=(i,0))
    view.setStyle({'stick': {'linewidth': 5}}, viewer=(i,0))
    view.zoomTo()
view.show()

# +
TarDir2='./SDF_R/'

for i in range(len(L_sdf)):
    NRP=TarDir2+os.path.splitext(L_sdf[i]['FN'])[0]+'_NS.sdf'
    with open(NRP, mode='w') as f:
        f.write(L_sdf[i]['sdf'].replace(' \n','\n'))
# -

#改行前の空白があるとRdkitで読み込めないので削除
SDF_R=L_sdf[1].replace(' \n','\n')
# print(L_sdf[0].replace(' \n','\n'))

testF='./PCCID0001_R.sdf'
with open(testF, mode='w') as f:
    f.write(SDF_R)

with open('./PCCID0001_R.sdf','rb') as f:
    suppl=Chem.ForwardSDMolSupplier(f)
    for m in suppl:
        print(list(m.GetPropNames()))

with open('./PCCID0001_R.sdf','rb') as f:
    suppl=Chem.ForwardSDMolSupplier(f)
    for m in suppl:
        print(m)

L_sdf[0].split()

L_FP[3]['RP']

with open('./PCCID_00000003.sdf','rb') as f:
    Suppl=Chem.ForwardSDMolSupplier(f,removeHs=False)
#     print(f.read().decode())
    for m in Suppl:
        print(list(m.GetPropNames()))

print(Chem.MolToMolBlock(m))

print(L_sdf[1])

list(x.GetPropNames())

with open('./SDF/PCCID_00000003.sdf', 'rb') as f:
    fsuppl=Chem.ForwardSDMolSupplier(f,removeHs=False)
    fmols=[]
    for x in fsuppl:
        print(x)

list(x.GetPropNames())

Suppl

Num=48
L_FP[Num]

FP=L_FP[Num]['RP']
suppl=Chem.SDMolSupplier(FP,removeHs=False)
mols = [x for x in suppl if x is not None]
len(mols) # 73

Chem.MolToSmiles(Chem.RemoveHs(mols[0]))

Draw.MolToImage(Chem.MolFromSmiles('NC(CS)C(=O)O'))

Chem.MolToSmiles(Chem.MolFromSmiles('NC(CS)C(=O)O'))

print(L_sdf[Num])

Draw.MolToImage(mols[0])

# +
MB=Chem.MolToMolBlock(mols[0])

view = p3d.view(width=300, height=250, linked=False)
view.addModel(MB,'sdf')
view.setStyle({'stick': {'linewidth': 5}}, viewer=(0,0))
view.zoomTo()
view.show()
