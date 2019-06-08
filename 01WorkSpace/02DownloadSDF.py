# -*- coding: utf-8 -*-
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

with open('./サリチル酸メチル.sdf','rb') as f:
    MSsuppl=Chem.ForwardSDMolSupplier(f,removeHs=False)
    MSs=[]
    for x in MSsuppl:
        if x is not None:
            MSs.append(x)

TarDir='./01MethylSalicylate/'

N_lim=5
for c in MSs[:N_lim]:
    Smile=Chem.MolToSmiles(c)
    L_comp=pcp.get_compounds(Smile,'smiles')
    comp=L_comp[0]
    comp_cid=comp.cid
    SDF_str='CID{}.sdf'.format(comp_cid)
    pcp.download('SDF',TarDir+SDF_str,comp_cid,'cid',record_type='3d',overwrite=True)

L_FP=[]
for f in os.listdir(TarDir):
    D_sdf={}
    if f.endswith('.sdf'):
        D_sdf['FN']=f
        D_sdf['RP']=os.path.join(TarDir,f)
        L_FP.append(D_sdf)

L_FP

L_sdf=[]
for D_sdf in L_FP:
    with open(D_sdf['RP'],'rb') as f:
        L_sdf.append(f.read().decode())

view = p3d.view(width=300, height=250*N_lim, linked=False,viewergrid=(N_lim,1))
for i in range(N_lim):    
    view.addModel(L_sdf[i],'sdf', viewer=(i,0))
    view.setStyle({'sphere': {'linewidth': 5}}, viewer=(i,0))
    view.zoomTo()
view.show()

print(L_sdf[0])

# + {"language": "bash"}
# pwd
