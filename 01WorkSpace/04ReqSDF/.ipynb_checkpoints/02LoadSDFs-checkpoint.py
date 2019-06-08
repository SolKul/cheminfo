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

L_FP=[]
for f in os.listdir(TarDir):
    D_sdf={}
    if f.endswith('.sdf'):
        D_sdf['FN']=f
        D_sdf['RP']=os.path.join(TarDir,f)
        L_FP.append(D_sdf)

FP='./SDF/PCCID_00000003.sdf'
suppl=Chem.SDMolSupplier(FP,removeHs=False)
mols = [x for x in suppl if x is not None]
len(mols) # 73

Draw.MolToImage(mols[0])

MB=Chem.MolToMolBlock(mols[0])

view = p3d.view(width=300, height=250, linked=False)
view.addModel(MB,'sdf')
view.setStyle({'sphere': {'linewidth': 5}}, viewer=(0,0))
view.zoomTo()
view.show()
