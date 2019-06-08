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

from rdkit import rdBase
print('rdkit version: {}'.format(rdBase.rdkitVersion))

# +
import pubchempy as pcp
from rdkit import Chem
 
benzene = pcp.get_compounds('benzene', 'name')
if len(benzene) == 1: benzene = benzene[0]
smiles = benzene.canonical_smiles
print(smiles) # 'C1=CC=CC=C1'
mol_ben = Chem.MolFromSmiles(smiles)
print(type(mol_ben)) # rdkit.Chem.rdchem.Mol
print(Chem.MolToMolBlock(mol_ben))
# -

suppl = Chem.SDMolSupplier('./サリチル酸メチル.sdf',removeHs=False)
mols=[]
for x in suppl:
    if x is not None:
        mols.append(x)
len(mols)

from rdkit.Chem import AllChem, Draw

print(Chem.MolToMolBlock(mols[0]))

Draw.MolsToGridImage(mols[:73], molsPerRow=3, subImgSize=(300,200),
                     legends=[x.GetProp('PRODUCT_NUMBER') for x in mols[:73]])

pcp.get_compounds('COC(=O)c1ccccc1O','smiles',record_type='3d')

suppl

with open('./サリチル酸メチル.sdf', 'rb') as f:
    print(f)

Chem.MolToSmiles(mols[63])

MA_4133_3dc=pcp.get_compounds(4133,'cid',record_type='3d')

MA_4133_3d=MA_4133_3dc[0]
MA_4133_3d.mmff94_energy_3d

type(MA_4133_3d)

pcp.download('SDF','./4133.sdf',4133,'cid',record_type='3d',overwrite=True)

with open('./4133.sdf', 'rb') as f:
    SDF_4133=f.read()

print(SDF_4133.decode())

with open('./4133.sdf', 'rb') as f:
    fsuppl=Chem.ForwardSDMolSupplier(f,removeHs=False)
    fmols=[]
    for x in fsuppl:
        if x is not None:
            fmols.append(x)
    

type(Chem.MolToMolBlock(fmols[0]))

print(Chem.MolToMolBlock(fmols[0]))

import py3Dmol as p3d

view = p3d.view(width=680, height=250, query='cid:5462328', viewergrid=(1,3), linked=False)
view.setStyle({'line': {'linewidth': 5}}, viewer=(0,0))
view.setStyle({'stick':{}}, viewer=(0,1))
view.setStyle({'sphere': {}}, viewer=(0,2))
view.setBackgroundColor('#ebf4fb', viewer=(0,0))
view.setBackgroundColor('#f9f4fb', viewer=(0,1))
view.setBackgroundColor('#e1e1e1', viewer=(0,2))
view.show()

view = p3d.view(width=300, height=250, linked=False)
view.addModel(Chem.MolToMolBlock(fmols[0]),'sdf')
view.setStyle({'sphere': {'linewidth': 5}}, viewer=(0,0))
view.zoomTo()
view.show()

with open('./サリチル酸メチル.sdf','rb') as f:
    MSsuppl=Chem.ForwardSDMolSupplier(f,removeHs=False)
    MSs=[]
    for x in MSsuppl:
        if x is not None:
            MSs.append(x)

for c in MSs:
    for p in c.GetPropNames():
         print('{}: {}'.format(p, c.GetProp(p)))
    print()

for c in MSs:
    print(Chem.MolToSmiles(c))

Chem.MolToSmiles(MSs[4])

MSdata=pcp.get_compounds(Chem.MolToSmiles(MSs[4]),'smiles')

MAc=MSdata[0]
MAc.cid

TarDir='./01MethylSalicylate/'

print('CID{}.sdf'.format(MAc.cid))
SDF_str='CID{}.sdf'.format(MAc.cid)

TarDir+SDF_str

MScid=MSdata[0].cid

pcp.download('SDF',TarDir+SDF_str,MScid,'cid',record_type='3d',overwrite=True)

with open(TarDir+SDF_str,'rb') as f:
    MAsdf=f.read().decode()

view = p3d.view(width=300, height=250, linked=False)
view.addModel(MAsdf,'sdf')
view.setStyle({'sphere': {'linewidth': 5}}, viewer=(0,0))
view.zoomTo()
view.show()
