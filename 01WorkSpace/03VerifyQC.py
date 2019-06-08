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

pcp.download('SDF','./03Verify/PubChem4133.sdf',4133,'cid',record_type='3d',overwrite=True)

suppl=Chem.SDMolSupplier('./03Verify/PubChem4133.sdf')
mols = [x for x in suppl if x is not None]
len(mols) # 73

print(Chem.MolToMolBlock(mols[0]))

Draw.MolToImage(mols[0])

Chem.MolToMolFile(mols[0],'./03Verify/PubChem4133.mol')

Chem.MolFromMolFile('./03Verify/PubChem4133.mol')

QCmol=Chem.MolFromMolFile('./03Verify/QC4133.mol')
print(Chem.MolToMolBlock(QCmol))

Chem.MolFromMolFile('./03Verify/pccdb_id86.mol')

PCC1smils='NCC(=O)COP(=O)(O)O'
L_comp=pcp.get_compounds(PCC1smils,'smiles')

TarDir='./03Verify/'
C_FN='CID{}'.format(C.cid)

FP=TarDir+C_FN+'.sdf'

C=L_comp[0]
FP=TarDir+C_FN+'.sdf'
print(C.cid)
pcp.download('SDF',FP,C.cid,record_type='3d',overwrite=True)

suppl=Chem.SDMolSupplier(FP,removeHs=False)
mols = [x for x in suppl if x is not None]
len(mols) # 73

# +
# PCC1mol=Chem.AddHs(mols[0])
# -

Draw.MolToImage(mols[0])

FP=TarDir+C_FN+'.mol'
Chem.MolToMolFile(PCC1mol,FP)

Chem.MolFromMolFile(FP)

PCCPmols=Chem.MolFromMolFile('./03Verify/pccdb_id5.mol',removeHs=False)
print(Chem.MolToMolBlock(PCCPmols))
