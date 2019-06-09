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

from rdkit import rdBase
from rdkit import Chem
import pubchempy as pcp
print('rdkit version: {}'.format(rdBase.rdkitVersion))

suppl = Chem.SDMolSupplier('../サリチル酸メチル.sdf',removeHs=False)
mols=[]
for x in suppl:
    if x is not None:
        mols.append(x)
len(mols)

from rdkit.Chem import MACCSkeys as maccs
from rdkit.Chem import DataStructs
from rdkit.Chem import AllChem

m_fps_1 = maccs.GenMACCSKeys(mols[0])
m_fps_2 = maccs.GenMACCSKeys(mols[1])

m_mfp_1=AllChem.GetMorganFingerprintAsBitVect(mols[0], 2, 2048)
m_mfp_2=AllChem.GetMorganFingerprintAsBitVect(mols[0], 2, 2048)

DataStructs.TanimotoSimilarity(m_fp_1, m_fp_2)

for prop in mols[0].GetPropNames():
    print(prop)

Draw.MolToImage(mols[2],legend=mols[2].GetProp('PRODUCT_NAME'))

len(m_fps_1.ToBinary())*8

len(m_fps_1.ToBitString())

from rdkit.Chem import Descriptors as des

# +
dec_list=des.descList

calc = {}
for i,j in dec_list:
    calc[i] = j
lipinski_list = ['NumHDonors', 'NumHAcceptors', 'MolWt', 'MolLogP']
# Ro5を判定する関数たち
def calc_lipinsk(mol):
    lipinski = {}
    for desc in lipinski_list:
        lipinski[desc] = calc[desc](mol)
    return lipinski

def check_lipinski(dic):
    if dic['MolWt'] <= 500 and dic['MolLogP'] <= 5 and dic['NumHDonors'] <=5 and dic['NumHAcceptors'] <=10:
        return True
    else:
        return False
    
def rule_of_five(mol):
    prop = calc_lipinsk(mol)
    if check_lipinski(prop):
        return mol


# +
lipinski_mols = []
bad_mols = []
for m in mols:
    if rule_of_five(m):
        lipinski_mols.append(m)
    else:
        bad_mols.append(m)

len(lipinski_mols), len(bad_mols)
# -

from rdkit.Chem import Draw

Draw.MolsToGridImage(bad_mols, molsPerRow=3, subImgSize=(300,200),
                     legends=[x.GetProp('PRODUCT_NAME') for x in bad_mols])

bad_mols

calc['MaxEStateIndex'](mols[0])
