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

import pandas as pd
from rdkit import rdBase, Chem, DataStructs
print(rdBase.rdkitVersion) # 2017.09.1
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import AllChem,Draw
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.AtomPairs import Pairs, Torsions

suppl = Chem.SDMolSupplier('../サリチル酸メチル.sdf',removeHs=False)
mols=[]
for x in suppl:
    if x is not None:
        mols.append(x)
len(mols)

# +
# 
# -

maccs_fps = [AllChem.GetMACCSKeysFingerprint(mol) for mol in mols]
maccs = DataStructs.BulkTanimotoSimilarity(maccs_fps[0], maccs_fps[1:])
