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

with open('./サリチル酸メチル.sdf','rb') as f:
    MSsuppl=Chem.ForwardSDMolSupplier(f,removeHs=False)
    MSs=[]
    for x in MSsuppl:
        if x is not None:
            MSs.append(x)

TarDir='./01MethylSalicylate/'

for c in MSs:
    Smile=Chem.MolToSmiles(c)
    pcp.get_compounds(Smile,'smiles')
    MAc=MSdata[0]
    MAc.cid
