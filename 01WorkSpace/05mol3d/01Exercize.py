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
print('rdkit version: {}'.format(rdBase.rdkitVersion))

import pubchempy as pcp

MA_4133_3dc=pcp.get_compounds(4133,'cid',record_type='3d')

MA_4133_3dc

#sdfをダウンロード
pcp.download('SDF','./mol_4133.sdf',4133,'cid',record_type='3d',overwrite=True)

with open('./4133.sdf', 'rb') as f:
    SDF_4133=f.read()

print(SDF_4133.decode())

import py3Dmol as p3d

view = p3d.view(width=300, height=250, linked=False)
view.addModel(SDF_4133.decode(),'sdf')
view.setStyle({'sphere': {'linewidth': 5}}, viewer=(0,0))
view.zoomTo()
view.show()
